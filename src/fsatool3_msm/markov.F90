module MARKOV
  use util
  use GlobalVariable
  use cluster
  use CoarseGrain, only: CoarseGrainAnalysis, clusterMapToMacroIndex
  use math, only: MathSolveEigenProblem
  implicit none
  integer :: tcmDim, lagOffset
  integer,allocatable :: reducedTCMIndex(:)
  real*8,dimension(:,:),allocatable :: tpm,tcm,eigenVec,macroTPM, impliedTimeScale
  real*8,dimension(:),allocatable :: limitpi,eigenVal,macropi
contains


  subroutine MarkovCheck(inputFile, resultFile)
    use mpi
    character(256) :: inputFile, resultFile, datafile, checkmethod
    integer :: i, ioFile, ierr, lagstart, lagend, nits
    namelist /check/ checkmethod, datafile, lagstart, lagend, nits, lagOffset

    nits = 5; lagOffset = 1.0d0
    call GetFreeUnit(ioFile)
    open(unit=ioFile, file=trim(inputFile), action="read")
    read(ioFile, nml=check, iostat=ierr)
    if(ierr < 0) call ErrorMessage("error in reading check namelist")
    close(ioFile)

    call GetFreeUnit(ioFile)
    call CheckFileAndOpen(datafile, ioFile)
    read(ioFile, *) cvFrameNumber, ncluster
    allocate(maptocluster(cvFrameNumber), trajindex(cvFrameNumber), ifTrimCluster(cvFrameNumber))
    ifTrimCluster = 0
    do i = 1, cvFrameNumber
      read(ioFile, *) maptocluster(i), trajindex(i)
    enddo
    call CheckFileAndClose(ioFile)

    if (trim(checkmethod) == "timescales") then
      call MarkovimpliedTimeScaleVsLagTime(resultFile, nits, lagstart, lagend)
    ! else if (trim(checkmethod) == "cktest")  then
    !   call makrov_task2_cktest()
    endif
  end subroutine

  subroutine MarkovAnalysis()
    use tram, only:TramAnalysis
    integer :: i, j, ioFile
    integer,allocatable :: trimmed_cluster(:)
    allocate(ifTrimCluster(ncluster));  ifTrimCluster = 0

    if(procId > 0) return
    if (cutNumber /= 0) then
      do i = 1, ncluster
          if (numberincluster(i) < cutNumber) ifTrimCluster(i) = 1
      enddo
    endif

    if(nlevel == 1) then
      call LogInfo("Build MSM from clusters")
      write(*, "(A, I6)") "Markov lag step:", lagStep
      call MarkovBuildTCMAndTPMFromClusters()
      allocate(eigenVal(tcmDim), eigenVec(tcmDim, tcmDim), limitpi(tcmDim))
      call MathSolveEigenProblem(tpm, eigenVal, eigenVec, limitpi)

      ! save reduced tcm matrix
      call GetFreeUnit(ioFile); open(unit=ioFile, file=trim(resultDir)//"/tcm_reduce.txt", action="write")
      do i = 1, tcmDim
        write(ioFile, "(10F6.0)") tcm(i, :)
      enddo
      close(ioFile)

      ! save tpm matrix
      call GetFreeUnit(ioFile); open(unit=ioFile, file=trim(resultDir)//"/tpm.txt", action="write")
      write(ioFile, "(10E12.5)") limitpi
      do i = 1, tcmDim
        write(ioFile, "(10E12.5)") tpm(i, 1:tcmDim)
      enddo
      close(ioFile)
    else
      call TramAnalysis(tpm, reducedTCMIndex)
      tcmDim = size(reducedTCMIndex, 1)
    endif

    if((ncluster-tcmDim) > 0) then
      allocate(trimmed_cluster(ncluster-tcmDim))
      j = 1
      do i = 1, ncluster
          if (.not. any(i == reducedTCMIndex(1:tcmDim))) then
            trimmed_cluster(j) = i
            j = j+1
          endif
      enddo
    endif
    write(*, "(a, I8, a, I5, a)") "System has ", cvFrameNumber, " microstates, corresponding  to", ncluster, " cluster number"
    if(cutNumber /= 0) then
      write(*,  "(a, I4)") "remove the cluster that has number of states which less than", cutNumber
    endif
    write(*, "(a)")"build tpm using the "//trim(tpmMethod)// "  method."
    write(*, "(a,i4)")"new dimension of TCM after reducing ",tcmDim
    write(*, "(a)") "cluster index which has been trimmed:"
    write(*, "(20I4)") trimmed_cluster
    call LogInfo()

    if (lumpingMethod == "bace") then
       call CoarseGrainAnalysis(tcm, nstate)
    else
       call CoarseGrainAnalysis(tpm, nstate)
    end if

    ! call LogInfo("Build MSM from Marco States")
    ! call MarkovBuildCoarseGrainTPM()
    ! write(*,"(a, I3, a, I5, a)") "The system has", nstate, " states, correspond to", tcmDim, " clusters."
    ! write(*, "(a)") "Marco tpm:"
    ! do i = 1,nstate; write(*, "(100f10.6)") macroTPM(i,:); enddo
    ! write(*, "(a)") "Pi of marco states:"
    ! write(*, "(10f8.3)")macropi

    ! call GetFreeUnit(ioFile); open(unit=ioFile, file=trim(resultDir)//"/macroinfo", action="write")
    !   write(ioFile, "(100ES15.8)") macropi
    !   write(ioFile, *)
    !   do i = 1,nstate; write(ioFile, "(100ES15.8)") macroTPM(i,:); enddo
    ! close(ioFile)
  end subroutine MarkovAnalysis

  subroutine MarkovBuildTCMAndTPMFromClusters()
    integer :: i,j,icluster,jcluster,isnap, nums_scc
    integer :: scc_label_array(ncluster)
    logical :: scc_filter(ncluster)
    real*8 :: temp,tpsum,tptcm(ncluster,ncluster)

    if ( allocated(tcm) ) deallocate(tcm);
    allocate(tcm(ncluster,ncluster)); tcm=0.0d0; tcmDim = 0

    do isnap=1,cvFrameNumber-lagStep
      i=maptocluster(isnap); j=maptocluster(isnap+lagStep)
      if ( (ifTrimCluster(i)==0) .and. (ifTrimCluster(j)==0) .and. (trajindex(isnap)==trajindex(isnap+lagStep))) then
          tcm(i,j)=tcm(i,j)+1.0d0
      end if
    end do

    tptcm = tcm
    call find_scc(tptcm, ncluster, scc_label_array, nums_scc)
    ! scc_filter is the logical array to keep the states at the largest scc
    do i = 1, ncluster
      if(scc_label_array(i) == 1) tcmDim = tcmDim + 1
    enddo

    if(allocated(tcm)) deallocate(tcm)
    if (allocated(reducedTCMIndex)) deallocate(reducedTCMIndex)
    if (allocated(tpm)) deallocate(tpm)
    allocate(reducedTCMIndex(tcmDim), tcm(tcmDim, tcmDim), tpm(tcmDim, tcmDim))

    j = 1
    scc_filter = .true.
    do i = 1, ncluster
      if(scc_label_array(i) == 1) then
        reducedTCMIndex(j) = i
        j = j+1
      else
        scc_filter(i) = .false.
      endif
    enddo

    j = 1;
    do i =1, ncluster
      if(scc_filter(i)) then
        tcm(j, :) = pack(tptcm(i, :), scc_filter)
        j = j+1
      endif
    enddo

    tpm=0d0
    select case (trim(tpmMethod))
    case("normal")
      tpsum = sum(tcm)
      do icluster=1,tcmDim
          temp = sum(tcm(icluster,:)); tpm(icluster,:)=tcm(icluster,:)/temp
      end do
    case("non-reversible")
      do icluster=1,tcmDim
          do jcluster=icluster+1,tcmDim
            temp = (tcm(icluster,jcluster)+tcm(jcluster,icluster))*0.5d0
            tcm(icluster,jcluster)=temp; tcm(jcluster,icluster)=temp
          end do
      end do
      tpsum = sum(tcm)
      do icluster=1,tcmDim
          temp = sum(tcm(icluster,:)); tpm(icluster,:)=tcm(icluster,:)/temp
          !limitpi(icluster) = temp/tpsum
      end do
    case("reversible")
      call MarkovSymmetrizeTCMAndNormalizeTPM()
    case default
      call MarkovSymmetrizeTCMAndNormalizeTPM()
   end select
  end subroutine MarkovBuildTCMAndTPMFromClusters

  subroutine MarkovBuildCoarseGrainTPM()
    integer :: i,icluster,jcluster, index

    if ( allocated(macroTPM) ) deallocate(macroTPM)
    if ( allocated(macropi) ) deallocate(macropi)
    allocate(macroTPM(nstate,nstate),macropi(nstate)); macroTPM=0.0d0; macropi=0.0d0
    do icluster=1,tcmDim
      do jcluster=1,tcmDim
          macroTPM(clusterMapToMacroIndex(icluster),clusterMapToMacroIndex(jcluster)) = &
          macroTPM(clusterMapToMacroIndex(icluster), clusterMapToMacroIndex(jcluster)) + tcm(icluster,jcluster)
       end do
    enddo
    forall(i=1:nstate)
      macroTPM(i, :) = macroTPM(i,:)/sum(macroTPM(i, :))
    end forall
    do i = 1, tcmDim
      index = clusterMapToMacroIndex(i)
      macropi(index) = numberincluster(reducedTCMIndex(i)) + macropi(index)
    enddo
    macropi = macropi/sum(macropi)
  end subroutine MarkovBuildCoarseGrainTPM


  subroutine MarkovTimeScale(tplagTime, tpimpliedTime, nlagTime)
    real*8,intent(in) :: tplagTime
    integer, intent(in) :: nlagTime
    real*8, intent(inout) :: tpimpliedTime(:)
    allocate(eigenVal(tcmDim), eigenVec(tcmDim, tcmDim), limitpi(tcmDim))
    call MathSolveEigenProblem(tpm, eigenVal, eigenVec, limitpi)
    tpimpliedTime(2:nlagTime) = (-tplagTime)/log(abs(eigenVal(2:nlagTime)))
    deallocate(eigenVal, eigenVec, limitpi)
  end subroutine MarkovTimeScale

  !Book:An Introduction to Markov State Models and Their Application to Long Timescale Molecular Simulation
  ! Chapter 4.6.1, page 50

  subroutine MarkovSymmetrizeTCMAndNormalizeTPM()
    integer :: i,j,icycle
    real*8 :: temp,xivec(tcmDim),civec(tcmDim),qlog,lastqlog
    real*8,dimension(:,:),allocatable :: tpmatrix

    allocate(tpmatrix(tcmDim,tcmDim)); tpmatrix=0.0d0
    do i=1,tcmDim
      do j=i,tcmDim
          call random_number(temp)
          tpmatrix(i,j)=temp
          tpmatrix(j,i)=temp
      end do
    end do

    do i=1,tcmDim
      civec(i)=sum(tcm(i,:))
    end do

    do
      icycle = icycle + 1
      do i=1,tcmDim;  xivec(i)=sum(tpmatrix(i,:)); end do
      qlog=0.0d0
      do i=1,tcmDim
          do j=i,tcmDim
            temp = (civec(i)/xivec(i)) + (civec(j)/xivec(j))
            tpmatrix(i,j)=(tcm(i,j)+tcm(j,i))/temp; tpmatrix(j,i)=tpmatrix(i,j)
            if ( tpmatrix(i,j) > 0.0d0 ) qlog = qlog + tcm(i,j)*log(tpmatrix(i,j)/xivec(i))
            if ( (j > i) .and. (tpmatrix(j,i) > 0.0d0) ) qlog = qlog + tcm(j,i)*log(tpmatrix(j,i)/xivec(j))
          end do
      end do
      if ( abs(qlog-lastqlog) < 1d-10) exit
      lastqlog=qlog
    end do
    do i=1,tcmDim
      xivec(i)=sum(tpmatrix(i,:)); tpm(i,:) = tpmatrix(i,:)/xivec(i)
    end do
    deallocate(tpmatrix)
  end subroutine

  subroutine MarkovimpliedTimeScaleVsLagTime(resultFile, nits, lagStart, lagEnd)
    character(MaxStringLen) :: resultFile
    integer :: i, j, ioImplied, nits, lagStart, lagEnd, lagTemp, ierr, lagBins
    integer, allocatable :: lagArray(:), subLagArray(:)
    real*8, allocatable :: subimpliedTimeScale(:, :)
    integer :: left, right, length
    integer :: displs(procNum), counts(procNum)

    if(procId == 0) then
      call LogInfo("Task1: calculate the impliedTimescale on lagTime")
      lagBins = (lagEnd - lagStart) / lagOffset
      allocate(lagArray(lagBins))
      do i = 1, lagBins
        lagArray(i) = lagStart + i * lagOffset
      enddo

      allocate(impliedTimeScale(nits, lagBins))
      call GetFreeUnit(ioImplied); open(unit=ioImplied,file=trim(resultFile),action="write")
      
      write(*, "(a, I7, a)") "Need to calculate ", lagbins, " timesteps"
      write(*, "(a, 10I5)") "lagStep:", lagArray
      write(ioImplied, *) "#lagStep #implied_timescale(nits-1)"
    endif

    call mpi_bcast(lagBins, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call PartitionCPUProcessors(procNum, lagBins, displs, counts, left, right)

    length = right - left + 1 

    allocate(subLagArray(length))
    allocate(subimpliedTimeScale(nits, length))
    call mpi_scatterv(lagArray, counts, displs, mpi_int, subLagArray, length, mpi_int, 0, mpi_comm_world, ierr)
    do i  = 1, size(subLagArray, 1)
      lagStep = dble(subLagArray(i))
      call MarkovBuildTCMAndTPMFromClusters()
      call MarkovTimeScale(dble(lagStep), subimpliedTimeScale(:, i), nits)
      write(*, "(a, I6, a)") "finished calculating", lagStep, " lagStep"
    enddo

    call mpi_gatherv(subimpliedTimeScale, length*nits, mpi_double,& 
      impliedTimeScale, counts*nits, displs*nits, mpi_double, 0, mpi_comm_world, ierr)

    if(procId == 0) then
      do i = 1, lagBins
        lagStep = lagArray(i)
        write(ioImplied,"(i4, 1x, 20f12.6)") lagStep, log(impliedTimeScale(2:nits, i))
      enddo
      close(ioImplied)
      call LogInfo()
    deallocate(impliedTimeScale, lagArray)
    endif
    deallocate(subLagArray, subimpliedTimeScale)
  endsubroutine MarkovimpliedTimeScaleVsLagTime
end module MARKOV
