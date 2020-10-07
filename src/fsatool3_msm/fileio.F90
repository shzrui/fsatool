module fileio
  use util
  use GlobalVariable
  use cluster, only: maptocluster, centersnaps, numberincluster
  use CoarseGrain, only: clusterMapToMacroIndex
  use markov, only: reducedTCMIndex, tcmDim
  use math, only: MathAtomRMSD
  implicit none
  integer, allocatable :: clusterindex(:), clusters(:), stateIndex(:)
  integer, allocatable :: macrocenter, states_to_cluster(:)
contains

  subroutine fileio_default_parameters()
    use GlobalVariable, only: clusterMethod, clusterCycle, cutNumber, &
                          tpmMethod, lumpingMethod, &
                          nstate, ifReadCoordinate, trajType
    clusterMethod="kmeans"
    clusterCycle=600
    cutNumber=0
    tpmMethod = "reversible"
    lumpingMethod = "pccaplus"
    nstate = 8
    ifReadCoordinate = .true.
    ifOnlyRunClustering = .false.
    trajType = 0
  end subroutine

  subroutine fileio_init_parameters()
    use mpi
    logical :: file_exist
    integer :: ierr, system, ioFile, i, j
    namelist /msm/  ncluster, clusterMethod, clusterCycle,cutNumber, &
         lagStep, ifReadCoordinate, lumpingMethod, tpmMethod, nstate, startState, endState, ifOnlyRunClustering

    call fileio_default_parameters()

    if(procId == 0) call LogInfo("Initialize parameters")
    call GetFreeUnit(ioFile); open(unit=ioFile,file=trim(inputFile),action="read")
    read(ioFile,nml=msm,iostat=ierr)
    if ( ierr < 0 ) call ErrorMessage("error in reading the msm namelist!")
    close(ioFile)

    if ( procId == 0 ) then
       if (ifReadCoordinate .eqv. .true.) then
          ierr = system("rm -rf "//trim(resultDir))
          ierr = system("mkdir -p "//trim(resultDir))
       endif
       write(*, "(a,a)") "result directory: ",trim(resultDir)
       write(*, *) "cluster method: ", trim(clusterMethod)
       write(*, *) "lumping method: ", trim(lumpingMethod)
    end if
    if ( procId == 0 ) call LogInfo()
  end subroutine fileio_init_parameters

  subroutine fileio_writeclusterinfo()
    integer :: iocluster,index,isnap,icluster,i,j,totalcount, atomDegree

    if(procId > 0) return
    allocate(clusterindex(ncluster+1),clusters(cvFrameNumber))
    clusterindex = 0
    totalcount = 0
    do i = 1, ncluster
       do j = 1, cvFrameNumber
          if (maptocluster(j) == i) then
             totalcount = totalcount + 1
             clusters(totalcount) = j ! clusters map cluster index to snap
          endif
       enddo
       clusterindex(i+1) = totalcount
    enddo

    call LogInfo("Writing cluster information")
    write(*, "(a)") "Cluster information has been putted into:" // trim(resultDir)
    write(*, "(a)") "cluster info has two files clusterindex.txt and clusters.txt"

    call GetFreeUnit(iocluster); open(unit=iocluster,file=trim(resultDir)//"/clusterindex.txt",action="write")
    write(iocluster,"(10i10)") ncluster,cvFrameNumber
    do icluster=1,ncluster
       write(iocluster,"(i6,3i10,10f10.5)") icluster,clusterindex(icluster+1),clusterindex(icluster+1)-clusterindex(icluster), &
            centersnaps(icluster), cvs(:, centersnaps(icluster))
    end do
    close(iocluster)
    
    call GetFreeUnit(iocluster); open(unit=iocluster,file=trim(resultDir)//"/cluster_forcheck.txt",action="write")
    write(iocluster, "(2I10)") cvFrameNumber, ncluster
    do i = 1, cvFrameNumber
       write(iocluster,"(i6,3i10,10f10.5)") maptocluster(i), trajindex(i)
    end do
    close(iocluster)

    call GetFreeUnit(iocluster); open(unit=iocluster,file=trim(resultDir)//"/clusters.txt",action="write")
    do icluster=1,ncluster
       do i=clusterindex(icluster)+1,clusterindex(icluster+1)
          isnap = clusters(i)
          write(iocluster,"(i10,2i8,10f10.5)") isnap, icluster, trajindex(isnap), cvs(:, isnap)
       end do
    end do
    close(iocluster)
    call LogInfo()
  end subroutine fileio_writeclusterinfo

  subroutine fileio_readclusterinfo()
    integer :: iocluster,isnap,icluster,i,j,k,ierr

    if(procId > 0) return

    call LogInfo("Reading cluster information")
    print "(a)", "Reading Cluster information from:" // trim(resultDir)
    print "(a)", "Reading two files clusterindex.txt and clusters.txt"

    call GetFreeUnit(iocluster) 
    open(unit=iocluster,file=trim(resultDir)//"/clusterindex.txt",action="read")
    read(iocluster,"(10i10)") ncluster,cvFrameNumber
    print "(5(a,i7,2x))","read cluster file, cluster number: ",ncluster," snap number:",cvFrameNumber
    if ( allocated(clusters) ) deallocate(clusters)
    if ( allocated(clusterindex) ) deallocate(clusterindex)
    if ( allocated(centersnaps) ) deallocate(centersnaps)
    allocate(clusters(cvFrameNumber),clusterindex(ncluster+1),centersnaps(ncluster)); 
    clusters=0; clusterindex=0; centersnaps=0

    if ( allocated(trajindex) ) deallocate(trajindex)
    if (allocated(numberincluster)) deallocate(numberincluster)
    if (allocated(maptocluster)) deallocate(maptocluster)
    allocate(trajindex(cvFrameNumber),cvs(cvDegree,cvFrameNumber),numberincluster(ncluster), maptocluster(cvFrameNumber))

    trajindex=0; cvs=0.0d0
    do icluster=1,ncluster
       read(iocluster,"(i6,3i10)") isnap,clusterindex(icluster+1),numberincluster(icluster),centersnaps(icluster)
    end do
    close(iocluster)

    call GetFreeUnit(iocluster); open(unit=iocluster,file=trim(resultDir)//"/clusters.txt",action="read")
    do icluster=1,ncluster
       do i=clusterindex(icluster)+1,clusterindex(icluster+1)
          read(iocluster,"(i10,2i8,10f10.6)") isnap,maptocluster(isnap),trajindex(isnap), cvs(1:cvDegree,isnap)
          clusters(i) = isnap
       end do
    end do
    close(iocluster)
    call LogInfo()
  end subroutine fileio_readclusterinfo

  subroutine fileio_writestateinfo()
    use math
    integer :: i, j, k, numstate, index, istate, totalcount, isnap, iostate
    integer :: numcluster_state(nstate), temp(tcmDim), state_center_snap(nstate), helparray(tcmDim)
    integer ::  state_center_cluster(nstate), numsnap_state(nstate)
    real*8 :: total, bignum, temp_double

    if(procId > 0) return
    call LogInfo("Writing state information")
    write(*, *) "The state info has been putted into " // trim(resultDir)
    allocate(states_to_cluster(tcmDim))

    totalcount = 0
    do i = 1, nstate
       do j = 1, tcmDim
          if(clusterMapToMacroIndex(j) == i) then
             totalcount = totalcount + 1
             states_to_cluster(totalcount) = j !map same continuous state index to cluster index
          endif
       enddo
       numcluster_state(i) = count(clusterMapToMacroIndex==i)
    enddo

    allocate(stateIndex(nstate+1)); stateIndex=0
    do i = 1, nstate
       stateIndex(i+1) = stateIndex(i) + numcluster_state(i)
    enddo

    numsnap_state = 0
    do i = 1, tcmDim
       numsnap_state(clusterMapToMacroIndex(i)) = numsnap_state(clusterMapToMacroIndex(i)) + numberincluster(reducedTCMIndex(i))
    enddo

    helparray = (/(i, i=1, tcmDim)/)

    ! find state center correspond to cluster center
    do i = 1, nstate
       numstate = numcluster_state(i)
       temp(1:numstate) = centersnaps(reducedTCMIndex(pack(helparray, clusterMapToMacroIndex==i)))
       bignum = 1e10
       do j = 1, numstate
          total = 0.0d0
          do k = 1, numstate
            if(trajType /= 0) then
              call MathAtomRMSD(traj(:, temp(j)), traj(:, temp(k)), temp_double, atomDegree)
            else
              call MathEuclideanDistanceSquare(cvs(:, temp(j)), cvs(:, temp(k)), temp_double, cvDegree)
            endif
             total = total + temp_double
          enddo
          if (total < bignum) then
             bignum = total
             state_center_snap(i) = temp(j)
             state_center_cluster(i) = maptocluster(temp(j))
          endif
       enddo
    enddo

    call GetFreeUnit(iostate); open(unit=iostate,file=trim(resultDir)//"/stateIndex.txt",action="write")
    write(iostate,"(i6,2i10)") nstate, tcmDim, sum(numsnap_state)
    do i = 1, nstate
       write(iostate,"(5i10, i10,10f8.3)") i, stateIndex(i+1),numcluster_state(i), numsnap_state(i), &
            state_center_cluster(i), state_center_snap(i), cvs(1:cvDegree, state_center_snap(i))
    enddo
    close(iostate)

    call GetFreeUnit(iostate); open(unit=iostate,file=trim(resultDir)//"/states.txt",action="write")
    index=0
    do istate=1,nstate
       do i=stateIndex(istate)+1,stateIndex(istate+1)
          index = index + 1
          isnap = centersnaps(reducedTCMIndex(states_to_cluster(i)))
         !  write(iostate,"(3i8, 1x, i10,10f8.3)") istate,states(i),reducedTCMIndex(states(i)), isnap, cvs(1:atomDegree, isnap)
          write(iostate,"(2i8, 1x, i10,10f8.3)") istate,states_to_cluster(i), isnap, cvs(1:cvDegree, isnap)
       end do
    end do
    close(iostate)
    call LogInfo()
  end subroutine fileio_writestateinfo

  subroutine mod_writetraj(inputFile, clusterresultDir)
    use netcdf_func
    use cluster, only: ClusterWriteCoordianteResultToNetcdfFile
    character(MaxStringLen) :: inputFile, clusterresultDir
    character(MaxStringLen) :: coorFiles(30)
    integer :: clusterIndex(100)
    integer :: numFile, ioFile, clusterNumber
    integer :: i, ierr, cvStride, atomNumber
    namelist /writeTraj/ coorFiles, cvStride, clusterIndex, clusterresultDir

    if(procId > 0) return
    !clusterresultDir = clusterDir
    cvStride = 1
    clusterIndex = 0
    coorFiles = ""
    call GetFreeUnit(ioFile)
    open(unit=ioFile, file=trim(inputFile), action="read")
    read(ioFile, nml=writeTraj, iostat=ierr)
    if(ierr<0) call ErrorMessage("error in reading writeTraj namelist")
    close(ioFile)

    call LogInfo("writing clustering result into netcdf file")
    call NetcdfInitParameter(coorFiles(1), tempAtomNumber=atomNumber)
    atomDegree = atomNumber * 3

    i = 1
    do 
       if(clusterIndex(i) == 0) exit
       i = i + 1
    enddo
    clusterNumber = i-1
    write(*, "(A, 10I5)") "writing cluster", clusterIndex(1:clusterNumber)


    call fileio_readclusterinfo() 
    allocate(traj(atomDegree, cvFrameNumber))

    call GetFileNumber(coorFiles, numFile)
    call NetcdfReadCoordinate(coorFiles, numFile, atomDegree, traj, cvStride)
    call ClusterWriteCoordianteResultToNetcdfFile(traj, atomDegree, cvFrameNumber, clusterIndex(1:clusterNumber), clusterresultDir)
    call LogInfo()
  end subroutine

end module fileio
