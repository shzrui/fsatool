module CoarseGrain
  use util
  use GlobalVariable, only: nstate, lumpingMethod
  use math
  implicit none
  integer, allocatable :: clusterMapToMacroIndex(:)
  real*8, allocatable :: eigenVector(:,:)
  real*8, allocatable :: membership(:, :)
  real*8, allocatable :: macroTPM(:, :)
contains

! Improved coarse-graining of Markov state models via explicit consideration of statistical uncertainty
! Gregory R. Bowman, J. Chem. Phys. 137, 134111 (2012);

  subroutine CoarseGrainAnalysis(tcmortpm, nstate)
    integer :: nstate, reduced_nstate
    integer :: i
    integer :: nmicro
    real*8, intent(in) :: tcmortpm(:, :)
    nmicro = size(tcmortpm, 1)

    if(lumpingMethod == "bace") then
       call LogInfo("BACE Coarse Grain")
       call CoarseGrainBACE(tcmortpm, nstate)
    elseif(lumpingMethod == "pcca") then
       call LogInfo("PCCA Coarse Grain")
       call CoarseGrainPCCA(tcmortpm, nstate)
    elseif (lumpingMethod == "pccaplus") then
       call LogInfo("PCCAplus Coarse Grain")
       call CoarseGrainPCCAplus(tcmortpm, nstate)
#ifdef DEBUG
       write(*, "(a)") "membership functions:"
       do i = 1, nmicro; print"(40f10.5)", membership(i, :);  enddo
#endif
       call CoarseGrainPCCAPlusTPM(tcmortpm, membership, nstate)
       write(*, "(a)") "PCCA Coarse Grain TPM:"
       do i = 1, nstate; print"(40f10.5)", macroTPM(i, :);  enddo
    else
       call ErrorMessage("lumpingMethod must be bace, pcca, or pccaplus")
    endif

    ! remap the number clusterMapToMacroIndex to 1~nstate
    i = 1; reduced_nstate = nstate;
    do while(.true.)
       if(i > reduced_nstate) exit
       if (.not. any(clusterMapToMacroIndex == i)) then
          where(clusterMapToMacroIndex>i) clusterMapToMacroIndex = clusterMapToMacroIndex-1
          reduced_nstate = reduced_nstate - 1;
       endif
       i = i +1
    enddo

    write(*, "(A)") "cluster index belong to macro states index:"
    write(*, "(20I4)") clusterMapToMacroIndex

    if (nstate > reduced_nstate)  then
       write(*, "(A, I4, 1x, A)") "After lumping only have", reduced_nstate, "macro states."
       write(*, "(A, I4, A, I3, A)")"Lumping algorithm can't coarse grain microstates to", &
       & nstate, " states, but has ", reduced_nstate, " states"
    endif
    nstate = reduced_nstate
    call LogInfo()
  end subroutine CoarseGrainAnalysis

  subroutine ModCoarseGrain(inputFile, resultFile)
   use GlobalVariable, only: ncluster, nstate, lumpingMethod
   real*8, allocatable :: tpm(:, :), macroTPM(:, :), pi(:)
   character(256) :: inputFile, resultFile, datafile
   integer :: ioFile, i, ierr
   namelist/lumping/ datafile, ncluster, nstate , lumpingMethod

   call GetFreeUnit(ioFile)
   open(unit=ioFile, file=trim(inputFile), action="read")
   read(ioFile, nml=lumping, iostat=ierr)
   if(ierr<0) call ErrorMessage("error in reading lumping namelist")
   close(ioFile)

   allocate(tpm(ncluster, ncluster), macroTPM(nstate, nstate), pi(ncluster))

   call GetFreeUnit(ioFile)
   open(unit=ioFile, file=trim(datafile), action="read")
   ! read(ioFile, *) pi
   do i = 1, ncluster
      read(ioFile, *) tpm(i, :)
   enddo
   close(ioFile)
   call CoarseGrainAnalysis(tpm, nstate)

   macroTPM = 0d0
   call GetFreeUnit(ioFile)

   ! if(lumpingMethod == "bace") then
   !    do i=1,ncluster
   !       do j=1,ncluster
   !          macroTPM(clusterMapToMacroIndex(i),clusterMapToMacroIndex(j)) = macroTPM(clusterMapToMacroIndex(i), clusterMapToMacroIndex(j)) &
   !            + tpm(i,j)
   !       end do
   !    enddo
   !    forall(i=1:nstate)
   !       macroTPM(i, :) = macroTPM(i,:)/sum(macroTPM(i, :))
   !    end forall
   !    open(unit=ioFile, file=trim(resultFile), action="write")
   !    do i = 1, nstate
   !       write(ioFile, *) macroTPM(i, :)
   !    enddo
   ! endif
   open(unit=ioFile, file=trim(resultFile), action="write")
   do i = 1, ncluster
      write(ioFile, *) clusterMapToMacroIndex(i)
   enddo
   close(ioFile)
  
  end subroutine

  subroutine CoarseGrainPCCAPlusTPM(tpm, membership, nstate)
    real*8, intent(in) :: tpm(:, :) , membership(:, :)
    integer, intent(in) :: nstate
    real*8 :: inverse_matrix(nstate, nstate)
    call MathMatrixInverse(matmul(transpose(membership), membership), inverse_matrix)
    macroTPM=matmul(inverse_matrix, matmul(matmul(transpose(membership), tpm), membership))
  end subroutine CoarseGrainPCCAPlusTPM

  subroutine CoarseGrainBACE(tcm, nstate)
    real*8, intent(in) :: tcm(:, :)
    integer, intent(in) :: nstate
    integer :: i,j,istate,jstate,iter,imin,jmin,stateleft, tcmDim, minxminy(2), helparray(nstate)
    real*8,dimension(:), allocatable :: tpqvec, sumci
    real*8 :: spiq,spjq
    real*8,dimension(:,:),allocatable :: tptcm, tpm, bayesmat
    integer,allocatable :: allowedcalc(:,:), statesold(:)

    tcmDim = size(tcm,1)
    allocate(tpqvec(tcmDim), clusterMapToMacroIndex(tcmDim), sumci(tcmDim), statesold(tcmDim))
    allocate(bayesmat(tcmDim, tcmDim), tpm(tcmDim, tcmDim), tptcm(tcmDim, tcmDim),allowedcalc(tcmDim, tcmDim))
    tptcm=tcm; bayesmat=0.0d0; stateleft=tcmDim; imin=0

    allowedcalc = 0
    where (tcm>1) allowedcalc = 1

    do i=1,tcmDim
       statesold(i) = i
       clusterMapToMacroIndex(i)=i
       tptcm(i, 1:tcmDim) = tptcm(i, 1:tcmDim) + 1.0d0/dble(tcmDim)
       sumci(i) = sum(tptcm(i,1:tcmDim))
       tpm(i,1:tcmDim) = tptcm(i,1:tcmDim)/sumci(i)
    end do

    do iter=tcmDim,nstate+1,-1
       if (iter==tcmDim) then
          do istate=1,stateleft
             do jstate=istate+1, stateleft
                if (allowedcalc(istate, jstate) > 0) then
                   tpqvec(1:stateleft) = (tptcm(istate,1:stateleft)+tptcm(jstate,1:stateleft))/(sumci(istate)+sumci(jstate))
                   spiq = dot_product(tptcm(istate, 1:stateleft), log(tpm(istate,1:stateleft)/tpqvec(1:stateleft)))
                   spjq = dot_product(tptcm(jstate, 1:stateleft), log(tpm(jstate,1:stateleft)/tpqvec(1:stateleft)))
                   bayesmat(istate,jstate) = 1.0d0/(spiq + spjq)
                endif
             enddo
          enddo
       else
          do istate=1,stateleft
             if (allowedcalc(imin, istate)>0) then
                tpqvec(1:stateleft) = (tptcm(imin,1:stateleft)+tptcm(istate,1:stateleft))/(sumci(imin)+sumci(istate))
                spiq = dot_product(tptcm(imin, 1:stateleft), log(tpm(imin,1:stateleft)/tpqvec(1:stateleft)))
                spjq = dot_product(tptcm(istate, 1:stateleft), log(tpm(istate,1:stateleft)/tpqvec(1:stateleft)))
                bayesmat(imin,istate) = 1.0/(spiq + spjq)
             endif
          enddo
       endif

       minxminy = maxloc(bayesmat(1:stateleft, 1:stateleft));
       if (minxminy(1) < minxminy(2)) then
          imin = minxminy(1); jmin=minxminy(2)
       else
          imin = minxminy(2); jmin=minxminy(1)
       endif

       tptcm(imin, 1:stateleft) = tptcm(imin, 1:stateleft) + tptcm(jmin, 1:stateleft)
       tptcm(1:stateleft, imin) = tptcm(1:stateleft, imin) + tptcm(1:stateleft, jmin)
       bayesmat(imin, 1:stateleft) = 0
       bayesmat(1:stateleft, imin) = 0

       where (clusterMapToMacroIndex == clusterMapToMacroIndex(statesold(jmin)))
          clusterMapToMacroIndex = clusterMapToMacroIndex(statesold(imin))
       end where

       do jstate=jmin,stateleft-1
          tptcm(jstate,1:stateleft) = tptcm(jstate+1,1:stateleft)
          tptcm(1:stateleft,jstate) = tptcm(1:stateleft,jstate+1)
          bayesmat(jstate,1:stateleft) = bayesmat(jstate+1,1:stateleft)
          bayesmat(1:stateleft,jstate) = bayesmat(1:stateleft,jstate+1)
       end do

       do jstate = jmin+1, stateleft
          statesold(jstate-1) = statesold(jstate)
       enddo

       stateleft = stateleft - 1
       do jstate=1,stateleft
          sumci(jstate) = sum(tptcm(jstate,1:stateleft))
          tpm(jstate,1:stateleft) = tptcm(jstate,1:stateleft)/sumci(jstate)
       end do

       allowedcalc(imin, :) = 0
       do jstate=1, stateleft
          if(tptcm(imin, jstate)>1 .and. jstate /= imin) then
             allowedcalc(imin, jstate) = 1
          endif
       enddo
    enddo
    ! renumber the clusterMapToMacroIndex to 1:n
    helparray=0; j=1
    do i =1, tcmDim
       if(any(helparray == clusterMapToMacroIndex(i))) then
          clusterMapToMacroIndex(i) = minloc(helparray, dim=1, mask=(helparray>=clusterMapToMacroIndex(i)))
       else
          helparray(j) = clusterMapToMacroIndex(i)
          clusterMapToMacroIndex(i) = j
          j = j + 1
       endif
    enddo
    deallocate(tpqvec, sumci, statesold, bayesmat, tpm, tptcm, allowedcalc)
  end subroutine CoarseGrainBACE

  subroutine CoarseGrainPCCA(tpm, nstate)
    ! CoarseGrain by pcca method
    real*8, dimension(:,:), intent(in) :: tpm
    integer, intent(in) :: nstate
    real*8, allocatable :: temparray(:), temparray2(:)
    integer :: nmicro,i,j,maxpos

    nmicro = size(tpm, 1)
    if (allocated(clusterMapToMacroIndex)) deallocate(clusterMapToMacroIndex)
    if (allocated(eigenVector)) deallocate(eigenVector)
    allocate(temparray(nmicro), eigenVector(nmicro, nmicro), temparray2(nmicro), clusterMapToMacroIndex(nmicro))
    call MathSolveEigenProblem(tpm, temparray, eigenVector, temparray2)
    clusterMapToMacroIndex=1
    do i = 1, nstate - 1
       temparray =  eigenVector(:, i+1)
       do j = 1, i
          temparray2(j) = maxval(pack(temparray, clusterMapToMacroIndex==j))-minval(pack(temparray, clusterMapToMacroIndex==j))
       enddo
       maxpos = maxloc(temparray2(1:j-1), dim=1)
       where((clusterMapToMacroIndex == maxpos) .and. (temparray>0.0d0))
          clusterMapToMacroIndex = i+1
       end where
    enddo
    deallocate(temparray, eigenVector, temparray2)
  end subroutine CoarseGrainPCCA

  subroutine CoarseGrainPCCAplus(tpm, nstate)
    real*8, intent(in) :: tpm(:, :)
    integer, intent(in) :: nstate
    integer :: nmicro, i, j, icount, index, nparameter, numres, ifault
    real*8 :: temp, dist2
    integer, allocatable :: stateIndex(:)
    real*8, allocatable :: a(:, :), ainv(:, :)
    real*8, allocatable :: pi(:), eigenValue(:), temparray(:)
    real*8, allocatable :: reducedeigvec(:,:), step(:), mincoor(:), xinit(:)

    nmicro = size(tpm, 1)
    allocate(pi(nmicro), eigenValue(nmicro), eigenVector(nmicro, nmicro),&
         stateIndex(nstate), temparray(nstate), a(nstate, nstate), ainv(nstate, nstate), &
         reducedeigvec(nmicro, nstate))
    call MathSolveEigenProblem(tpm, eigenValue, eigenVector, pi)

    !normalize the eigenValue by the weight pi
    do i = 1, nmicro
       temp = dot_product(eigenVector(:, i)**2, pi)
       eigenVector(:, i) = eigenVector(:, i) / sqrt(temp)
    enddo

    eigenVector(:, 1) = abs(eigenVector(:, 1))
    reducedeigvec = eigenVector(:, 1:nstate)

    ! Need to find the first index, select the most spread index
    ! https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
    stateIndex = 0; dist2 = 0
    do i = 1, nmicro
       temp = dot_product(reducedeigvec(i,:), reducedeigvec(i, :))
       if (temp > dist2) then
          dist2 = temp
          stateIndex(1) = i
       endif
    enddo
    index = stateIndex(1)
    temparray = reducedeigvec(index, :)

    do i=1, nmicro
       reducedeigvec(i, :) = reducedeigvec(i, :) - temparray
    enddo

    do i = 1, nstate - 1
       dist2 = 0.0d0
       temparray = reducedeigvec(stateIndex(i), :)
       do j = 1, nmicro
          reducedeigvec(j,:) = reducedeigvec(j, :) - dot_product(reducedeigvec(j, :), temparray)*temparray
          temp = dot_product(reducedeigvec(j, :), reducedeigvec(j, :))
          if ((temp > dist2) .and. (.not. any(stateIndex == j))) then
             dist2  = temp
             stateIndex(i+1) = j
          endif
       enddo
       temp = sqrt(dot_product(reducedeigvec(stateIndex(i+1), :), reducedeigvec(stateIndex(i+1), :)))
       reducedeigvec = reducedeigvec / temp
    enddo

    a =  eigenVector(stateIndex, 1:nstate)
    call MathMatrixInverse(a, ainv) ! get the inverse matrix of a
    nparameter = (nstate-1)*(nstate-1)
    allocate(mincoor(nparameter), step(nparameter), xinit(nparameter))
    xinit = pack(transpose(ainv(2:nstate, 2:nstate)), .true.)

    step = 0.05
    ! call Neled-Mead method to find mincoor

    membership = matmul(eigenVector(:, 1:nstate), ainv)
    call MathOptimazationNelderMean(CoarseGrainPCCAPlusObjectFunction, nparameter, xinit, mincoor, temp, & 
                           0.00005d0, step, 10, 2000000, icount, numres, ifault)
    print*, "call Neled-Mead method to minimize the func"
    if (ifault == 2) then
       print*, "iteration terminated because kcount was exceed without convergence, maybe increase max steps"
    else
       print*, "icount, numres:", icount, numres
    endif

    ! rescale the mincoor to ainv
    ainv(2:nstate, 2:nstate) = reshape(mincoor, (/nstate-1, nstate-1/), order=[2,1])
    do i = 2, nstate
       ainv(i, 1) = -sum(ainv(i,2:nstate))
    enddo
    ainv(1, :) = -minval(matmul(eigenVector(:, 2:nstate), ainv(2:nstate, :)), dim=1)
    ainv = ainv / sum(ainv(1,:))

    ! At last, we can get membership function, normalize the membership function
    membership = matmul(eigenVector(:, 1:nstate), ainv)
    where (membership<0) membership = 0.0d0
    where (membership>1) membership = 1.0d0
    do i = 1, nmicro
       membership(i, :) = membership(i, :) / sum(membership(i, :))
    enddo

    ! get the clusterMapToMacroIndex function
    clusterMapToMacroIndex = maxloc(membership, dim=2)
    deallocate(pi, eigenValue, eigenVector, stateIndex, &
         temparray, a, ainv, reducedeigvec)
  end subroutine CoarseGrainPCCAplus

  real*8 function CoarseGrainPCCAPlusObjectFunction(ainv, nparameter)
    integer :: nparameter
    real*8,intent(in) :: ainv(nparameter)
    real*8,allocatable :: expandainv(:, :)
    integer :: n, i, j

    i = int(sqrt(dble(nparameter)))
    n = i + 1
    allocate(expandainv(n, n))
    expandainv(2:n, 2:n) = reshape(ainv, (/i, i/), order=[2,1])
    do i = 2, n
       expandainv(i, 1) = -sum(expandainv(i,2:n))
    enddo
    expandainv(1, :) = -minval(matmul(eigenVector(:, 2:n), expandainv(2:n, :)), dim=1)
    expandainv = expandainv / sum(expandainv(1,:))
    CoarseGrainPCCAPlusObjectFunction = 0.0d0
    do i= 1, n
       do j=1, n
          CoarseGrainPCCAPlusObjectFunction = CoarseGrainPCCAPlusObjectFunction + (expandainv(j,i)**2) / expandainv(1,i)
       enddo
    enddo
    CoarseGrainPCCAPlusObjectFunction = -CoarseGrainPCCAPlusObjectFunction
    deallocate(expandainv)
    return
  end function CoarseGrainPCCAPlusObjectFunction

end module CoarseGrain
