module tpt
  use util
  use GlobalVariable, only: procId, nstate, startState, endState, resultDir, MaxStringLen
  use fileio, only: stateIndex, clusterMapToMacroIndex => states_to_cluster
  use math, only: MathGetTPMPi, MathMatrixInverse
  use markov, only: tpm, pi=>limitpi, tcmDim
  implicit none
  integer :: numberPathWay
  integer, parameter :: maxStateNumber = 10
  integer, allocatable :: Astate(:), Bstate(:)
  integer, allocatable :: strongestPathways(:)
  integer, pointer :: sendStatePointer(:), recvStatePointer(:)
  integer :: numStatesInA, numStatesInB
  real*8 :: totalFlux,kAB
  real*8,dimension(:,:),allocatable :: mfpt, netFlux, flux, macroNetFlux
  real*8,dimension(:),allocatable :: forwardCommittor, backwardCommittor
  type :: pathwayinfo
     real*8 :: flux
     integer :: nodes
     integer, allocatable :: pathway(:)
     real*8 :: fluxRatio
  end type pathwayinfo
  TYPE(pathwayinfo), allocatable, dimension(:) :: pathways
contains

  subroutine ModTPT(inputFile, outputfile)
    character(MaxStringLen) :: inputFile, outputfile
    if(procId > 0) return
    startState = 0; endState = 0
    call TPTModReadFile(inputFile)
    call TPTAssignStates()
    allocate (flux(size(tpm, 1), size(tpm, 1)), macroNetFlux(nstate, nstate))
    resultDir = "./"
    call TPTCalculatePathway(outputfile)
  end subroutine ModTPT

  subroutine TPTAnalysis()
    integer :: i
    if(procId>0) return
    stateIndex(1) = 1
    do i = 2, nstate + 1
      stateIndex(i) = stateIndex(i) + 1
    enddo
    call TPTAssignStates()
    allocate(flux(size(tpm, 1), size(tpm, 1)), macroNetFlux(nstate, nstate))
    call TPTCalculatePathway()
  end subroutine

  subroutine TPTModReadFile(inputFile)
    character(MaxStringLen) :: inputFile, tpmFile, statesFile
    logical :: readPi, calculateMFPT
    integer :: i, j, prevIntegerTemp, ioFile, ierr
    integer :: temp 
    namelist /tpt/ tpmFile, statesFile, startState, endState, nstate, readPi, calculateMFPT

    statesFile = ""
    readPi = .true.
    calculateMFPT = .true.

    if (procId  > 0) return
    call GetFreeUnit(ioFile)
    open(unit=ioFile, file=trim(inputFile), action="read")
    read(ioFile, nml=tpt,iostat=ierr)
    if (ierr<0) call ErrorMessage("error in reading the tpt namelist")
    close(ioFile)

    if (statesFile == "") then
      tcmDim = nstate
    else
      call GetFreeUnit(ioFile)
      open(unit=ioFile, file=trim(statesFile), action="read")
      i = 0
      do
        i = i + 1
        read(ioFile, *, iostat=ierr)
        if (ierr < 0) exit
      enddo
      tcmDim = i - 1
      close(ioFile)
    endif

    allocate(tpm(tcmDim, tcmDim), pi(tcmDim), stateIndex(nstate+1), clusterMapToMacroIndex(tcmDim))
    if(allocated(netFlux)) deallocate(netFlux)
    allocate(netFlux(tcmDim, tcmDim))

    call GetFreeUnit(ioFile)
    open(unit=ioFile, file=trim(tpmFile), action="read")

    if (readPi) read(ioFile, *) pi
    do i = 1, tcmDim
      read(ioFile, *) tpm(i, :)
    enddo
    close(ioFile)

    if (.not. readPi) call MathGetTPMPi(tpm, tcmDim, pi)

    if (statesFile /= "") then
      call GetFreeUnit(ioFile)
      open(unit=ioFile, file=trim(statesFile), action="read")
      prevIntegerTemp = 1; j = 1; stateIndex(j) = 1
      do i = 1, tcmDim
        read(ioFile, *) temp, clusterMapToMacroIndex(i)
        if (temp /= prevIntegerTemp) then 
          prevIntegerTemp = temp
          j = j + 1
          stateIndex(j) = i
        endif
      enddo
      stateIndex(j+1) = tcmDim+1
      close(ioFile)
    else
      stateIndex(1) = 1
      do i = 1, tcmDim
        clusterMapToMacroIndex(i) = i
        stateIndex(i+1) = i+1
      enddo
    endif
  end subroutine

  subroutine TPTAssignStates()
    integer :: i, j, k, count
    count = 0
    ! count number of cluster in Astate
    do i = 1,  maxStateNumber
      if (startState(i) /= 0) then
        count = count + stateIndex(startState(i)+1) - stateIndex(startState(i))
      else
        exit
      endif
    enddo
    allocate(Astate(count))
    sendStatePointer => startState(1:i-1)
    numStatesInA = i - 1

    !count number of cluster in Bstate
    count = 0
    do i = 1, maxStateNumber
      if (endState(i) /= 0) then
        count = count + stateIndex(endState(i)+1) - stateIndex(endState(i))
      else
        exit
      endif
    enddo
    allocate(Bstate(count))
    recvStatePointer => endState(1:i-1)
    numStatesInB = i - 1

    ! Given the c
    j = 1
    do k = 1, maxStateNumber
      if(startState(k) /= 0) then
        do i = stateIndex(startState(k)), stateIndex(startState(k)+1) - 1
          Astate(j) = clusterMapToMacroIndex(i)
          j = j + 1
        enddo
      endif
    enddo

    j = 1
    do k = 1, maxStateNumber
      if(endState(k) /= 0) then
        do i = stateIndex(endState(k)), stateIndex(endState(k)+1)-1
          Bstate(j) = clusterMapToMacroIndex(i)
          j = j + 1
        enddo
      endif
    enddo
  end subroutine

  subroutine TPTCalculatePathway(outputfile)
    character(MaxStringLen), optional :: outputfile
    integer :: i, ioFile
    character(MaxStringLen) :: resultFile = "macroNetFlux.txt"
    
    if (present(outputfile)) resultFile = outputfile
    call LogInfo("TPT Analysis")
    call TPTGetFlux(tpm, size(tpm, 1), Astate, Bstate, flux, reversible=.false., pi=pi)
    call TPTCoarseGrainFlux(flux, macroNetFlux)
    call TPTExtractPathWays(macroNetFlux, numberPathWay, sendStatePointer, recvStatePointer, totalFlux)
    call TPTPrintResult(numberPathWay)
    call GetFreeUnit(ioFile)

    open(unit=ioFile, file= trim(resultDir) //"/" // resultFile, action="write")
    do i = 1, size(macroNetFlux ,1)
      write(ioFile, "(100f12.9)") macroNetFlux(i, :)
    enddo
    close(ioFile)
    call LogInfo()
  end subroutine

  subroutine TPTPrintResult(numberPathWay)
    integer :: numberPathWay 
    integer :: i
    real*8 :: accumulatedfluxRatio
    real*8,allocatable :: tempDoubleArray(:)
    integer, allocatable :: sortIndex(:)

    write(*, "(a, I8)") "number of micro state: ", tcmDim
    write(*, "(a, I5)")"number of macro state: ", nstate
    write(*, "(a, 10I4)") "int state: ", startState(1:numStatesInA)
    write(*, "(a, 10I4)")"end state: ", endState(1:numStatesInB)
    write(*, "(a, 50f12.8)")"forward committor: ",forwardCommittor(:)
    write(*, "(a, 500f12.8)")"backward committor: ",backwardCommittor(:)
    write(*, "(a, 50f10.6)")"pi: ", pi
    write(*, "(a, f12.8)")"totalFlux: ", totalFlux
    write(*, "(a, f12.8)")"kAB: ", kAB
    write(*, "(a, I3, a)")"There are total ", numberPathWay, " pathways in the network"
    allocate(tempDoubleArray(numberPathWay), sortIndex(numberPathWay))
    do i = 1, numberPathWay
        tempDoubleArray(i) = pathways(i) % fluxRatio
        sortIndex(i) = i
    enddo
   ! numpath is not so big, so using simple sort to return sortIndex
    accumulatedfluxRatio = 0d0
    call TPTGetFeatureIndex(tempDoubleArray, numberPathWay, sortIndex)
    write(*, "(a, 5x, a, 2x, a, 2x, a)") "flux", "ratio", "accumulatedfluxRatio", "pathway"
    do i = 1, numberPathWay
      accumulatedfluxRatio = accumulatedfluxRatio + pathways(sortIndex(i))%fluxRatio
      write(*, "(E12.4, 2f8.4, 3x, 20I3)") pathways(sortIndex(i))%flux, pathways(sortIndex(i)) %fluxRatio, &
      & accumulatedfluxRatio, pathways(sortIndex(i))%pathway
    enddo
  end subroutine

  subroutine TPTGetFlux(p, numberOfState, fromState, toState, flux, reversible, pi)
    integer,intent(in) :: numberOfState
    real*8, intent(in):: p(numberOfState, numberOfState)
    real*8, intent(out) :: flux(numberOfState, numberOfState)
    logical, optional :: reversible
    real*8, dimension(:), optional :: pi
    integer, dimension(:) :: fromState, toState
    real*8 :: lhs(numberOfState, numberOfState), rhs(numberOfState)
    real*8 :: subpi(numberOfState)
    integer :: i, j, ipiv(numberOfState), info

    lhs=p; rhs = 0.0d0
    ! calculate forward committor
    forall(i=1:numberOfState) lhs(i,i) = lhs(i,i) - 1.0d0
    do i = 1, size(fromState, 1)
      lhs(fromState(i), :)=0.0d0
      lhs(:, fromState(i)) = 0.0d0
      lhs(fromState(i), fromState(i)) = 1.0d0
    enddo
    do i =1, size(toState, 1)
      lhs(toState(i),:) = 0.0d0
      lhs(toState(i), toState(i)) = 1.0d0
      rhs(toState(i)) = 1.0d0
    enddo

    call DGESV(numberOfState, 1, lhs, numberOfState, ipiv, rhs, numberOfState, info)
    forwardCommittor = rhs
    if (present(reversible) .and. reversible) then
      backwardCommittor = 1.0d0 - forwardCommittor
      subpi = pi
    else
      ! calculate backward committor
      ! call MathGetTPMPi(p, numberOfState, subpi)
      rhs = 0.0d0
      do i = 1, numberOfState
         do j = 1, numberOfState
            lhs(i,j) = pi(j)* p(j,i) / pi(i)
         enddo
      enddo
      forall(i=1:numberOfState) lhs(i,i) = lhs(i,i) - 1.0d0
      do i = 1, size(fromState, 1)
         lhs(fromState(i), :)=0.0d0
         lhs(fromState(i), fromState(i)) = 1.0d0
         rhs(fromState(i)) = 1.0d0
      enddo
      do i = 1, size(toState, 1)
         lhs(toState(i),:) = 0.0d0
         lhs(:, toState(i)) = 0.0d0
         lhs(toState(i), toState(i)) = 1.0d0
      enddo
      call DGESV(numberOfState, 1, lhs, numberOfState, ipiv, rhs, numberOfState, info)
      backwardCommittor = rhs
    endif

    flux = 0.0d0
    forall(i=1:numberOfState, j=1:numberOfState, i/=j) flux(i,j)=pi(i)*p(i,j)*backwardCommittor(i)*forwardCommittor(j)
  end subroutine 

  subroutine TPTExtractPathWays(netFlux, numpath, fromState, toState, totalFlux)
    real*8 :: minFlux, totalFlux
    real*8,dimension(:, :) :: netFlux
    integer :: numpath
    integer, dimension(:), pointer :: fromState, toState
    real*8, allocatable :: netFluxcopy(:, :)
    integer :: i, nodes, index, j, k

   numpath = 0; allocate(pathways(100))
    do k = 1, size(fromState, 1)
      do j = 1, size(toState, 1)
         netFluxcopy = netFlux
         do
            call dijkstra(netFluxcopy, size(netFluxcopy, 1), fromState(k), toState(j))
            nodes = size(strongestPathways, 1)
            if(nodes == 1) exit
            numpath = numpath + 1
            minFlux = 1e10
            do i = 1, nodes-1
               if (minFlux > netFluxcopy(strongestPathways(i), strongestPathways(i+1)))  then
                  minFlux = netFluxcopy(strongestPathways(i), strongestPathways(i+1))
                  index = i
               endif
            enddo
            do i = 1, nodes-1
               netFluxcopy(strongestPathways(i), strongestPathways(i+1)) = &
                     netFluxcopy(strongestPathways(i), strongestPathways(i+1)) - &
                     minFlux
            enddo
            netFluxcopy(strongestPathways(index), strongestPathways(index+1)) = 0.0d0
            pathways(numpath)%flux = minFlux
            pathways(numpath)%nodes = nodes
            pathways(numpath)%pathway = strongestPathways
            pathways(numpath)%fluxRatio =  minFlux / totalFlux
         enddo
      enddo
   enddo
   if (numpath > 100) write(*,"(A)") "numpath is larger than allocated number of pathway, allocate pathway larger"
  end subroutine TPTExtractPathWays

  subroutine dijkstra(network, numberOfState, fromState, toState)
    integer :: numberOfState, fromState, toState
    real*8,intent(in) :: network(numberOfState, numberOfState)
    real*8 ::  dist(numberOfState), temp
    integer :: visitedState(numberOfState), previousState(numberOfState)
    integer :: j, newState, oldState, tempArray(numberOfState)

    visitedState = 0; dist = -huge(1.0); oldState = toState
    newState = fromState; visitedState(newState) = 1
    do
       if(newState == toState) exit ! if find end, then exit
       if(newState == oldState) then ! if can't find other node next to intistate
          if (allocated(strongestPathways)) deallocate(strongestPathways)
          allocate(strongestPathways(1))
          return
       endif
       oldState = newState
       do  j = 1, numberOfState
          if(network(oldState, j) > 0d0 .and. visitedState(j) == 0) then
             if(network(oldState, j) > dist(j))  then
                dist(j) = network(oldState, j)
                previousState(j) = oldState
             endif
          endif
       enddo
       temp = 0
       do j = 1, numberOfState
          if(visitedState(j) == 0) then
             if(dist(j) > temp) then
                temp = dist(j)
                newState = j
             endif
          endif
       enddo
       visitedState(newState) = 1
    enddo
    newState = toState; j = 1; tempArray(j) = toState
    do
       if(newState == fromState) exit
       newState = previousState(newState)
       j = j + 1
       tempArray(j) = newState
    enddo
    if (allocated(strongestPathways)) deallocate(strongestPathways)
    allocate(strongestPathways(j)); strongestPathways = tempArray(j:1:-1)
  end subroutine dijkstra


  subroutine TPTCoarseGrainFlux(flux, macroNetFlux)
    real*8, intent(in), dimension(:, :) :: flux
    real*8, intent(out), dimension(nstate, nstate) ::  macroNetFlux
    real*8, DIMENSION(nstate, nstate) :: macroFlux 
    real*8, dimension(nstate) :: pFold, macroStationaryDistribution, qMinus
    integer :: i, j, m, n
    real*8 :: tempFlux, tempDoubleSum, temp, tempDoubleSum2

    macroFlux = 0.d0
    do i = 1, nstate
      do j = 1, nstate
        tempFlux = 0d0
        if(i/=j) then
          do m = stateIndex(i), stateIndex(i+1) - 1
            do n = stateIndex(j), stateIndex(j+1) - 1
              tempFlux = tempFlux + flux(clusterMapToMacroIndex(m), clusterMapToMacroIndex(n))
            enddo
          enddo
          macroFlux(i, j) = tempFlux
        endif
      enddo
    enddo

    ! forall(i=1:nstate, j=1:nstate, i/=j) 
    forall(i=1:nstate, j=1:nstate) 
      macroNetFlux(i,j)=max(macroFlux(i,j)-macroFlux(j,i), 0.0d0)
    endforall

    totalFlux = 0d0
    do i = 1, maxStateNumber
      if(endState(i)/= 0) then
        totalFlux = sum(macroNetFlux(:, endState(i)))
      else
        exit
      endif
    enddo

    do i = 1, nstate
      temp = 0d0; tempDoubleSum=0d0; tempDoubleSum2 = 0d0
      do j = stateIndex(i), stateIndex(i+1) - 1
        temp = temp + pi(clusterMapToMacroIndex(j))
      enddo
      macroStationaryDistribution(i) = temp
      do j = stateIndex(i), stateIndex(i + 1) - 1
        tempDoubleSum = tempDoubleSum + forwardCommittor(clusterMapToMacroIndex(j)) * pi(clusterMapToMacroIndex(j))
        tempDoubleSum2 = tempDoubleSum2 + backwardCommittor(clusterMapToMacroIndex(j)) * pi(clusterMapToMacroIndex(j))
      enddo
      pFold(i) = tempDoubleSum/temp
      qMinus(i) = tempDoubleSum2/temp
    enddo

    if (allocated(forwardCommittor)) then
      DEALLOCATE(forwardCommittor, pi, backwardCommittor)
      allocate(forwardCommittor(nstate), pi(nstate), backwardCommittor(nstate))
      forwardCommittor = pFold
      pi = macroStationaryDistribution
      backwardCommittor = qMinus
    endif
    kAB = totalFlux / sum(pi * backwardCommittor)
  end subroutine

  subroutine TPTGetFeatureIndex(tempArray, features, sortIndex)
    integer :: features
    real*8 :: temp, tempArray(features)
    integer :: tempint, sortIndex(features)
    integer :: i, j

    do i = 1, features
      sortIndex(i) = i
    enddo
    do i = 1, features
      do j = i+1, features
          if (tempArray(j) > tempArray(i)) then
            temp = tempArray(i)
            tempArray(i) = tempArray(j)
            tempArray(j)  = temp

            tempint = sortIndex(i)
            sortIndex(i) = sortIndex(j)
            sortIndex(j) = tempint
          endif
      enddo
    enddo
  end subroutine

  subroutine TPTMfpt(macroTPM)
    real*8, intent(in) :: macroTPM(:, :)
    real*8 :: fundamentalMatrix(nstate, nstate)
    real*8 :: tempDoubleMatrix(nstate, nstate)
    integer :: i, j

    allocate(mfpt(nstate, nstate))
    if (allocated(pi)) deallocate(pi)
    allocate(pi(nstate))
    call MathGetTPMPi(tpm, nstate, pi)
    tempDoubleMatrix = - macroTPM
    forall(i=1:nstate)
       tempDoubleMatrix(i, i) = 1 + tempDoubleMatrix(i, i)
       tempDoubleMatrix(i, :) = tempDoubleMatrix(i, :) + pi(:)
    end forall
    call MathMatrixInverse(tempDoubleMatrix, fundamentalMatrix)
    mfpt = 0
    forall(i=1:nstate, j=1:nstate, i/=j)
       mfpt(i, j) = (fundamentalMatrix(j,j) - fundamentalMatrix(i, j)) / pi(j)
    end forall
    deallocate(pi)
  end subroutine TPTMfpt

end module
