module GlobalVariable
  use util
  implicit none
  integer, parameter :: MaxStringLen = 256
  integer, parameter :: MaxLevel=20
  integer, parameter :: MaxCVNumber = 100
  integer, parameter :: MaxClusterAtomNumber = 1000
  real*8,parameter :: GasConstant=1.9872065d-3
  real*8 :: kelvinArray(MaxLevel)
  integer :: numFilePerLevel(MaxLevel), clusterAtomIndex(MaxClusterAtomNumber)
  integer :: atomDegree, cvFrameNumber, procId, procNum, cvDegree, trajType, nlevel, ncluster, clusterCycle
  integer :: cutNumber, lagStep, nstate, randomSeed
  integer, dimension(:), allocatable :: levelIndex, ifTrimCluster, trajIndex
  real*8, dimension(:,:), pointer :: traj, cvs
  real*8, dimension(:), allocatable :: potentialPerSnap
  real*8 :: lagTime,kelvin
  logical :: ifReadCoordinate, ifOnlyRunClustering
  character(MaxStringLen),save :: resultDir, tpmMethod, clusterMethod, lumpingMethod, inputFile
  integer, dimension(MaxLevel), target :: startState, endState
contains

subroutine GlobalVariableMPIInitialize()
  use mpi
  integer :: ierr

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,procId,ierr)
  call mpi_comm_size(mpi_comm_world,procNum,ierr)
end subroutine

subroutine GlobalVariableReadTrajs()
  use mpi
  use MPISharedMemory
  use netcdf_func
  integer, parameter :: maxnfile=30
  character(MaxStringLen) :: cvFiles(maxnfile),coorFiles(maxnfile)
  integer :: numFile
  integer, allocatable :: numFramePerFile(:)
  integer :: netcdfFrameNumber, cvStride, natom, i
  integer :: ioFile, ierr, shared_cv_win, shared_traj_win
  namelist /trajs/ kelvinArray, trajType, cvFiles, coorFiles, numFilePerLevel, clusterAtomIndex

  cvFiles=""; cvStride = 1; kelvinArray=0; numFilePerLevel=0; clusterAtomIndex = 0; trajType=0

  if ( procId == 0 ) then
    call LogInfo("Read Collective Variables File")

    call GetFreeUnit(ioFile)
    open(unit=ioFile, file=inputFile, action="read")
    read(ioFile, nml=trajs, iostat=ierr)
    if ( ierr < 0 ) call ErrorMessage("error in reading trajs namelist")
    close(ioFile)

    call GetFileNumber(cvFiles, numFile)
    allocate(numFramePerFile(numFile))
    call GetCvFileInfo(cvFiles, numFile, cvFrameNumber, cvDegree, 1, numFilePerLevel)
    cvDegree = cvDegree - 1 ! the first column is trajIndex

    write(*, "(a, I8)")     "numer of files           = ", numFile
    write(*, "(a, I8)")     "number of cv             = ", cvDegree
    write(*, "(a, I8)")     "number of CV Frame       = ", cvFrameNumber

    if (trajType == 2) then 
      call GetNetCdfInfo(coorFiles, numFile, netcdfFrameNumber, atomDegree)
      cvStride = cvFrameNumber / netcdfFrameNumber
      cvFrameNumber = netcdfFrameNumber
    endif
  endif

  call mpi_bcast(trajType, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(numFile, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(cvFiles, numFile*MaxStringLen, mpi_char, 0, mpi_comm_world, ierr)
  call mpi_bcast(cvDegree, 1, mpi_int, 0, mpi_comm_world, ierr)
  call mpi_bcast(cvStride, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(cvFrameNumber, 1, mpi_integer, 0, mpi_comm_world, ierr)
  
  if (trajType /= 0) then 
    call mpi_bcast(coorFiles, numFile*MaxStringLen, mpi_char, 0, mpi_comm_world, ierr)
    call mpi_bcast(atomDegree, 1, mpi_integer, 0, mpi_comm_world, ierr)
  endif

  allocate(trajIndex(cvFrameNumber), potentialPerSnap(cvFrameNumber), levelIndex(cvFrameNumber))

  call InitalizeSharedMemory(mpi_comm_world)
  call MPISharedMemoryAllocate2DimDoubleArray([cvDegree, cvFrameNumber], cvs, shared_cv_win)

  if (trajType /= 0) then
    call MPISharedMemoryAllocate2DimDoubleArray([atomDegree, netcdfFrameNumber], traj, shared_traj_win)
  endif

  if(procId == 0) then
    call GetLevelFromKelvinArray()
    write(*, "(a, I8)")     "dimension of traj        = ", atomDegree
    write(*, "(a, I8)")     "Trajectory type          = ", trajType
    write(*, "(a, I8)")     "number of NetCDF Frame   = ", cvFrameNumber
    write(*, "(a, I8)")     "Frame stride             = ", cvStride
    write(*, "(a, I8)")     "Number of temperature    = ", nlevel
    write(*, "(a, 10f8.3)") "Kelvin array             = ", kelvinArray(1:nlevel)
  endif

  if(shared_id == 0) then
    call ReadCvFile(cvFiles, numFile, cvs, cvStride, trajIndex)
    if(trajType == 1) then
      call ReadAmberTrajFile(coorFiles, numFile, atomDegree, traj, cvStride)
    else if(trajType == 2) then
      call NetcdfReadCoordinate(coorFiles, numFile, atomDegree, traj, 1)
    endif
  endif
  call mpi_win_fence(0, shared_cv_win, ierr)
  if(trajType /= 0) call mpi_win_fence(0, shared_traj_win, ierr)
  if(procId==0) call LogInfo()
end subroutine

subroutine GetFileNumber(files, numFile)
  character(MaxStringLen), intent(in) :: files(:)
  integer, intent(out) :: numFile
  integer :: i

  i = 1
  do 
    if(files(i) == "" ) exit
    i = i + 1
  enddo
  numFile = i - 1
end subroutine

subroutine GetCvFrameNumber(files, numFile, stride, frameLine, numFramePerFile)
  character(MaxStringLen), intent(in) :: files(:)
  integer, intent(in) :: numFile, stride, frameLine
  integer, intent(out), optional :: numFramePerFile(:)
  integer :: i, j, ierr, ioFile

  cvFrameNumber = 0
  do i = 1, numFile
    call checkFileAndOpen(files(i), ioFile)
    outer: do
      inner: do j = 1, stride * frameLine 
        read(ioFile, *, iostat=ierr)
        if (ierr < 0) exit outer
      enddo inner
      cvFrameNumber = cvFrameNumber + 1
      if(present(numFramePerFile)) numFramePerFile(i) = numFramePerFile(i) + 1
    enddo outer
    call CheckFileAndClose(ioFile)
  enddo
end subroutine

subroutine GetNetcdfInfo(files, numFile, numFrame, numDimension)
  use NETCDF_FUNC
  character(MaxStringLen), intent(in) :: files(:)
  integer, intent(in):: numFile
  integer, intent(out) :: numFrame, numDimension

  logical :: ifExist
  integer :: i, ierr, ioFile, temp, numAtom

  numFrame = 0
  do i = 1, numFile
    inquire(file=files(i), exist=ifExist)
    if(.not. ifExist) call ErrorMessage("File " // trim(files(i)) // " does not exist")
    call NetcdfInitParameter(files(i), temp, numAtom)
    numFrame = numFrame + temp
  enddo
  numDimension = numAtom * 3
end subroutine

subroutine ReadCvFile(files, numFile, cvs, cvStride, trajIndex)
  character(MaxStringLen), intent(in) :: files(:)
  integer, intent(in) :: numFile, cvStride
  integer, intent(out) :: trajIndex(:)
  real*8, intent(out) :: cvs(:, :)

  real*8 ::  cv_temp(cvDegree)
  integer :: snapIndex = 1
  integer :: frameLine, ioFile
  integer :: numFrames(numFile)
  integer :: i, j, k, ierr, ilevel, accumFileNum

  ilevel = 1; accumFileNum = 0;
  inner: do i = 1, numFile
    if(i - accumFileNum > numFilePerLevel(ilevel)) then 
      accumFileNum = accumFileNum + numFilePerLevel(ilevel)
      ilevel = ilevel + 1
    endif
    call GetFreeUnit(ioFile)
    open(unit = ioFile, file=files(i), action="read")
    outer: do 
      do  j = 1, (cvStride - 1) ! skip using cvStride
        read(ioFile, *, iostat=ierr)
        if(ierr < 0) exit outer
      enddo
      if (nlevel == 1) then ! we don't need read potential
        read(ioFile, *, iostat=ierr) trajIndex(snapIndex), cvs(:, snapIndex)
      else ! need to specify the potential 
        write(*, "(a)") "PLEASE MAKE SURE THAT CVSFILE SHOULD CONTAIN POTENTIAL AFTER TRAJINDEX"
        read(ioFile, *, iostat=ierr) trajIndex(snapIndex), potentialPerSnap(snapIndex), cvs(:, snapIndex)
      endif
      levelIndex(snapIndex) = ilevel
      if (ierr < 0) exit
      snapIndex = snapIndex + 1
    enddo outer
  enddo inner
end subroutine

subroutine ReadAmberTrajFile(files, numFile, atomDegree, trajs, cvStride)
  character(MaxStringLen), intent(in) :: files(:)
  integer, intent(in) :: numFile, atomDegree, cvStride
  real*8, intent(out) :: trajs(:, :)

  integer :: snapIndex = 0
  integer :: frameLine, ioFile
  integer :: numFrames(numFile)
  integer :: i, j

  frameLine = ceiling( dble(atomDegree) / 10)

  do i = 1, numFile
    call GetFreeUnit(ioFile)
    open(unit = ioFile, file=files(i), action="read")
    read(ioFile, *)
    do  j = 1, frameLine * (cvStride - 1) ! skip using cvStride
      read(ioFile, *)
    enddo
    snapIndex = snapIndex + 1
    if(snapIndex > cvFrameNumber) exit
    read(ioFile, "(10f8.3)") trajs(:, snapIndex)
  enddo
end subroutine

subroutine GetLevelFromKelvinArray()
  integer :: i
  nlevel = 0
  do i = 1, MaxLevel
    if(kelvinArray(i) > 0) nlevel = nlevel + 1
  enddo
end subroutine

subroutine GetCvFileInfo(files, numFile, cvNum, cvDimension, stride, numFilePerLevel) 
  character(MaxStringLen), intent(in) :: files(:)
  integer, intent(in) :: numFile
  integer, intent(in) :: stride
  integer, optional, intent(out) :: numFilePerLevel(:)
  integer, intent(out) :: cvNum, cvDimension

  real*8 :: tempArray(MaxCVNumber)
  character(MaxStringLen) :: tempStr
  integer :: i, ierr

  integer :: ioFile
  call CheckFileAndOpen(files(1), ioFile)
  read(ioFile, '(A)') tempStr
  do i = 1, MaxCVNumber
    read(tempStr, *, iostat = ierr) tempArray(1:i)
    if (ierr < 0 ) exit
  enddo
  call CheckFileAndClose(ioFile)
  cvDimension = i - 1
  if (present(numFilePerLevel))  then
    call getCvFrameNumber(files, numFile, stride, 1, numFilePerLevel)
  else 
    call getCvFrameNumber(files, numFile, stride, 1, numFilePerLevel)
  endif
  cvNum = cvFrameNumber
end subroutine

subroutine CheckFileAndOpen(file, ioFile)
  character(MaxStringLen), INTENT(IN) :: file
  integer, intent(out) :: ioFile

  logical :: ifExist
  inquire(file=file, exist=ifExist)
  if(.not. ifExist) call ErrorMessage("File " // trim(file) // " does not exist")
  call GetFreeUnit(ioFile)
  open(file=file, unit=ioFile, action="read")
  write(*, '(A)') "Open and Read " // trim(file) 
end subroutine

subroutine CheckFileAndClose(ioFile)
  integer, intent(in) :: ioFile
  logical :: ifUsed  
  inquire(unit = ioFile, opened = ifUsed)
  if(.not. ifUsed) call ErrorMessage("File unit is not used, cannot close it")
  close(ioFile)
end subroutine

subroutine PartitionCPUProcessors(procNum, cvFrameNumber, displs, counts, left, right)
  integer, intent(in) :: procNum, cvFrameNumber
  integer, intent(out) :: displs(:), counts(:)
  integer, intent(out) :: left, right
  integer :: i, length
  do i = 0, procNum-1
    left = i * cvFrameNumber / procNum + 1
    right = (i+1) * cvFrameNumber / procNum
    length = right - left + 1
    displs(i+1) = left - 1
    counts(i+1) = length
  enddo
  left = procId * cvFrameNumber / procNum + 1
  right = (procId+1) * cvFrameNumber / procNum
end subroutine

subroutine RestoreRandomSeed(randomSeed)
  integer, intent(in) :: randomSeed
  integer, allocatable :: seedArray(:)
  integer :: i

  call RANDOM_SEED(size = i)
  allocate(seedArray(i))
  seedArray = randomSeed
  call RANDOM_SEED(put=seedArray)
  deallocate(seedArray)
end subroutine

subroutine GetRMSDAtomNumber(num)
  integer, intent(out) :: num
  integer :: i
  num = atomDegree / 3
  if(clusterAtomIndex(1) /= 0) then
    do i = 1, MaxClusterAtomNumber
      if(clusterAtomIndex(i) == 0) exit
    enddo
    num = i - 1
  endif
end subroutine

end module