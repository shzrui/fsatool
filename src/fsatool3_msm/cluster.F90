module cluster
  use math, only: MathAtomRMSD, MathEuclideanDistanceSquare, MathNormalization
  use MPISharedMemory
  use util
  use GlobalVariable
  use netcdf_func
  implicit none
  integer,dimension(:),allocatable :: maptocluster, centersnaps, numberincluster
  integer :: rmsdAtomNumber
contains

  subroutine ClusterAnalysis()
    use mpi
    integer :: stat, i, j, ierr
    real*8, pointer :: normalized_cvs(:, :)
    integer, ALLOCATABLE :: subsetTrajIndex(:)
    integer :: win_normalized_shared

    call MPISharedMemoryAllocate2DimDoubleArray([cvDegree, cvFrameNumber], normalized_cvs, win_normalized_shared)

    if (trajType == 0) then
      if (shared_id == 0) then
        call MathNormalization(cvs, normalized_cvs, cvDegree, cvFrameNumber)
      endif
    endif

    call mpi_win_fence(0, win_normalized_shared, ierr)

    if (trim(clusterMethod) == "kmeans") then
       if(procId == 0) call LogInfo("Running Kmeans++ Cluster Algorithm")
       call ClusterKmeans(normalized_cvs, cvDegree, ncluster, cvFrameNumber, clusterCycle, stat, MathEuclideanDistanceSquare)
    elseif (trim(clusterMethod)=="kmedoids") then
      if(procId == 0) call LogInfo("Running Kmedoids Cluster Algorithm")
      call ClusterKmedoids(normalized_cvs, cvDegree, ncluster, cvFrameNumber, clusterCycle, stat, MathEuclideanDistanceSquare)
    elseif (trajType /= 0 .and. trim(clusterMethod) == "coordinate") then
      if(procId == 0) call LogInfo("Running Kmedoids Coordinate Cluster Algorithm")
      call GetRMSDAtomNumber(rmsdAtomNumber)
      if(atomDegree /= rmsdAtomNumber * 3) then
        write(*, "(A, I4, A)") "Using subset Atom Index for RMSD Calcualtion, contains ", rmsdAtomNumber, " atoms: "
        write(*, "(10I6)") clusterAtomIndex(1:rmsdAtomNumber)
        allocate(subsetTrajIndex(rmsdAtomNumber*3))
        do i = 1, rmsdAtomNumber
            subSetTrajIndex(3*i - 2) = 3 * clusterAtomIndex(i) - 2
            subSetTrajIndex(3*i - 1) = 3 * clusterAtomIndex(i) - 1
            subSetTrajIndex(3*i) = 3 * clusterAtomIndex(i)
        enddo
        call ClusterKmedoids(traj(subsetTrajIndex, :), rmsdAtomNumber*3, ncluster, cvFrameNumber, clusterCycle, stat, MathAtomRMSD)
        deallocate(subsetTrajIndex)
      else 
        call ClusterKmedoids(traj , atomDegree, ncluster, cvFrameNumber, clusterCycle, stat, MathAtomRMSD)
      endif
    else
       if(procId == 0) call ErrorMessage("clusterMethod must choose be kmeans, kmedoids or coordinate")
    endif

    if(procId == 0) then
      write(*, "(a, I4,2x,a)") trim(clusterMethod) // " has ", ncluster, "clusters"
      if (stat >= clusterCycle) then
          write(*, "(a,I4,2x,a)")"clustering has run", stat, "iterations, maybe increase clusterCycle larger"
      else
          write(*, "(a,I4,2x,a)")"clustering  has run", stat, "iterations and converged"
      end if
      call LogInfo()
    endif
    call mpi_win_fence(0, win_normalized_shared, ierr)
    call mpi_win_free(win_normalized_shared, ierr)
    ! call MPISharedMemoryDeallocate()
  end subroutine ClusterAnalysis

  subroutine ModCluster(inputFile, resultFile)
    use mpi

    integer :: ioFile, ierr, i, numFile, shared_win
    real*8 :: start, finish
    character(256) :: inputFile, resultFile
    character(256), dimension(100) :: dataFile
    namelist /cluster/ ncluster, clusterCycle, clusterMethod, dataFile, trajType, randomSeed, clusterAtomIndex

    clusterCycle = 400
    trajType = 0
    dataFile = ""
    randomSeed = -1
    clusterAtomIndex = 0

    if (procId == 0) then
      call CPU_TIME(start)
      call GetFreeUnit(ioFile)
      OPEN(unit = ioFile, file=trim(inputFile), action="read")
      READ(ioFile, nml=cluster, iostat=ierr)
      if (ierr < 0 ) call ErrorMessage("error in reading the cluster namelists")
      close(ioFile)

      if (randomSeed == -1) then
        call init_random_seed()
      else
        call RestoreRandomSeed(randomSeed)
      endif

      call GetFileNumber(dataFile, numFile)

      if (trajType == 0) then
        call GetCvFileInfo(dataFile, numFile, cvFrameNumber, cvDegree, 1)        
      else if (trajType == 2) then
        call GetNetcdfInfo(dataFile, numFile, cvFrameNumber, atomDegree)
        call GetRMSDAtomNumber(rmsdAtomNumber)
      endif
    endif

    call mpi_bcast(ncluster, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(cvDegree, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(atomDegree, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(cvFrameNumber, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(clusterMethod, len(clusterMethod), mpi_char, 0, mpi_comm_world, ierr)
    call mpi_bcast(clusterCycle, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(trajType, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(numFile, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(dataFile, 256*numFile, mpi_char, 0, mpi_comm_world, ierr)

    if(ncluster < procNum) STOP "procNum must be less than the ncluster"
    if(cvFrameNumber < procNum) STOP "procNum must be less than the cvFrameNumber"

    call InitalizeSharedMemory(mpi_comm_world)

    if (trajType == 0) then
      call MPISharedMemoryAllocate2DimDoubleArray((/cvDegree, cvFrameNumber/), cvs, shared_win)
    else
      call MPISharedMemoryAllocate2DimDoubleArray((/atomDegree, cvFrameNumber/), traj, shared_win)
    endif

    if (shared_id == 0) then
      if (trajType == 0) then
        call CheckFileAndOpen(dataFile(1), ioFile)
        read(ioFile, *) cvs
        call CheckFileAndClose(ioFile)
      else if(trajType == 2) then
        call LogInfo("Running Coordinate Cluster Algorithm")
        call NetcdfReadCoordinate(dataFile, numFile, atomDegree, traj, 1)
      endif
    endif

    call mpi_win_fence(0, shared_win, ierr)
    call ClusterAnalysis()
    if (procId == 0) then
      call GetFreeUnit(ioFile)
      open(unit = ioFile, file=trim(resultFile), action="write")
      do i = 1, cvFrameNumber
          write(ioFile, *) maptocluster(i), centersnaps(maptocluster(i))
      enddo
      close(ioFile)
      resultDir="./"
      if (trajType==2) call ClusterWriteCoordianteResultToNetcdfFile(traj, atomDegree, cvFrameNumber, [(i,i=1,ncluster)], resultDir)
      call CPU_TIME(finish)
      write(*, "(a, f8.2, a)") "The cluster program has run ", finish-start, " seconds"
    endif
  end subroutine

  subroutine ClusterKmeansPlusPlusInitCentroids(coor, cvFrameNumber, cvDegree, ncluster, DistanceFunction, clusterCenterIndex)
    use mpi
    integer, intent(in) :: cvFrameNumber, cvDegree, ncluster
    real*8, intent(in) :: coor(cvDegree, cvFrameNumber)
    integer,intent(out)  :: clusterCenterIndex(ncluster)
    external :: DistanceFunction

    integer :: i, j, left_snap, right_snap, snap_length, ierr
    real*8 :: weights(cvFrameNumber)
    real*8 :: ran, tot, dist2, sub_tot
    real*8 , allocatable :: sub_dist2s(:), sub_weights(:)
    integer, dimension(procNum) :: snap_displs, snap_counts

    if(procId == 0) then
      write(*, "(a)") "Using kmeans++ to select initial center "
    endif

    call PartitionCPUProcessors(procNum, cvFrameNumber, snap_displs, snap_counts, left_snap, right_snap)
    snap_length = right_snap - left_snap + 1
    allocate(sub_dist2s(snap_length), sub_weights(snap_length))

    if(procId == 0) then
      call random_number(ran)
      j = int(ran*cvFrameNumber) + 1
    endif

    call mpi_bcast(j, 1, mpi_integer, 0, mpi_comm_world, ierr)
    clusterCenterIndex(1) = j
    sub_weights = huge(1.0)

    do i = 2, ncluster
      sub_tot = 0
      do j = left_snap, right_snap
        call DistanceFunction(coor(:, j), coor(:, clusterCenterIndex(i-1)), dist2, cvDegree)
        if (dist2 < sub_weights(j-left_snap+1)) sub_weights(j-left_snap+1) = dist2
        sub_tot = sub_tot + sub_weights(j-left_snap+1)
      enddo
      call mpi_gatherv(sub_weights, snap_length, mpi_double, weights, &
                      snap_counts, snap_displs, mpi_double, 0, mpi_comm_world, ierr)
      call mpi_reduce(sub_tot, tot, 1, mpi_double, mpi_sum, 0, mpi_comm_world, ierr)
      
      if (procId == 0) then
        call random_number(ran)
        ran = ran * tot
        tot = 0
        do j = 1, cvFrameNumber
          tot = tot + weights(j)
          if (tot > ran) exit
        enddo
      endif
      call mpi_bcast(j, 1, mpi_integer, 0, mpi_comm_world, ierr)
      clusterCenterIndex(i) = j
      if(procId == 0) then
        if( mod (i, 20) == 0) write(*, *) "Kmeans++ init subroutine generates ", i, " centers"
      endif
    enddo
    deallocate(sub_weights, sub_dist2s)
  end subroutine

  subroutine ClusterKmeans(coor, cvDegree, ncluster, cvFrameNumber, niter, stat, DistanceFunction)
    use mpi
    integer,intent(in) :: ncluster, cvFrameNumber, niter, cvDegree
    external :: DistanceFunction
    integer :: i, j, k, last_center, now_center, stat
    real*8 :: dist2, maxnum, minnum, tolerance, center_shift
    real*8 :: coor(cvDegree, cvFrameNumber), coor_center(cvDegree, ncluster), temp_dist(ncluster)
    real*8 :: old_coor_center(cvDegree, ncluster)
    real*8, allocatable :: sub_coor(:, :)
    integer, allocatable :: sub_maptocluster(:)
    real*8 :: subsums(cvDegree, ncluster), sums(cvDegree, ncluster)
    integer :: subnumbers(ncluster), numbers(ncluster)
    logical :: continue_run
    integer, DIMENSION(procNum) :: displs, counts
    integer :: left, right, length
    integer :: ierr

    tolerance = 1e-5

    allocate(maptocluster(cvFrameNumber), centersnaps(ncluster), numberincluster(ncluster))
    ! allocate(centersnaps(ncluster))
    call ClusterKmeansPlusPlusInitCentroids(coor, cvFrameNumber, cvDegree, ncluster, DistanceFunction, centersnaps)

    if (procId == 0) then
      do i = 1, ncluster
        coor_center(:, i) = coor(:, centersnaps(i))
      enddo
    endif

    call PartitionCPUProcessors(procNum, cvFrameNumber, displs, counts, left, right)

    length = right - left + 1
    allocate(sub_coor(cvDegree, length), sub_maptocluster(length))

    ! broadcast the coordinate and scatter the data
    call mpi_bcast(coor_center, ncluster*cvDegree, mpi_double, 0, mpi_comm_world, ierr)
    call mpi_scatterv(coor, cvDegree*counts, cvDegree*displs, mpi_double, sub_coor, & 
                     length*cvDegree, mpi_double, 0, mpi_comm_world, ierr)

    old_coor_center = coor_center
    j = 0 

    do while (j < niter)
      j = j + 1
      call AssignSnapsToKmeansCentroids(sub_coor, coor_center, sub_maptocluster, subSums, subnumbers, DistanceFunction)
      call mpi_gatherv(sub_maptocluster, length, mpi_int, maptocluster, counts, displs, mpi_int, 0, mpi_comm_world, ierr)
      call mpi_reduce(subsums, sums, ncluster*cvDegree, mpi_double, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(subnumbers, numberincluster, ncluster, mpi_int, mpi_sum, 0, mpi_comm_world, ierr)
      if (procId == 0) then
        center_shift = 0d0
        do i = 1, ncluster
          coor_center(:, i) = sums(:, i) / numberincluster(i)
          call DistanceFunction(coor_center(:, i), old_coor_center(:, i), dist2, cvDegree)
          center_shift = center_shift + sqrt(dist2)
        enddo
        if(mod(j, 20) == 0) then
          write(*, "(a, I5, a, F8.5)") "kmeans has run", j ," iterations, center shift is ", center_shift**2
        endif
      endif
      call mpi_bcast(coor_center, cvDegree*ncluster, mpi_double, 0, mpi_comm_world, ierr)
      call mpi_bcast(center_shift, 1, mpi_double, 0, mpi_comm_world, ierr)
      old_coor_center = coor_center
      if(center_shift**2 < tolerance) exit
    enddo

    ! finalize the kmeans cluster
    if(procId == 0) then
      temp_dist = 1e9; stat = j
      do i = 1, cvFrameNumber
        call DistanceFunction(coor(:, i), coor_center(:, maptocluster(i)), dist2, cvDegree)
        if (dist2 < temp_dist(maptocluster(i))) then
            centersnaps(maptocluster(i)) = i
            temp_dist(maptocluster(i)) = dist2
        end if
      enddo
    endif
  end subroutine

  subroutine AssignSnapsToKmeansCentroids(data, center, assignment, sums, number, DistanceFunction)
    real*8 :: data(:, :), center(: ,:), sums(:, :)
    external :: DistanceFunction
    integer :: number(:)
    integer :: assignment(:), i, j
    real*8 :: min_dist
    real*8 :: dist2

    sums = 0d0
    number = 0
    do i = 1, size(data, 2)
        min_dist = 1e9
        do j = 1, size(center, 2)
            call DistanceFunction(data(:, i), center(:, j), dist2, cvDegree)
            if (dist2 < min_dist) then
                min_dist = dist2
                assignment(i) = j
            endif
        enddo
        sums(:, assignment(i)) = sums(:, assignment(i)) + data(:, i)
        number(assignment(i)) = number(assignment(i)) + 1
    enddo
  end subroutine

  subroutine ClusterKmedoids(coor, cvDegree, ncluster, cvFrameNumber, niter, stat, DistanceFunction)
    use mpi
    integer, intent(in) :: ncluster, cvFrameNumber, niter, cvDegree 
    real*8,  intent(in) :: coor(cvDegree, cvFrameNumber)
    integer, intent(out) :: stat
    external :: DistanceFunction

    integer :: i
    real*8 ::  percent

    allocate(maptocluster(cvFrameNumber), centersnaps(ncluster), numberincluster(ncluster))

    if (cvFrameNumber<1000) then ! calculate distance matrix
      call ClusterKmedoidsSmallData(coor, cvDegree, ncluster, cvFrameNumber, niter, stat, DistanceFunction)
    else
      if (procId == 0) then
        write(*, *)"Frame Number is too large, unable to construct distance matrix, Using CLARANS method"
      endif
      percent = 500 / dble(cvFrameNumber)
      !percent = 1
      if(procId == 0) then 
        print*, "Maximum Generate number of test point is :" , int(percent * cvFrameNumber)
      endif
      call  ClusterKmedoidsLargeData(coor, cvDegree, ncluster, cvFrameNumber, niter, percent, stat, DistanceFunction)
    end if
    do i = 1, ncluster
      numberincluster(i) = count(maptocluster == i)
    enddo
  end subroutine ClusterKmedoids

  subroutine ClusterKmedoidsSmallData(coor, cvDegree, ncluster, cvFrameNumber, niter, stat, DistanceFunction)
    use util
    use mpi
    integer, intent(in):: cvDegree, ncluster, cvFrameNumber, niter
    integer, intent(out) :: stat
    real*8, intent(in) :: coor(cvDegree, cvFrameNumber)
    external :: DistanceFunction

    real*8, pointer :: distmat(:, :)
    real*8 :: min_value
    integer :: shared_distmat_win, left, right, ierr
    integer, dimension(procNum) :: counts, displs
    integer, allocatable :: oneclusterindex(:)
    integer :: i, j, k, temparray(cvFrameNumber), temparray2(cvFrameNumber), nowcentersnap(ncluster)

    call PartitionCPUProcessors(procNum, cvFrameNumber, displs, counts, left, right)
    call MPISharedMemoryAllocate2DimDoubleArray([cvFrameNumber, cvFrameNumber], distmat, shared_distmat_win)

    do i = left, right
      do j = i + 1, cvFrameNumber
        call DistanceFunction(coor(:, i), coor(:, j), distmat(i, j), cvDegree)
        distmat(j, i) = distmat(i, j)
      enddo
    enddo

    call mpi_win_fence(0, shared_distmat_win, ierr)
    call ClusterKmeansPlusPlusInitCentroids(coor, cvFrameNumber, cvDegree, ncluster, DistanceFunction, centersnaps)
    if(procId ==  0) then
      temparray = (/(i, i=1,cvFrameNumber)/)
      stat=niter
      do i=1, niter
        maptocluster = minloc(distmat(centersnaps, :), dim=1)
        do j=1, ncluster
          allocate(oneclusterindex(count(maptocluster==j)))
          oneclusterindex = pack(temparray, maptocluster==j)
          nowcentersnap(j) = oneclusterindex(minloc(sum(distmat(oneclusterindex, oneclusterindex), dim=2), dim=1))
          deallocate(oneclusterindex)
        enddo
        call isort(nowcentersnap, temparray2, 1, ncluster)
        call isort(centersnaps, temparray2, 1, ncluster)
        if(all(nowcentersnap == centersnaps)) then
            stat = i
            exit
        else
            centersnaps = nowcentersnap
        endif
      enddo
    endif
    ! deallocate(distmat)
    call mpi_win_fence(0, shared_distmat_win, ierr)
    call mpi_win_free(shared_distmat_win, ierr)
  end subroutine

  subroutine ClusterKmedoidsLargeData(coor, cvDegree, ncluster, cvFrameNumber, niter, percent, stat, DistanceFunction)
    use mpi
    integer, intent(in) :: ncluster, cvFrameNumber, cvDegree, niter
    real*8,  intent(in) :: coor(cvDegree, cvFrameNumber), percent
    integer, intent(out) :: stat
    external :: DistanceFunction

    integer :: newCenterSnaps(ncluster)
    integer :: i, ierr
    real*8 :: time

    call ClusterKmeansPlusPlusInitCentroids(coor, cvFrameNumber, cvDegree, ncluster, DistanceFunction, centersnaps)

    i = 0; stat=0
    do while(i < niter)
      if(procId == 0 .and. i > 0) write(*, *) "kmedoids has run ", i, "/", niter, " steps" 
      call AssignSnapToKmedoidsCentroids(coor, ncluster, cvFrameNumber, cvDegree, DistanceFunction)
      call KmedoidsCalculateCentroids(coor, ncluster, cvFrameNumber, cvDegree, newCenterSnaps, percent, DistanceFunction)
      if (all(centersnaps == newCenterSnaps)) exit
      centersnaps = newCenterSnaps
      i = i + 1
      stat = i
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
  end subroutine
  
  subroutine AssignSnapToKmedoidsCentroids(coor, ncluster, cvFrameNumber, cvDegree, DistanceFunction)
    use mpi
    integer, intent(in) :: ncluster, cvFrameNumber, cvDegree
    real*8, intent(in) :: coor(cvDegree, cvFrameNumber)
    external :: DistanceFunction

    integer, allocatable :: sub_maptocluster(:)
    integer :: i, j, ierr
    integer :: left_snap, right_snap, snap_length
    integer, dimension(procNum) :: snap_displs, snap_counts
    real*8 :: min_value, dist2s


    ! partition snap
    call PartitionCPUProcessors(procNum, cvFrameNumber, snap_displs, snap_counts, left_snap, right_snap)
    snap_length = right_snap - left_snap + 1

    allocate(sub_maptocluster(snap_length))

    do i = left_snap, right_snap
      min_value = huge(1.0)
      do j = 1, ncluster
        call DistanceFunction(coor(:, i), coor(:, centersnaps(j)), dist2s, cvDegree)
        if(dist2s < min_value)  then
          min_value = dist2s
          sub_maptocluster(i-left_snap+1) = j
        endif
      enddo
    enddo
    call mpi_allgatherv(sub_maptocluster, snap_length, mpi_integer, maptocluster, snap_counts, &
                        snap_displs, mpi_integer, mpi_comm_world, ierr)
    deallocate(sub_maptocluster)
  end subroutine

  subroutine KmedoidsCalculateCentroids(coor, ncluster, cvFrameNumber, cvDegree, &
    newCentersnaps, percent, DistanceFunction)
    use mpi
    real*8, intent(in) :: coor(cvDegree, cvFrameNumber), percent
    integer, intent(in) :: ncluster, cvFrameNumber, cvDegree
    integer, intent(out) :: newCenterSnaps(ncluster)
    external :: DistanceFunction
    integer, allocatable :: sub_centersnaps(:)
    integer :: cluster_displs(ncluster), cluster_counts(ncluster)
    integer :: i, j, ierr
    integer :: left_cluster, right_cluster, cluster_length

    call PartitionCPUProcessors(procNum, ncluster, cluster_displs, cluster_counts, left_cluster, right_cluster)
    cluster_length = right_cluster - left_cluster + 1

    allocate(sub_centersnaps(cluster_length))

    do i = left_cluster, right_cluster
      call KmedoidsCalculateSingleClusterCentroids(i, coor, cvFrameNumber, &
      sub_centersnaps(i-left_cluster+1), cvDegree, percent, DistanceFunction)
    enddo

    call mpi_allgatherv(sub_centersnaps, cluster_length, mpi_integer, newCenterSnaps, &
                        cluster_counts, cluster_displs, mpi_integer, mpi_comm_world, ierr)
    
    deallocate(sub_centersnaps)
  end subroutine

  subroutine KmedoidsCalculateSingleClusterCentroids(icluster, coor, cvFrameNumber, &
      icentersnap, cvDegree, percent, DistanceFunction)
    integer, intent(in) :: icluster,  cvFrameNumber, cvDegree
    real*8, intent(in) :: percent
    real*8, intent(in) :: coor(cvDegree, cvFrameNumber)
    integer, intent(out) :: icentersnap
    external :: DistanceFunction

    real*8 :: testscore, score
    integer, allocatable :: testPointsIndex(:)
    integer :: j, swapIndex
    integer :: nTestPoints

    nTestPoints = int(percent*cvFrameNumber)
    allocate(testPointsIndex(nTestPoints))
    ! do i = 1, ncluster
    call GenerateNeighbors(nTestPoints, icluster, testPointsIndex, coor, cvFrameNumber, cvDegree)
    call DistanceScore(coor, cvFrameNumber, cvDegree, icluster, centersnaps(icluster), score, DistanceFunction)
    swapIndex = 0

    do j = 1, nTestPoints
      call DistanceScore(coor, cvFrameNumber, cvDegree, icluster, testPointsIndex(j), testscore, DistanceFunction)
      if(testScore < score) then
        score = testScore
        swapIndex = j
      endif
    enddo
    if(swapIndex /= 0) then
      icentersnap = testPointsIndex(swapIndex)
    else
      icentersnap = centersnaps(icluster)
    endif
    deallocate(testpointsIndex)
  end subroutine

  subroutine GenerateNeighbors(nPoints, clusterIndex, pointsIndex, coor, cvFrameNumber, cvDegree)
    integer, intent(inout) :: nPoints
    integer, intent(in) ::  cvFrameNumber, cvDegree, clusterIndex
    real*8,  intent(in) ::  coor(cvDegree, cvFrameNumber)
    integer, intent(out) :: pointsIndex(nPoints)

    integer, allocatable :: subIndex(:)
    real*8 :: temp_real
    integer :: i, j, count, temp_int

    count = 0
    do i = 1, cvFrameNumber
      if (maptocluster(i) == clusterindex) then
        count= count+ 1
      endif
    enddo
    allocate(subIndex(count))
    if (count < nPoints) nPoints = count

    j = 1
    do i = 1, cvFrameNumber
        if (maptocluster(i) == clusterindex) then
            subIndex(j) = i
            j = j + 1
        endif
    enddo

    ! shuffle array
    do i = count, 1, -1
      call RANDOM_NUMBER(temp_real)
      j = int(ceiling(temp_real * i))
      temp_int = subIndex(i)
      subIndex(i) = subIndex(j)
      subIndex(j) = temp_int
    enddo
    pointsIndex(1:nPoints) = subIndex(1:nPoints)
   deallocate(subIndex)

  end subroutine

  subroutine DistanceScore(coor, cvFrameNumber, cvDegree, clusterindex, testPointIndex, score, DistanceFunction)
    real*8, intent(in) :: coor(cvDegree, cvFrameNumber)
    integer, intent(in) :: cvFrameNumber, cvDegree
    real*8, intent(out) :: score
    external :: DistanceFunction
    integer :: i, testPointIndex, clusterindex
    real*8 :: dist2

    score = 0d0
    do i = 1, cvFrameNumber
      if(maptocluster(i) == clusterindex) then
        call DistanceFunction(coor(:, i), coor(:, testPointIndex), dist2, cvDegree)
        score = score + dist2
      endif
    enddo
  end subroutine

  subroutine ClusterWriteCoordianteResultToNetcdfFile(traj, cvDegree, cvFrameNumber, clusterIndex, resultDir)
    real*8, intent(in) :: traj(cvDegree, cvFrameNumber)
    integer, intent(in) :: cvDegree, cvFrameNumber
    integer, intent(in) :: clusterIndex(:)
    character(MaxStringLen) :: temp_char, resultDir
    integer :: w_ncid, w_coorDVID, traj_index_counts
    integer :: i, j

    do i = 1, size(clusterIndex)
      write(temp_char, "(I6)") clusterIndex(i)
      temp_char = trim(resultDir)//"cluster_"//trim(adjustl(temp_char)) // ".nc"
      write(*, *) "writing cluster ", clusterIndex(i), "into " // trim(temp_char)
      call NetcdfCreateFile(temp_char, w_ncid, w_coorDVID)
      traj_index_counts = 1
      do j = 1, cvFrameNumber
         if(maptocluster(j) == clusterIndex(i)) then
            call NetcdfWriteCoordinate(w_ncid, w_coorDVID, reshape(traj(:, j), (/3, cvDegree/3/)), traj_index_counts)
            traj_index_counts = traj_index_counts + 1
         endif
      enddo
      call close_NCfile(w_ncid)
   enddo
  end subroutine
end module