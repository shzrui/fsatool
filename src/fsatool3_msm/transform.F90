module transform
  use util
  implicit none
contains

  subroutine mod_transform(inputFile, resultFile)
    character(256) :: datafile, transformmethod, inputFile, resultFile
    integer :: cvFrameNumber, nfeature, ncomponent, ioFile, ierr, lagTime, i
    real*8, allocatable :: data(:, :), maparray(:, :)

    namelist /transform/ datafile, transformmethod, ncomponent, lagTime

    call LogInfo("Transformation")

    call GetFreeUnit(ioFile)
    open(unit=ioFile, file=trim(inputFile), action="read")
    read(ioFile, nml=transform, iostat=ierr)
    if(ierr < 0) call ErrorMessage("error in reading the transform namelist")
    close(ioFile)


    call GetFreeUnit(ioFile)
    open(unit=ioFile, file=trim(datafile), action="read")
    read(ioFile, *) cvFrameNumber, nfeature
    if(ncomponent > nfeature) STOP "ncomponent must not big than the features in the datafile"
    allocate(data(cvFrameNumber, nfeature), maparray(cvFrameNumber, ncomponent))
    do i = 1, cvFrameNumber
      read(ioFile, *) data(i, :)
    enddo
    close(ioFile)

    if(trim(transformmethod) == "pca") then
      call pca(data, cvFrameNumber, nfeature, ncomponent, maparray)
    else if(trim(transformmethod) == "tica") then
      call tica(data, cvFrameNumber, nfeature, ncomponent, maparray, lagTime)
    else
      call ErrorMessage("transformmethod must be pca or tica")
    endif

    open(unit=ioFile, file=trim(resultFile), action="write")
    do i = 1, cvFrameNumber
      write(ioFile, *) maparray(i, :)
    enddo
    close(ioFile)

    deallocate(data, maparray)
    call LogInfo()
  end subroutine

  subroutine pca(tempArray, cvFrameNumber, nfeature, ncomponent, maparray)
    integer, intent(in) :: cvFrameNumber, nfeature, ncomponent
    real*8, intent(inout) :: tempArray(cvFrameNumber, nfeature)
    real*8, intent(out) :: maparray(cvFrameNumber, ncomponent)

    integer :: i, j, info
    real*8 :: average(nfeature), cov(nfeature, nfeature)
    real*8 :: eigenVector(nfeature, nfeature)
    real*8 :: eigenValue_real(nfeature), eigenValue_img(nfeature), work(4*nfeature), tempmatrix(nfeature,nfeature)
    integer :: component_index(nfeature)

    do i = 1, nfeature
      average(i) = sum(tempArray(:,i))/cvFrameNumber
      tempArray(:, i) = tempArray(:, i) - average(i)
    enddo

    do i = 1, nfeature
      do j = 1, nfeature
          cov(i,j) = sum((tempArray(:, i) - average(i)) * (tempArray(:, j) - average(j)))/(cvFrameNumber - 1)
      enddo
    enddo
    call dgeev("N","V", nfeature, cov, nfeature, eigenValue_real, eigenValue_img, &
        tempmatrix, nfeature, eigenVector, nfeature, work, 4*nfeature, info)
    call TPTGetFeatureIndex(eigenValue_real, nfeature, component_index)
    do i = 1, ncomponent
      tempmatrix(:, i) = eigenVector(:, component_index(i))
    enddo
    maparray = matmul(tempArray, tempmatrix(:, 1:ncomponent))
    write(*, "(a,10f10.5)")"PCA eigenValue", eigenValue_real
    write(*, "(a)")"PCA eigenVector:"
    do i = 1, nfeature
      write(*, "(10f10.5)") eigenVector(i, :)
    enddo
  end subroutine pca

  subroutine tica(data, npoint, nfeature, ncomponent, maparray, lagTime)
    use math, only: MathMatrixInverse
    integer :: lagTime, npoint, nfeature, ncomponent, i, j, offset, info, component_index(nfeature)
    real*8 :: data(npoint, nfeature), cov(nfeature, nfeature), cov_tau(nfeature, nfeature), exx(nfeature, nfeature)
    real*8 :: data_mean(nfeature, 1), cov_inv(nfeature, nfeature)
    real*8 :: eigenValue_real(nfeature), eigenValue_img(nfeature), tempmatrix(nfeature, nfeature)
    real*8 :: eigenVector(nfeature,nfeature), work(4*nfeature), maparray(npoint, ncomponent)

    offset = npoint - lagTime
    exx = matmul(transpose(data(:offset, :)), data(:offset, :)) &
        + matmul(transpose(data(lagTime+1:, :)), data(lagTime+1:, :))
    exx = exx / dble(offset * 2)
    do i=1, nfeature
      data_mean(i, 1) = sum(data(:offset, i) + data(lagTime+1:, i)) / dble(offset * 2)
    enddo
    cov  = exx - matmul(data_mean, transpose(data_mean))
    exx = matmul(transpose(data(:offset, :)), data(lagTime+1:, :)) / dble(offset)
    exx = (transpose(exx) + exx) /2.0d0
    cov_tau = exx - matmul(data_mean, transpose(data_mean))

    call dpotrf("L", nfeature, cov, nfeature, info)
    do i = 1, nfeature
      do j = i+1, nfeature
          cov(i,j) = 0.0d0
      enddo
    enddo
    call MathMatrixInverse(cov, cov_inv)
    !call inv(cov, cov_inv, nfeature)
    !cov_inv = cov
    tempmatrix = matmul(cov_inv, matmul(cov_tau, transpose(cov_inv)))
    call dgeev("N","V", nfeature, tempmatrix, nfeature, eigenValue_real, eigenValue_img, &
        tempmatrix, nfeature, eigenVector, nfeature, work, 4*nfeature, info)
    eigenVector = matmul(transpose(cov_inv), eigenVector)

    do i = 1,nfeature
      eigenVector(:,i) = eigenValue_real(i) * eigenVector(:, i)
    enddo
    call TPTGetFeatureIndex(eigenVector, nfeature, component_index)
    do i = 1, ncomponent
      tempmatrix(:, i) = eigenVector(:, component_index(i))
    enddo
    do i = 1, nfeature
      data(:, i) = data(:, i) - data_mean(i, 1)
    end do
    maparray = matmul(data, tempmatrix(:, 1:ncomponent))
    write(*, "(a,10f10.5)")"tica eigenValue", eigenValue_real
    do i = 1, nfeature
      write(*, "(10f10.5)") eigenVector(i, :)
    enddo
  end subroutine tica

subroutine TPTGetFeatureIndex(tempArray, nfeature, sortIndex)
  integer :: nfeature
  real*8 :: temp, tempArray(nfeature)
  integer :: tempint, sortIndex(nfeature)
  integer :: i, j

  do i = 1, nfeature
     sortIndex(i) = i
  enddo
  do i = 1, nfeature
     do j = i+1, nfeature
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

end module