module NETCDF_FUNC
    use netcdf
    use AmberNetcdf_mod
    implicit none
    integer :: atom_number
    integer :: frameNumber
    logical :: has_box
contains

    subroutine NetcdfInitParameter(file, tempFrameNumber, tempAtomNumber)
        character(*), intent(in) :: file
        integer, intent(out), optional :: tempAtomNumber, tempFrameNumber
        character(256) :: title
        integer :: ncid, coordVID, velocityVID, timeVID, cellLengthVID, cellAngleVID, TempVID

        if (NC_openRead(file, ncid)) stop "error open"
        if (NC_setupMdcrd(ncid, title, frameNumber, atom_number, coordVID, &
            velocityVID, timeVID, cellLengthVID, cellAngleVID, TempVID)) stop "error setup  mdcrd"
        if (cellLengthVID == -1) has_box = .false.
!        write(*, "(A, I8, A)") trim(file) // " has ",  frameNumber , " frames"
        if (present(tempAtomNumber)) tempAtomNumber = atom_number
        if (present(tempFrameNumber)) tempFrameNumber = frameNumber
    end subroutine

    subroutine read_NCfile(file, ncid, coorDVID)
        character(*),intent(in) :: file
        integer, intent(out) :: ncid, coorDVID
        character(256) :: title
        integer :: timeVID, cellLengthVID, cellAngleVID, tempVID
        integer :: velocityVID, frameDID
        if (NC_openRead(file, ncid)) stop "error open"
        if (NC_setupMdcrd(ncid, title, frameNumber, atom_number, coordVID, &
            velocityVID, timeVID, cellLengthVID, cellAngleVID, TempVID)) stop "error setup  mdcrd"
    end subroutine
        

    subroutine NetcdfCreateFile(file, ncid, coorDVID)
        character(*),intent(in) :: file
        integer, intent(out) :: ncid, coorDVID
        integer :: timeVID, velocityVID, frcVID, cellLengthVID, cellAngleVID, tempVID
        integer :: frameDID
        if (NC_create(file, 'O', .false., atom_number, .true.,  &
                       .false., has_box, .false., .true., &
                       .false., "level file", &
                       ncid, timeVID, coordVID, velocityVID, frcVID,&
                       cellLengthVID, cellAngleVID, TempVID, &
                       frameDID)) stop "file create file"
    end subroutine
    
    subroutine NetcdfWriteCoordinate(ncid, coordVID, crd, mdcrd_frame)
        integer, intent(in) :: ncid, coordVID
        real*8, intent(in):: crd(:, :)
        integer, intent(in) :: mdcrd_frame
        call checkNCerror(nf90_put_var(ncid, coordVID, crd(:, :), &
                                       start=(/1,1,mdcrd_frame/), &
                                       count=(/3, atom_number, 1/)), "write atom coords")
    end subroutine
    
    subroutine NetcdfWriteCellInfo(ncid, cellLengthVID, cellAngleVID, length, Angle, mdcrd_frame)
        integer :: ncid, cellLengthVID, cellAngleVID 
        real*8 :: length(3), Angle(3)
        integer :: mdcrd_frame
        call checkNCerror(nf90_put_var(ncid, cellLengthVID, length, &
                              start = (/1, mdcrd_frame/), count = (/3,1/)), "write cell length")
        call checkNCerror(nf90_put_var(ncid, cellAngleVID, angle, &
                              start = (/1, mdcrd_frame/), count = (/3,1/)), "write cell length")
    end subroutine

    subroutine NetcdfGetCoordinate(ncid, coordVID, crd, frame)
        integer, intent(in) :: ncid, coordVID
        integer, intent(in) :: frame
        real*8, dimension(:, :), intent(out) :: crd
        call checkNCerror(nf90_get_var(ncid, coordVID, crd, start=(/1,1,frame/)), &
                    "get coor variable")
    end subroutine

    subroutine NetcdfFinishWriting(ncid)
        integer :: ncid
        call checkNCerror(nf90_sync(ncid))
    end subroutine

    subroutine close_NCfile(ncid)
        integer :: ncid
        call NC_close(ncid)
    end subroutine

    subroutine NetcdfReadCoordinate(datafile, num_file, atomDegree, traj, stride)
        integer,intent(in) :: num_file, atomDegree
        integer, intent(in) :: stride
        character(*), intent(in) :: datafile(num_file)
        real*8, intent(inout) :: traj(:, :)
        real*8, allocatable :: coor(:, :)
        integer :: frame_number_per_file(num_file)
        integer :: i, j, ncid, coordvid
        integer ::  cur_index, total_frame
  
        cur_index = 0
        total_frame = 0

        do i = 1, num_file
           call NetcdfInitParameter(datafile(i))
           frame_number_per_file(i)  =  frameNumber
           total_frame = total_frame + frameNumber
        enddo

        ! if (total_frame/cvStride /= cvFrameNumber) then 
        !     write(*, *) "frame number: ", total_frame, "cvFrameNumber:", cvFrameNumber
        !     STOP "Trajectorys frame number is not equal to cvFrameNumber"
        ! endif
        if (atom_number*3 /= atomDegree) then
            write(*, *) "atom_number * 3: ", atom_number * 3, "atomDegree:", atomDegree
            STOP "Trajectorys atom number * 3 is not equal to atomDegree"
        endif
        allocate(coor(3, atom_number))
        do i = 1, num_file
           call read_ncfile(datafile(i), ncid, coordvid)
           do  j = 1,  frame_number_per_file(i)
              cur_index = cur_index + 1
              if(mod(cur_index , stride) == 0)  then
                call NetcdfGetCoordinate(ncid, coordvid, coor, j)
                traj(:, cur_index/stride) = pack(coor, .true.)
              endif
           enddo
           write(*, "(a)") "finished reading file  " // trim(datafile(i))
        enddo
        deallocate(coor)
    end subroutine
end module