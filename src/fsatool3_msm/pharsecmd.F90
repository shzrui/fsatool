module getcmd
    use GlobalVariable, only: inputFile, resultDir, MaxStringLen
    implicit none
    character(MaxStringLen) :: argname, resultFile
    integer, save :: iarg = 0
    public :: get_program
contains
    subroutine get_program(msm)
        logical  :: msm
        msm = .false.
        !if(iargc() == 0) return;
        iarg = iarg + 1;
        call get_command_argument(iarg, argname)
        select case (trim(argname))
            case("cluster")
                call run_cluster()
            case("lumping")
                call run_lumping()
            case("tpt")
                call run_tpt()
            case("check")
                call run_check()
            case("tram")
                call run_tram()
            case ("transform")
                call run_transform()
            case("write")
                call run_writeTraj()
            case("msm")
                call run_msm()
                msm = .true.
            case default
                call run_msm()
                msm = .true.
        end select
    end subroutine get_program

    subroutine get_nextcmd()
        integer :: narg
        narg = iargc()
        do
            if(iarg >= narg) exit
            iarg = iarg + 1
            call get_command_argument(iarg, argname)
            if (trim(argname) == "-i") then
                iarg = iarg + 1
                call get_command_argument(iarg, argname)
                inputFile = argname
            elseif  (trim(argname) == "-o") then
                iarg = iarg + 1
                call get_command_argument(iarg, argname)
                resultFile = argname
            else
                continue
            endif
        enddo
    end subroutine

    subroutine run_cluster()
        use cluster, only: ModCluster
        inputFile = "cluster.in"
        resultFile = "cluster.out"
        call get_nextcmd()
        call check_file_exist()
        call ModCluster(inputFile, resultFile)
    end subroutine

    subroutine run_lumping()
        use CoarseGrain, only: ModCoarseGrain
        inputFile = "lumping.in"
        resultFile = "lumping.out"
        call get_nextcmd()
        call ModCoarseGrain(inputFile, resultFile)
    end subroutine

    subroutine run_tpt()
        use tpt, only: ModTPT
        inputFile = "tpt.in"
        resultFile = "tpt.out"
        call get_nextcmd()
        call check_file_exist()
        call ModTPT(inputFile, resultFile)
    end subroutine

    subroutine run_check()
        use markov, only: MarkovCheck
        inputFile = "check.in"
        resultFile = "check.out"
        call get_nextcmd()
        call check_file_exist()
        call MarkovCheck(inputFile, resultFile)
    end subroutine

    subroutine run_tram()
        use tram, only: mod_tram
        inputFile = "tram.in"
        resultFile = "tram.out"
        call get_nextcmd()
        call check_file_exist()
        call mod_tram(inputFile, resultFile)
    end subroutine

    subroutine run_transform()
        use transform, only: mod_transform
        inputFile = "transform.in"
        resultFile = "transform.out"
        call get_nextcmd()
        call check_file_exist()
        call mod_transform(inputFile, resultFile)
    end subroutine

    subroutine run_writeTraj()
        use fileio, only: mod_writeTraj
        inputFile = "writeTraj.in"
        resultDir = "info/"
        resultFile = ""
        call get_nextcmd()
        call check_file_exist()
        call mod_writeTraj(inputFile, resultDir)
    end subroutine

    subroutine run_msm()
        inputFile = "msm.in"
        resultDir = "info/"
        resultFile = ""
        call get_nextcmd()
        call check_file_exist()
        if (resultFile /= "") resultDir = resultFile
    end subroutine


    subroutine check_file_exist()
        logical :: file_exist
        inquire(file=inputFile, exist=file_exist)
        if( .not. file_exist) stop "input file not exist"
    end subroutine

end module getcmd
