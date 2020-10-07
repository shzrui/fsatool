program main
  use GlobalVariable
  use cluster
  use fileio
  use cluster
  use markov
  use tpt
  use getcmd
  use mpi

  implicit none
  integer :: ierr
  logical :: msm

  call GlobalVariableMPIInitialize()
  call get_program(msm)

  if (msm) then
    call fileio_init_parameters()
    call GlobalVariableReadTrajs()
    if (ifReadCoordinate .eqv. .false.) then
        call fileio_readclusterinfo()
    else
        call ClusterAnalysis()
        call fileio_writeclusterinfo()
    endif
    if(.not. ifOnlyRunClustering) then
      call MarkovAnalysis()
      call fileio_writestateinfo()
    endif
    call TPTAnalysis()
  endif
  call mpi_finalize(ierr)
end program