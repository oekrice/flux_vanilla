PROGRAM lare3d

  USE shared_data
  USE setup
  USE boundary
  USE diagnostics
  USE lagran
  USE remap
  USE mpi_routines
  USE observables
  USE emergence
  USE initial_conditions
  USE read_data
  !USE deck
  USE welcome
  USE oliver_output

  IMPLICIT NONE

  INTEGER :: i = 0
  REAL(num) :: walltime_current,dwalltime
  CHARACTER(LEN=64):: input_value

  CALL minimal_init !setup.f90
  CALL mpi_minimal_init !mpi_routines.f90
  CALL welcome_message !welcome.f90

  call get_command_argument(1, input_value)
  read(unit=input_value,fmt=*) run_id

  if (rank == 0) print*, 'Run id: ', run_id
  CALL import_parameters

  CALL set_boundary_conditions !boundary.f90
  CALL mpi_initialise !mpi_routines.f90
  CALL After_Control !setup.f90
  CALL grid !setup.f90
  time = 0.0_num
  CALL equilibrium    ! define initial profiles !setup.f90
  CALL set_boundary_conditions !boundary.f90
  CALL boundary_conditions !boundary.f90

  if (machine_flag < 0.5_num) then
    data_dir = '/home/grads/trcn27/rdata/flux_emergence/'
  else if (machine_flag  < 1.5_num) then
    data_dir = '/nobackup/trcn27/flux_emergence/'
  else
    data_dir = './Data'
  end if

  if (rank == 0) print*, 'Data directory: ', data_dir
  WHERE(rho.LT.1.0e-9)rho=1.0e-9
! !
   IF (rank .EQ. 0) THEN
      PRINT *,"Equilibrium run OK, starting main loop"
   ENDIF
   walltime_start=MPI_WTIME()
   diag_num = 0
   snap_num = 0

   DO
      if (time .ge. t_end*float(diag_num)/float(ndiags)) then
      !print*, 'Diagnostic number', diag_num, 'at time', time
      ! CALL output_diags(diag_num)
       diag_num = diag_num + 1
      end if

      if (time .ge. t_end*float(snap_num)/float(nsnaps)) then

       CALL save_snap(snap_num)
       snap_num = snap_num + 1
     end if

      IF (((i >= nsteps) .AND. (nsteps >= 0)) .OR. (time >= t_end)) EXIT
      i = i + 1

      CALL eta_calc !lagran.f90
!      IF (include_neutrals) CALL neutral_fraction !emergence.f90
!      IF (partially_ionized) CALL perpendicular_resistivity !emergence.f90

      CALL set_dt !diagnostics.f90
      CALL lagrangian_step !lagran.f90
      CALL eulerian_remap(i) !remap.f90
      IF (rke) CALL energy_correction !diagnostics.f90

      energy = energy +  correction_factor*dt*(energy_reference-energy)  !Energy correction factor for the Newton Cooling

!      IF (any_open) THEN
!         CALL open_bcs !openboundaries.f90
!         CALL boundary_conditions !boundaries.f90
!      END IF
!      IF (include_neutrals) CALL neutral_fraction !emergence.f90
!      !CALL output_routines(i) !diagnostics.f90
!      IF (newton_cooling) CALL newton_relax !emergence.f90
! !
      WHERE(rho.LT.1.0e-9)rho=1.0e-9
! !
   END DO
!
   CALL mpi_close !mpiroutines.f90
   CALL close_files !setup.f90
  CALL MPI_BARRIER(comm,errcode)
  CALL MPI_FINALIZE(errcode)    

END PROGRAM lare3d

