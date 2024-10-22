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
  USE deck
  USE welcome
  USE oliver_output

  IMPLICIT NONE

  INTEGER :: i = 0
  REAL(num) :: walltime_current,dwalltime

  CALL minimal_init !setup.f90
  CALL mpi_minimal_init !mpi_routines.f90
  CALL welcome_message !welcome.f90
  CALL Read_Deck !deck.f90
  CALL set_boundary_conditions !boundary.f90
  CALL mpi_initialise !mpi_routines.f90
  CALL After_Control !setup.f90
  CALL open_files    !setup.f90
  CALL grid !setup.f90
  time = 0.0_num
  IF (restart) THEN
    ! CALL equilibrium    ! define initial profiles !setup.f90
    ! bx = 0.0_num
    ! by = 0.0_num
    ! bz = 0.0_num
     grav = 1.0_num
     CALL restart_data !setup.f90
    ! CALL flux_rope    ! below equilibrium subroutine
  ELSE
     CALL equilibrium    ! define initial profiles !setup.f90
  ENDIF
  CALL set_boundary_conditions !boundary.f90
  CALL boundary_conditions !boundary.f90
!
  WHERE(rho.LT.1.0e-9)rho=1.0e-9
!
  CALL output_routines(i) !diagnostics.f90
  IF (obs) CALL obs_setup !observables.f90
  IF (rank .EQ. 0) THEN
     PRINT *,"Equilibrium run OK, starting main loop"
  ENDIF
  walltime_start=MPI_WTIME()

  ndiags = 500
  nsnaps = 500
  diag_num = 0
  snap_num = 0
  DO

     if (time .ge. t_end*float(diag_num)/float(ndiags)) then   
      
      print*, 'Diag number', diag_num, 'at time', time
      !CALL output_diags(diag_num)
      diag_num = diag_num + 1
  
    end if


     if (time .ge. t_end*float(snap_num)/float(nsnaps)) then   
      
      !print*, 'Diag number', diag_num, 'at time', time
      CALL save_snap(snap_num)
      snap_num = snap_num + 1
    end if
 
     IF (((i >= nsteps) .AND. (nsteps >= 0)) .OR. (time >= t_end)) EXIT
     i = i + 1
     CALL eta_calc !lagran.f90
     IF (include_neutrals) CALL neutral_fraction !emergence.f90
     IF (partially_ionized) CALL perpendicular_resistivity !emergence.f90
     CALL set_dt !diagnostics.f90
     CALL lagrangian_step !lagran.f90
     CALL eulerian_remap(i) !remap.f90
     IF (rke) CALL energy_correction !diagnostics.f90
     IF (any_open) THEN
        CALL open_bcs !openboundaries.f90
        CALL boundary_conditions !boundaries.f90
     END IF
     IF (include_neutrals) CALL neutral_fraction !emergence.f90
     !CALL output_routines(i) !diagnostics.f90
     IF (obs) CALL obs_store !observables.f90
     IF (newton_cooling) CALL newton_relax !emergence.f90
!
     WHERE(rho.LT.1.0e-9)rho=1.0e-9
!
  END DO

  IF (obs) CALL obs_close !observables.f90
  CALL mpi_close !mpiroutines.f90
  CALL close_files !setup.f90
  CALL MPI_BARRIER(comm,errcode)
  CALL MPI_FINALIZE(errcode)    

END PROGRAM lare3d

