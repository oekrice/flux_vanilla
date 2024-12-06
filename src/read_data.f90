!My module to read in the parameters saved out by run.py, to replace the old input deck version which doesn't really work.
!Need to establish ALL of the parameters that the deck gives -- there are probably some that I haven't considered, and some which need hard-coding

MODULE read_data

  USE shared_data

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE import_parameters

  REAL(num), DIMENSION(30):: parameters(0:29)

  CHARACTER(LEN=64):: parameter_filename, parameter_root
  CHARACTER(LEN=3) :: run_string

  parameter_root = './parameters/variables'
  write (run_string,'(I3.3)') run_id
  parameter_filename = trim(trim(parameter_root)//trim(run_string)//'.txt')

  OPEN(1, FILE = parameter_filename)
  READ(1, *) parameters
  CLOSE(1)

  if (rank == 0) print*, 'Parameters read-in from filename ', parameter_filename

  !Lift all the parameters from the original deck code. Things are in a mixed-up order and some are hard-coded but that's fine for now.
       !nx
       nx_global = int(parameters(13))
       !ny
       ny_global = int(parameters(14))
       !nz
       nz_global = int(parameters(15))
       !nprocx
       nprocx = 0
       !nprocy
       nprocy = 0
       !nprocz
       nprocz = 0

       x_stretch = .false.
       y_stretch = .false.
       z_stretch = .false.

       nsteps = -1

       t_end = parameters(2)

       dt_snapshots = 1.0

       x_start = parameters(7)
       x_end = parameters(8)
       y_start = parameters(9)
       y_end = parameters(10)
       z_start = parameters(11)
       z_end = parameters(12)

       visc1 = 0.1_num
       visc2 = 0.5_num
       visc3 = 0.00001_num

       ResistiveMHD = .true.
       j_max = 0.0_num
       eta0 = 0.0_num
       eta_background = 0.0001_num

       dt_multiplier = 0.8_num
       gamma = 5.0_num/ 3.0_num

       obs = .false.

       include_neutrals = .false.
       partially_ionized = .false.
       newton_cooling = .false.
       rke = .true.

       data_dir = './Data/'
       restart = .false.
       restart_snapshot = 0

       !Boundaries

       xbc_left = periodic
       xbc_right  = periodic
       ybc_up = periodic
       ybc_down  = periodic
       zbc_front = other
       zbc_back = other
       farfield = .false.
       damping = .true.

       ndiags = int(parameters(4))
       nsnaps = int(parameters(3))

       machine_flag = int(parameters(16))

  END SUBROUTINE import_parameters


END MODULE read_data







