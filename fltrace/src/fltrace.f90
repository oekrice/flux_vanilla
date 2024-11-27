!*******************************************************************************
PROGRAM fltrace
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code.
!*******************************************************************************
    USE grid
    USE shared_data

    IMPLICIT NONE

    CHARACTER(LEN=64):: input_value
    CHARACTER(LEN=64):: parameter_filename
    CHARACTER(LEN=4):: run_id


    call get_command_argument(1, input_value)
    read(unit=input_value,fmt=*) run

    write (run_id,'(I3.3)') run
    parameter_filename = trim('./fl_data/flparameters'//trim(run_id)//'.txt')

    print*, parameter_filename
    !Import parameters from text file (saved out by python in the fltrace directory)
    !Need to update this so it can do different numbers
    open(1,file= parameter_filename)
    read(1, *) flparameters
    close(1)

    !##########################################
    !DATA ROOT HERE. NEED TO READ IN AS STRING REALLY
    data_root = '../Data/'

    ! Put some of the major variables in here - things that can be changed occasionally but not in a series of runs

    run = int(flparameters(0))
    nx = int(flparameters(1))
    ny = int(flparameters(2))
    nz = int(flparameters(3))
    x0 = flparameters(4)
    x1 = flparameters(5)
    y0 = flparameters(6)
    y1 = flparameters(7)
    z0 = flparameters(8)
    z1 = flparameters(9)

    snap = int(flparameters(10))
    nstarts = int(flparameters(11))
    print_flag = int(flparameters(12))

    max_line_length = int(flparameters(13))
    ds_factor =  flparameters(14)
    weakness_limit =  flparameters(15)

    data_source = int(flparameters(16))

    if (data_source > 0) then
        if (data_source == 15) data_root = '../Data_15/'
        if (data_source == 150) data_root = '../Data_150/'
    end if
    call establish_grid()  !Establish the grid and read in the magnetic field

    call establish_starts()  !Import the location of the start points

    call bfield_boundary() !Populate ghost points so the interpolator works

    call integrate_fieldlines()

    call export_fieldlines()

    contains

    SUBROUTINE establish_starts

    !Read in the start points from the provided 'starts.txt'
    IMPLICIT NONE

    INTEGER:: i
    REAL(num), DIMENSION(:):: starts_import(0:nstarts*3-1)
    CHARACTER(LEN=64):: starts_filename
    CHARACTER(LEN=4):: run_id


    ALLOCATE(starts(0:nstarts-1,0:2))

    write (run_id,'(I3.3)') run

    print*, run, run_id
    starts_filename = trim('./fl_data/starts'//trim(run_id)//'.txt')

    open(1,file= starts_filename)
    read(1, *) starts_import
    close(1)

    do i = 0, nstarts-1
        starts(i,:) = starts_import(3*i:3*i+2)
    end do

    END SUBROUTINE establish_starts

    SUBROUTINE interpolate_bfield(x_point, y_point, z_point)
    !Finds the magnetic field vector at a given point, taking advantage of the staggered grid.
    IMPLICIT NONE

    REAL(num):: x_point, y_point, z_point !Coordinates of point to interpolate
    REAL(num):: xp, yp, zp !Coordinates in the respective dimentions
    INTEGER:: xi, yi, zi !Cell indices in each dimension
    REAL(num):: xf, yf, zf !Distance 'up' each cell

    b1 = 0.0_num
    !Establish ratios
    xp = nx*(x_point - x0)/(x1 - x0)
    yp = ny*(y_point - y0)/(y1 - y0)
    zp = nz*(z_point - z0)/(z1 - z0)

    !Interpolate bx
    xi = int(xp); yi = int(yp + 0.5_num); zi = int(zp + 0.5_num)
    xf = xp - xi; yf = yp + 0.5_num - yi; zf = zp + 0.5_num - zi

    b1(0) = b1(0) + bx(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + bx(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    b1(0) = b1(0) + bx(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + bx(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    b1(0) = b1(0) + bx(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + bx(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    b1(0) = b1(0) + bx(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + bx(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    !Interpolate by
    xi = int(xp+0.5_num); yi = int(yp); zi = int(zp + 0.5_num)
    xf = xp + 0.5_num - xi; yf = yp - yi; zf = zp + 0.5_num - zi

    b1(1) = b1(1) + by(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + by(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    b1(1) = b1(1) + by(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + by(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    b1(1) = b1(1) + by(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + by(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    b1(1) = b1(1) + by(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + by(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    !Interpolate bz
    xi = int(xp+0.5_num); yi = int(yp+0.5_num); zi = int(zp)
    xf = xp + 0.5_num - xi; yf = yp + 0.5_num - yi; zf = zp - zi

    b1(2) = b1(2) + bz(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + bz(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    b1(2) = b1(2) + bz(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + bz(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    b1(2) = b1(2) + bz(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + bz(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    b1(2) = b1(2) + bz(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + bz(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    !print*, bz(xi+1,yi+1,zi), bz(xi+1,yi+1,zi+1)
    !print*, xf, yf, zf,  b1
    !print*, '______________________'
    if (sqrt(sum(b1**2)) < weakness_limit) then
        null_point = 1
    end if

    b1 = b1/sqrt(sum(b1**2))

    END SUBROUTINE interpolate_bfield

    SUBROUTINE bfield_boundary

    !Populates the ghost points so the interpolator works. Actual values aren't really that important
    IMPLICIT NONE

    by(0:nx+1,0:ny,0) = by(0:nx+1,0:ny,1) - dz*(bz(0:nx+1,1:ny+1,0) - bz(0:nx+1, 0:ny,0))/dy
    bx(0:nx, 0:ny+1,0) = bx(0:nx,0:ny+1,1) - dz*(bz(1:nx+1,0:ny+1,0) - bz(0:nx,0:ny+1,0))/dx

    !UPPER BOUNDARY (Zero Current)
    by(0:nx+1,0:ny,nz+1) = by(0:nx+1,0:ny,nz) + dz*(bz(0:nx+1,1:ny+1,nz) - bz(0:nx+1, 0:ny,nz))/dy
    bx(0:nx, 0:ny+1,nz+1) = bx(0:nx,0:ny+1,nz) + dz*(bz(1:nx+1,0:ny+1,nz) - bz(0:nx,0:ny+1,nz))/dx

    !x boundaries (Zero current, and zero flux)
    bz(0,0:ny+1,0:nz) = bz(1,0:ny+1,0:nz) - dx*(bx(0,0:ny+1,1:nz+1) - bx(0, 0:ny+1,0:nz))/dz
    by(0,0:ny,0:nz+1) = by(1, 0:ny,0:nz+1) - dx*(bx(0,1:ny+1,0:nz+1) - bx(0,0:ny,0:nz+1))/dy

    bz(nx+1,0:ny+1,0:nz) = bz(nx,0:ny+1,0:nz) + dx*(bx(nx,0:ny+1,1:nz+1) - bx(nx, 0:ny+1,0:nz))/dz
    by(nx+1,0:ny,0:nz+1) = by(nx, 0:ny,0:nz+1) + dx*(bx(nx,1:ny+1,0:nz+1) - bx(nx,0:ny,0:nz+1))/dy

    !y boundaries (Zero current, and zero flux)
    bz(0:nx+1,0,0:nz) = bz(0:nx+1, 1,0:nz) - dy*(by(0:nx+1,0,1:nz+1) - by(0:nx+1,0,0:nz))/dz
    bx(0:nx,0,0:nz+1) = bx(0:nx,1,0:nz+1) - dy*(by(1:nx+1,0,0:nz+1) - by(0:nx, 0,0:nz+1))/dx

    bz(0:nx+1,ny+1,0:nz) = bz(0:nx+1, ny,0:nz) + dy*(by(0:nx+1,ny,1:nz+1) - by(0:nx+1,ny,0:nz))/dz
    bx(0:nx,ny+1,0:nz+1) = bx(0:nx,ny,0:nz+1) + dy*(by(1:nx+1,ny,0:nz+1) - by(0:nx, ny,0:nz+1))/dx

    END SUBROUTINE bfield_boundary

    SUBROUTINE integrate_line(start, updown, line_number)

    IMPLICIT NONE

    REAL(num), DIMENSION(:):: start(0:2), pt(0:2)
    REAL(num), DIMENSION(:,:):: line(0:max_line_length-1,0:2)
    INTEGER:: lcount, updown, line_number
    line = 1e6; lcount = 0

    pt = start; null_point = 0
    do while (.true.)
        if ((pt(0) < x0) .or. (pt(0) > x1) .or. (pt(1) < y0) .or. (pt(1) > x1) .or. &
        (pt(2) < z0) .or. (pt(2) > z1) .or. (lcount > max_line_length-1)) then
            exit
        end if
        if (null_point > 0.5_num) then
            line(1:max_line_length-1,0:2) = 1e6
            exit
        end if
        line(lcount,0:2) = pt
        call interpolate_bfield(pt(0), pt(1), pt(2))
        pt = pt + updown*ds*b1
        lcount = lcount + 1
    end do

    all_lines(line_number,:,:) = line

    END SUBROUTINE

    SUBROUTINE integrate_fieldlines()
    !Uses the 'starts' array to output an array of all the fieldlines
    !Traces in both directions
    IMPLICIT NONE

    integer:: start_index, line_number

    allocate(all_lines(0:2*nstarts-1, 0:max_line_length-1, 0:2))

    line_number = 0
    do start_index = 0, nstarts-1
        !Integrate opposite to magnetic field
        call integrate_line(starts(start_index,:),-1,line_number)
        line_number = line_number + 1
        !Integrate with magnetic field
        call integrate_line(starts(start_index,:),1,line_number)
        line_number = line_number + 1
    end do

    if (print_flag > 0.5_num) print*, nstarts*2, 'Field lines integrated'

    END SUBROUTINE integrate_fieldlines








END PROGRAM fltrace







