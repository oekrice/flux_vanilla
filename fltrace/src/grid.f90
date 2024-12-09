!*******************************************************************************
MODULE grid
!*******************************************************************************
! Generates all the shared grid data, as grid3d does in the python code. Maybe additionally put shared grid in arrays in a different file - not sure yet. Most (if not all) of the grid arrays should not depend on the process, but I guess we'll see about that.
!*******************************************************************************
    USE shared_data
    USE netcdf

    IMPLICIT NONE

    contains

    subroutine establish_grid()
        !Establishes the grid arrays, which hopefully can be read in mainly from the netcdf. Perhaps.
        !Also read in the magnetic field here because may as well
        IMPLICIT NONE
        INTEGER:: ncid, vid
        character(len=4) :: snap_id
        call get_coordinates

        allocate(bx(0:nx,-1:ny+1,-1:nz+1))
        allocate(by(-1:nx+1,0:ny,-1:nz+1))
        allocate(bz(-1:nx+1,-1:ny+1,0:nz))

        write (snap_id,'(I4.4)') snap
        bfield_filename = trim(trim(data_root)//trim(snap_id)//'.nc')

        if (print_flag > 0.5) print*, 'Using magnetic field from file ', bfield_filename

        call try(nf90_open(trim(bfield_filename), nf90_nowrite, ncid))

        call try(nf90_inq_varid(ncid, 'bx', vid))
        call try(nf90_get_var(ncid, vid, bx(0:nx,1:ny,1:nz)))

        call try(nf90_inq_varid(ncid, 'by', vid))
        call try(nf90_get_var(ncid, vid, by(1:nx,0:ny,1:nz)))

        call try(nf90_inq_varid(ncid, 'bz', vid))
        call try(nf90_get_var(ncid, vid, bz(1:nx,1:ny,0:nz)))

        call try(nf90_close(ncid))

        !call check_divergence

        if (print_flag > 0.5) print*, 'Grid established and magnetic field read-in'

        ds = ds_factor*min(dx, dy, dz)  !Tracing 'timestep'

    end subroutine establish_grid

    subroutine get_coordinates
    !Allocate dimension arrays
    IMPLICIT NONE
    INTEGER:: i,j,k
    allocate(xs(0:nx),ys(0:ny),zs(0:nz))
    allocate(xc(1:nx),yc(1:ny),zc(1:nz))

    dx = (x1 - x0)/nx;     dy = (y1 - y0)/ny;     dz = (z1 - z0)/nz

    xs(0) = x0; ys(0) = ny; zs(0) = nz
    do i = 1, nx
        xs(i) = xs(i-1) + dx
    end do
    do j = 1, ny
        ys(j) = ys(j-1) + dy
    end do
    do k = 1, nz
        zs(k) = zs(k-1) + dz
    end do

    xc(1:nx) = 0.5_num*(xs(0:nx-1) + xs(1:nx))
    yc(1:ny) = 0.5_num*(ys(0:ny-1) + ys(1:ny))
    zc(1:nz) = 0.5_num*(zs(0:nz-1) + zs(1:nz))

    dx = sum((xs(1:nx) - xs(0:nx-1)))/ nx
    dy = sum((ys(1:ny) - ys(0:ny-1)))/ ny
    dz = sum((zs(1:nz) - zs(0:nz-1)))/ nz

    end subroutine get_coordinates

    subroutine check_divergence
    !Checks that the solenoidal condition is fine in every cell
    IMPLICIT none
    real(num):: div(1:nx,1:ny,1:nz)

    div = 0.0_num
    div(1:nx,1:ny,1:nz) = div(1:nx,1:ny,1:nz) + (bx(1:nx,1:ny,1:nz) - bx(0:nx-1,1:ny,1:nz))/dx
    div(1:nx,1:ny,1:nz) = div(1:nx,1:ny,1:nz) + (by(1:nx,1:ny,1:nz) - by(1:nx,0:ny-1,1:nz))/dy
    div(1:nx,1:ny,1:nz) = div(1:nx,1:ny,1:nz) + (bz(1:nx,1:ny,1:nz) - bz(1:nx,1:ny,0:nz-1))/dz

    if (maxval(abs(div)) > 1d-10) then
        print*, 'Read-in magnetic field is not divergence free, stopping'
        STOP
    end if

    end subroutine check_divergence

    subroutine export_fieldlines
    !Output the calculated magnetic field lines as a netcdf file flines.nc
    IMPLICIT NONE

    character(len=64):: filename
    integer:: aid, bid, cid, vid, ncid
    CHARACTER(LEN=4):: snap_id

    write (snap_id,'(I3.3)') snap

    filename = trim('./fl_data/flines'//trim(snap_id)//'.nc')

    call try(nf90_create(trim(filename), nf90_clobber, ncid))

    !Define variables
    call try(nf90_def_dim(ncid, 'a', 3, aid))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'b', max_line_length, bid))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'c', nstarts*2, cid))  !Make up fake dimensions here

    call try(nf90_def_var(ncid, 'lines', nf90_double, (/cid,bid,aid/), vid))
    call try(nf90_enddef(ncid))

    !Write variables
    call try(nf90_put_var(ncid, vid, all_lines))
    call try(nf90_close(ncid))

    if (print_flag > 0.5_num) print*, 'Field lines exported to file ', filename

    return

    end subroutine export_fieldlines

    SUBROUTINE try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        print*, 'Ensure data directory is correct in fltrace.f90'
        stop
    end if

    END SUBROUTINE try

END MODULE grid
