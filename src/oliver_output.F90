
MODULE oliver_output

  USE shared_data
  USE netcdf

  IMPLICIT NONE

  
CONTAINS

  
SUBROUTINE establish_grid()

    ! Imports the initial condition from the inits file.
    ! Will try to be smart with this, such that the vector potential is only read in in individual processes

    CHARACTER(LEN =64):: init_filename
    INTEGER:: ncid, vid, status

    write (init_filename, "(A12, A2, I1, A3)") './inits/init', '00', 0, '.nc'
   
    print*, trim(init_filename)

    status = nf90_open(trim(init_filename), nf90_nowrite, ncid)
    !call try(nf90_open(trim(init_filename), nf90_nowrite, ncid))

    !call try(nf90_close(ncid))

END SUBROUTINE establish_grid


SUBROUTINE save_snap(snap_num)
    !Exports the magnetic field at this plot_num to an appropriate netcdf file
    IMPLICIT NONE

    CHARACTER(LEN =64):: output_filename
    INTEGER:: snap_num, proc_write
    INTEGER:: ncid, vid
    INTEGER:: xs_id, ys_id, zs_id
    INTEGER:: xc_id, yc_id, zc_id
    INTEGER:: bx_id, by_id, bz_id
    INTEGER:: en_id, rho_id
    INTEGER:: vx_id, vy_id, vz_id


    print*, 'lbound flux', sum(abs(bz(1:nx,1:ny,0)))*dxb(1)*dyb(1)
    if (snap_num < 10) then
        write (output_filename, "(A7,A3,I1,A3)") "./Data/", "000", snap_num, ".nc"
    else if (snap_num < 100) then
        write (output_filename, "(A7,A2,I2,A3)") "./Data/", "00", snap_num, ".nc"
    else if (snap_num < 1000) then
        write (output_filename, "(A7,A1,I3,A3)")  "./Data/", "0", snap_num, ".nc"
    else if (snap_num < 10000) then
        write (output_filename, "(A7,I4,A3)")  "./Data/", snap_num, ".nc"
    end if

    if (rank == 0) then
    call try(nf90_create(trim(output_filename), nf90_clobber, ncid))

    call try(nf90_def_dim(ncid, 'xs', nx_global+1, xs_id))
    call try(nf90_def_dim(ncid, 'ys', ny_global+1, ys_id))
    call try(nf90_def_dim(ncid, 'zs', nz_global+1, zs_id))

    call try(nf90_def_dim(ncid, 'xc', nx_global, xc_id))
    call try(nf90_def_dim(ncid, 'yc', ny_global, yc_id))
    call try(nf90_def_dim(ncid, 'zc', nz_global, zc_id))

    call try(nf90_def_var(ncid, 'bx', nf90_double, (/xs_id ,yc_id, zc_id/), bx_id))
    call try(nf90_def_var(ncid, 'by', nf90_double, (/xc_id ,ys_id, zc_id/), by_id))
    call try(nf90_def_var(ncid, 'bz', nf90_double, (/xc_id ,yc_id, zs_id/), bz_id))

    call try(nf90_def_var(ncid, 'en', nf90_double, (/xc_id ,yc_id, zc_id/), en_id))
    call try(nf90_def_var(ncid, 'rho', nf90_double, (/xc_id ,yc_id, zc_id/), rho_id))

    call try(nf90_def_var(ncid, 'vx', nf90_double, (/xs_id ,ys_id, zs_id/), vx_id))
    call try(nf90_def_var(ncid, 'vy', nf90_double, (/xs_id ,ys_id, zs_id/), vy_id))
    call try(nf90_def_var(ncid, 'vz', nf90_double, (/xs_id ,ys_id, zs_id/), vz_id))

    
    call try(nf90_enddef(ncid))
    call try(nf90_close(ncid))

    end if
    call MPI_BARRIER(comm,errcode)

    !Each process writes data in turn

    do proc_write = 0 ,nproc-1
        call MPI_BARRIER(comm,errcode)

        if (rank == proc_write) then
            call try(nf90_open(trim(output_filename), nf90_write, ncid))

            call try(nf90_inq_varid(ncid, 'bx', vid))
            call try(nf90_put_var(ncid, vid, bx(0:nx,1:ny,1:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny,nz/)))

            call try(nf90_inq_varid(ncid, 'by', vid))
            call try(nf90_put_var(ncid, vid, by(1:nx,0:ny,1:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny+1,nz/)))

            call try(nf90_inq_varid(ncid, 'bz', vid))
            call try(nf90_put_var(ncid, vid, bz(1:nx,1:ny,0:nz), &
                 start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny,nz+1/)))
            
            call try(nf90_inq_varid(ncid, 'en', vid))
            call try(nf90_put_var(ncid, vid, energy(1:nx,1:ny,1:nz), &
                 start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny,nz/)))
            
            call try(nf90_inq_varid(ncid, 'rho', vid))
            call try(nf90_put_var(ncid, vid, rho(1:nx,1:ny,1:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny,nz/)))

            call try(nf90_inq_varid(ncid, 'vx', vid))
            call try(nf90_put_var(ncid, vid, vx(0:nx,0:ny,0:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny+1,nz+1/)))

            call try(nf90_inq_varid(ncid, 'vy', vid))
            call try(nf90_put_var(ncid, vid, vy(0:nx,0:ny,0:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny+1,nz+1/)))

            call try(nf90_inq_varid(ncid, 'vz', vid))
            call try(nf90_put_var(ncid, vid, vz(0:nx,0:ny,0:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny+1,nz+1/)))

            call try(nf90_close(ncid))

        end if
        call MPI_BARRIER(comm,errcode)

    end do


    call mpi_barrier(comm, errcode)
    if (rank == 0) print*, 'Saved snapshot number', snap_num, ' at time', time, 'to file ', output_filename

    return


END SUBROUTINE save_snap

subroutine try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        STOP
    end if

end subroutine try

END MODULE oliver_output







