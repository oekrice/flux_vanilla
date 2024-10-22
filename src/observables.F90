MODULE observables
  USE shared_data
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: obs_store, obs_setup, obs_close

  REAL(num), DIMENSION(:,:,:),ALLOCATABLE :: trace171, trace195, trace284, vxavg,vyavg,vzavg, rhosq
  REAL(num), DIMENSION(2) :: trace171_temp, trace195_temp, trace284_temp
  REAL(num) :: temp0, cadence, temp_min, temp_max
  REAL(num), DIMENSION(81) ::  trace_temperature = (/10000.0_num, 11220.2_num, 12589.3_num, 14125.4_num, &
       15848.9_num, 17782.8_num, &
       19952.6_num, 22387.2_num, 25118.9_num, 28183.8_num, 31622.8_num, 35481.4_num, &
       39810.7_num, 44668.4_num, 50118.7_num, 56234.1_num, 63095.8_num, 70794.6_num, &
       79432.8_num, 89125.1_num, 100000._num, 112202._num, 125893._num, 141254._num, &
       158489._num, 177828._num, 199526._num, 223872._num, 251189._num, 281838._num, &
       316228._num, 354814._num, 398107._num, 446684._num, 501187._num, 562341._num, &
       630958._num, 707946._num, 794328._num, 891251._num, 1.00000e+06_num,  1.12202e+06_num, &
       1.25893e+06_num,  1.41254e+06_num,  1.58489e+06_num,  1.77828e+06_num,  1.99526e+06_num,  2.23872e+06_num, &
       2.51189e+06_num,  2.81838e+06_num,  3.16228e+06_num,  3.54814e+06_num,  3.98108e+06_num,  4.46684e+06_num, &
       5.01187e+06_num,  5.62341e+06_num,  6.30958e+06_num,  7.07946e+06_num,  7.94328e+06_num,  8.91250e+06_num, &
       1.00000e+07_num,  1.12202e+07_num,  1.25893e+07_num,  1.41254e+07_num,  1.58489e+07_num,  1.77828e+07_num, &
       1.99526e+07_num,  2.23872e+07_num,  2.51189e+07_num,  2.81838e+07_num,  3.16228e+07_num,  3.54814e+07_num, &
       3.98108e+07_num,  4.46684e+07_num,  5.01187e+07_num,  5.62341e+07_num,  6.30958e+07_num,  7.07946e+07_num, &
       7.94328e+07_num,  8.91251e+07_num,  1.00000e+08_num/)
  REAL(num), DIMENSION(81) ::  responce_171 = (/&
       5.45892e-37_num, 2.96296e-36_num, 1.32689e-35_num, 1.40409e-34_num, 2.09738e-34_num, 6.67642e-34_num, &
       2.63632e-33_num, 1.03653e-32_num, 2.50294e-32_num, 5.91061e-32_num, 1.25398e-31_num, 2.16015e-31_num, &
       3.54893e-31_num, 7.26868e-31_num, 2.16902e-30_num, 5.84424e-30_num, 1.29420e-29_num, 2.43707e-29_num, &
       3.44712e-29_num, 3.94201e-29_num, 4.02970e-29_num, 4.11173e-29_num, 4.25963e-29_num, 4.54789e-29_num, &
       5.14511e-29_num, 6.52981e-29_num, 1.08785e-28_num, 2.25969e-28_num, 4.23163e-28_num, 7.00099e-28_num, &
       9.90517e-28_num, 1.25046e-27_num, 1.59339e-27_num, 2.14750e-27_num, 2.89122e-27_num, 3.82969e-27_num, &
       5.13403e-27_num, 6.95003e-27_num, 8.95450e-27_num, 1.06023e-26_num, 1.07646e-26_num, 8.53764e-27_num, &
       5.35122e-27_num, 2.77166e-27_num, 1.11679e-27_num, 4.33286e-28_num, 2.53220e-28_num, 1.98387e-28_num, &
       1.74851e-28_num, 1.11579e-28_num, 5.00072e-29_num, 3.58847e-29_num, 3.05370e-29_num, 2.48073e-29_num, &
       1.59738e-29_num, 1.19557e-29_num, 1.01835e-29_num, 9.67795e-30_num, 9.63824e-30_num, 9.46016e-30_num, &
       8.97847e-30_num, 8.17961e-30_num, 7.34915e-30_num, 6.76283e-30_num, 6.33254e-30_num, 5.92946e-30_num, &
       5.54267e-30_num, 5.18732e-30_num, 4.86737e-30_num, 4.58599e-30_num, 4.33676e-30_num, 4.11449e-30_num, &
       3.91329e-30_num, 3.72897e-30_num, 3.55953e-30_num, 3.40456e-30_num, 3.26201e-30_num, 3.13050e-30_num, &
       3.00831e-30_num, 2.89459e-30_num, 2.75796e-30_num /) /  1.0764600E-26_num


  REAL(num), DIMENSION(81) ::  responce_195 =(/&
       0.00000_num, 0.00000_num, 0.00000_num, 0.0_num, 0.0_num, 0.0_num, &
       0.0_num, 0.0_num, 1.77852e-37_num, 2.25485e-36_num, 2.03320e-35_num, 5.30076e-34_num, &
       4.36306e-33_num, 2.81844e-32_num, 1.46507e-31_num, 6.27947e-31_num, 2.26517e-30_num, 6.62970e-30_num, &
       1.03783e-29_num, 1.44101e-29_num, 1.57295e-29_num, 1.62683e-29_num, 1.72390e-29_num, 2.05525e-29_num, &
       3.14863e-29_num, 6.11858e-29_num, 1.38518e-28_num, 2.76311e-28_num, 3.75032e-28_num, 3.40593e-28_num, &
       2.35060e-28_num, 1.35958e-28_num, 6.84562e-29_num, 4.27816e-29_num, 4.08908e-29_num, 4.59713e-29_num, &
       6.28437e-29_num, 1.06146e-28_num, 3.25445e-28_num, 9.80586e-28_num, 2.54526e-27_num, 5.25317e-27_num, &
       7.60421e-27_num, 8.04732e-27_num, 6.73427e-27_num, 4.18724e-27_num, 1.73652e-27_num, 5.45464e-28_num, &
       1.96170e-28_num, 1.00155e-28_num, 7.42104e-29_num, 5.82088e-29_num, 4.40128e-29_num, 3.47070e-29_num, &
       2.80539e-29_num, 2.09966e-29_num, 1.45175e-29_num, 1.03065e-29_num, 8.07076e-30_num, 7.55943e-30_num, &
       1.08439e-29_num, 2.04995e-29_num, 3.41079e-29_num, 4.82806e-29_num, 5.83860e-29_num, 6.05102e-29_num, &
       5.72019e-29_num, 5.15207e-29_num, 4.51063e-29_num, 3.89888e-29_num, 3.33803e-29_num, 2.83352e-29_num, &
       2.37829e-29_num, 1.96955e-29_num, 1.60984e-29_num, 1.30576e-29_num, 1.05049e-29_num, 8.39515e-30_num, &
       6.66081e-30_num, 5.27824e-30_num, 2.08609e-30_num /) / 8.0473203E-27_num


  REAL(num), DIMENSION(81) ::  responce_284 =(/&
       0.0_num, 0.0_num, 0.0_num, 8.37951e-35_num, 9.33354e-38_num, 9.04757e-37_num, &
       9.65942e-36_num, 1.40779e-33_num, 5.03407e-34_num, 2.25085e-33_num, 1.35895e-32_num, 7.02053e-32_num, &
       1.84190e-31_num, 5.07800e-31_num, 1.75099e-30_num, 4.78710e-30_num, 1.01460e-29_num, 1.76627e-29_num, &
       2.32726e-29_num, 2.34422e-29_num, 2.17872e-29_num, 2.28581e-29_num, 2.85665e-29_num, 3.97122e-29_num, &
       5.20891e-29_num, 6.14373e-29_num, 6.74460e-29_num, 6.99926e-29_num, 6.75119e-29_num, 5.96136e-29_num, &
       5.54531e-29_num, 6.52556e-29_num, 9.12838e-29_num, 1.33515e-28_num, 1.81445e-28_num, 2.22340e-28_num, &
       2.43745e-28_num, 2.36674e-28_num, 2.14989e-28_num, 1.93995e-28_num, 1.80218e-28_num, 1.75707e-28_num, &
       1.85153e-28_num, 2.29116e-28_num, 3.56458e-28_num, 6.03612e-28_num, 7.85972e-28_num, 7.28284e-28_num, &
       5.40238e-28_num, 3.59186e-28_num, 2.20558e-28_num, 1.34941e-28_num, 8.46400e-29_num, 5.24126e-29_num, &
       3.13507e-29_num, 1.84155e-29_num, 1.06357e-29_num, 5.94984e-30_num, 3.33008e-30_num, 2.20294e-30_num, &
       1.84440e-30_num, 1.73975e-30_num, 1.71672e-30_num, 1.60171e-30_num, 1.43832e-30_num, 1.28647e-30_num, &
       1.16300e-30_num, 1.07440e-30_num, 1.00969e-30_num, 9.57979e-31_num, 9.14655e-31_num, 8.77543e-31_num, &
       8.44559e-31_num, 8.14229e-31_num, 7.85862e-31_num, 7.59100e-31_num, 7.33642e-31_num, 7.09277e-31_num, &
       6.85839e-31_num, 6.63222e-31_num, 6.41373e-31_num /) / 7.8597201E-28_num


CONTAINS


  SUBROUTINE obs_setup
    CHARACTER(LEN=20) :: logfile
    REAL(num) :: average_temp
    INTEGER :: ndims
    LOGICAL :: periods(3), reorder
    INTEGER :: starts(3), sizes(3), subsizes(3)


    ! Set up the normalised perameters.
    temp0 = 1.0e8_num

    ! As a check calculate the average temperature in the domain.
    average_temp = SUM(energy(1:nx,1:ny,1:nz) * (gamma - 1.0_num)) / REAL(nx * ny * nz,num)

    ! set up the instrument parameters in normalised units.
    cadence = 5.0_num

    ! Sort out the trace temperature responce, convert to Log \epsilon
    trace_temperature = LOG10(trace_temperature / temp0 / (gamma - 1.0_num))
    temp_min = trace_temperature(1)
    temp_max = trace_temperature(81) - temp_min

    ! set the really small reponce values to zero to stop underflow
    WHERE(responce_171 < 1.0e-4_num) responce_171 = 0.0_num 
    WHERE(responce_195 < 1.0e-4_num) responce_195 = 0.0_num 
    WHERE(responce_284 < 1.0e-4_num) responce_284 = 0.0_num 

    ! These aren't used anymore but are a handy guide so I've left them
    trace171_temp = (/4.0e5_num, 1.5e6_num/) / (gamma - 1.0_num) / temp0
    trace195_temp = (/9.0e5_num, 2.5e6_num/) / (gamma - 1.0_num) / temp0
    trace284_temp = (/1.0e6_num, 5.0e6_num/) / (gamma - 1.0_num) / temp0


    ! Write out the log 
    IF (rank ==0) THEN
       WRITE(logfile, '(a,"/observables.log")') data_dir
       OPEN(unit=60, STATUS='REPLACE',FILE=logfile)
       WRITE(60,*) 'Observables'
       WRITE(60,*) '-----------'
       WRITE(60,*)
       WRITE(60,*) 'The following shows normalised and SI'
       WRITE(60,*) 'Normalising temperature (t0) = ', temp0
       WRITE(60,*)
       WRITE(60,*) 'Initial average temp:',average_temp * temp0,'Kelvin'
       WRITE(60,*)
       WRITE(60,*) 'Trace 171 range from:',trace171_temp(1),'(',trace171_temp(1)*temp0*(gamma-1.0_num),'K)'
       WRITE(60,*) '                  to:',trace171_temp(2),'(',trace171_temp(2)*temp0*(gamma-1.0_num),'K)'
       WRITE(60,*) 'Trace 195 range from:',trace195_temp(1),'(',trace195_temp(1)*temp0*(gamma-1.0_num),'K)'
       WRITE(60,*) '                  to:',trace195_temp(2),'(',trace195_temp(2)*temp0*(gamma-1.0_num),'K)'
       WRITE(60,*) 'Trace 284 range from:',trace284_temp(1),'(',trace284_temp(1)*temp0*(gamma-1.0_num),'K)'
       WRITE(60,*) '                  to:',trace284_temp(2),'(',trace284_temp(2)*temp0*(gamma-1.0_num),'K)'

    ENDIF

    ! Set up arrays
    ALLOCATE(trace171(1:nx,1:ny,1:nz),trace195(1:nx,1:ny,1:nz),trace284(1:nx,1:ny,1:nz))
    ALLOCATE(vxavg(1:nx,1:ny,1:nz),vyavg(1:nx,1:ny,1:nz),vzavg(1:nx,1:ny,1:nz),rhosq(1:nx,1:ny,1:nz))
    trace171 = 0.0_num
    trace195 = 0.0_num
    trace284 = 0.0_num
    vxavg = 0.0_num
    vyavg = 0.0_num
    vzavg = 0.0_num

    ! Set up the MPI output subarray
    ndims = 3
    sizes = (/nx_global, ny_global, nz_global/)
    subsizes = (/nx,ny,nz/)

    starts(1) = coordinates(3) * nx
    starts(2) = coordinates(2) * ny 
    starts(3) = coordinates(1) * nz 
    
    
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, mpireal, obstype, errcode)

    CALL MPI_TYPE_COMMIT(obstype,errcode)



  END SUBROUTINE obs_setup




  SUBROUTINE obs_store
    REAL(num), SAVE :: next_snap_time = 0.0_num, previous_snap_time = 0.0_num
    INTEGER, SAVE :: snapnum = 1
    CHARACTER(LEN=17) :: filename
    LOGICAL,SAVE :: firstcall=.TRUE.
    INTEGER :: responce_number, localcellcount, filehandle
    REAL(num) :: rho_sq_dt
    LOGICAL, PARAMETER :: usempiwrite = .FALSE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! This subroutine won't restart properly

    IF (firstcall) THEN
       next_snap_time = cadence
       firstcall = .FALSE.
    ENDIF

    ! store the relevent n^2 values
    DO iz = 1, nz
       DO iy = 1, ny
          DO ix = 1, nx
             responce_number = NINT((LOG10(energy(ix,iy,iz)) - temp_min) / temp_max * 81.0_num)
             IF (responce_number > 81) responce_number = 81
             IF (responce_number < 1) responce_number = 1
             rho_sq_dt = rho(ix,iy,iz)**2 * dt
             trace171(ix,iy,iz) = trace171(ix,iy,iz) + rho_sq_dt * responce_171(responce_number)
             trace195(ix,iy,iz) = trace195(ix,iy,iz) + rho_sq_dt * responce_195(responce_number)
             trace284(ix,iy,iz) = trace284(ix,iy,iz) + rho_sq_dt * responce_284(responce_number)
             vxavg(ix,iy,iz) = vxavg(ix,iy,iz) + vx(ix,iy,iz) * rho_sq_dt
             vyavg(ix,iy,iz) = vyavg(ix,iy,iz) + vy(ix,iy,iz) * rho_sq_dt
             vzavg(ix,iy,iz) = vzavg(ix,iy,iz) + vz(ix,iy,iz) * rho_sq_dt
             rhosq(ix,iy,iz) = rhosq(ix,iy,iz) + rho_sq_dt
          END DO
       END DO
    END DO


    IF (time >= next_snap_time) THEN

       w1 = time - previous_snap_time
       trace171 = trace171 / w1
       trace195 = trace195 / w1
       trace284 = trace284 / w1
       vxavg = vxavg / rhosq
       vyavg = vyavg / rhosq
       vzavg = vzavg / rhosq
       rhosq = rhosq / w1

       IF (rank ==0) THEN
          WRITE(60,*) 'Time: ', time, 'Snapnum:',snapnum
       ENDIF

       ! Write out the data.
#ifndef NONMPIIO
       ! Create the filename
       WRITE(filename, '(a,"/",i4.4,".l3d_obs")') data_dir, snapnum

       ! Rank 0 deletes the file as there is no replace option on the file open,
       ! not deleting causes literal overwriting and potential confusion.
       IF (rank == 0) CALL MPI_FILE_DELETE(filename, MPI_INFO_NULL, errcode)

       ! Wait for the deletetion before we open the file...
       CALL MPI_BARRIER(comm,errcode)

       ! Open the file, filehandle is the MPI file unit number.
       CALL MPI_FILE_OPEN(comm, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, &
            MPI_INFO_NULL, filehandle, errcode)

       ! If rank 0 then write the file header.
       IF (rank == 0) THEN
          CALL MPI_FILE_WRITE(filehandle,(/nx_global, ny_global, nz_global /), 3, MPI_INTEGER, &
               status, errcode)
          CALL MPI_FILE_WRITE(filehandle,KIND(1.0_num),1,MPI_INTEGER,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,time,1,mpireal,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,xb_global(1:nx_global),nx_global,mpireal,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,yb_global(1:ny_global),ny_global,mpireal,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,zb_global(1:nz_global),nz_global,mpireal,status,errcode)
       ENDIF

       ! Set my view of the file
       CALL MPI_FILE_SET_VIEW(filehandle, initialdisp, mpireal, obstype,&
            "native", MPI_INFO_NULL, errcode)

       localcellcount = nx * ny * nz

       CALL MPI_FILE_WRITE_ALL(filehandle,trace171, localcellcount, mpireal, status, errcode)
       CALL MPI_FILE_WRITE_ALL(filehandle,trace195, localcellcount, mpireal, status, errcode)
       CALL MPI_FILE_WRITE_ALL(filehandle,trace284, localcellcount, mpireal, status, errcode)
       CALL MPI_FILE_WRITE_ALL(filehandle,vxavg, localcellcount, mpireal, status, errcode)
       CALL MPI_FILE_WRITE_ALL(filehandle,vyavg, localcellcount, mpireal, status, errcode)
       CALL MPI_FILE_WRITE_ALL(filehandle,vzavg, localcellcount, mpireal, status, errcode)
       CALL MPI_FILE_WRITE_ALL(filehandle,rhosq, localcellcount, mpireal, status, errcode)
       CALL MPI_FILE_CLOSE(filehandle, errcode)
#else
       WRITE(filename, '(a,"/",i3.3,i4.4,".obs")') data_dir, rank, snapnum
       OPEN(unit=61,FORM='UNFORMATTED',STATUS='REPLACE',FILE=filename)
       WRITE(61) REAL(nx_global,num)-0.1_num,REAL(ny_global,num)-0.1_num, &
            REAL(nz_global,num)-0.1_num
       WRITE(61)  REAL(nx,num)-0.1_num,REAL(ny,num)-0.1_num, REAL(nz,num)-0.1_num
       WRITE(61) time
       WRITE(61) trace171(1:nx,1:ny,1:nz)
       WRITE(61) trace195(1:nx,1:ny,1:nz)
       WRITE(61) trace284(1:nx,1:ny,1:nz)
       WRITE(61) vxavg(1:nx,1:ny,1:nz)
       WRITE(61) vyavg(1:nx,1:ny,1:nz)
       WRITE(61) vzavg(1:nx,1:ny,1:nz)
       WRITE(61) rhosq(1:nx,1:ny,1:nz)
       WRITE(61) xc(1:nx),yc(1:ny),zc(1:nz)
       CLOSE(61)
#endif

       ! Prepare for the next snapshot
       snapnum = snapnum + 1
       next_snap_time = next_snap_time + cadence
       previous_snap_time = time
       trace171 = 0.0_num
       trace195 = 0.0_num
       trace284 = 0.0_num
       vxavg = 0.0_num
       vyavg = 0.0_num
       vzavg = 0.0_num


    ENDIF

  END SUBROUTINE obs_store



  SUBROUTINE obs_close

    DEALLOCATE(trace171,trace195,trace284)
    DEALLOCATE(vxavg,vyavg,vzavg,rhosq)

    IF (rank ==0) THEN
       CLOSE(60)
    ENDIF


  END SUBROUTINE obs_close


END MODULE observables
