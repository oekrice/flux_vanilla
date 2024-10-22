MODULE initialise
  USE shared_data
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: control_variables, grid, equilibrium
  PUBLIC :: open_files, close_files

  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global, dxnew
  REAL(num), DIMENSION(:), ALLOCATABLE :: yb_global, dynew
  REAL(num), DIMENSION(:), ALLOCATABLE :: zb_global, dznew

CONTAINS


  SUBROUTINE control_variables

    x_stretch = .FALSE.      ! if false then use uniform grid
    y_stretch = .FALSE.      ! otherwise stretch grid based on stretch_x
    z_stretch = .FALSE.     

    nsteps = 10000000 ! maximum number of steps allowed 
    t_end = 300.0_num         ! final computational time
    dt_snapshots = 5.0_num

    length_x = 200.0_num 
    length_y = 200.0_num
    length_z = 150.0_num          

    tensor_shock_visc = .TRUE.
    visc1 = 0.2_num       ! linear viscosity 
    visc2 = 0.5_num       ! quadratic viscosity
    visc3 = 0.02_num       ! background shear viscosity

    resistiveMHD = .TRUE. 
    j_max = 0.0_num     ! maximum allowed current density
    eta0 = 0.0001_num        ! scales plasma resistivity
    eta_background = 0.0_num

    dt_multiplier = 0.8_num   ! fraction of CLF limit used for dt        

    gamma = 5.0_num / 3.0_num  ! ratio of specific heats
    grav = 0.0_num             ! normalised gravity

    rke = .TRUE.   ! add kinetic energy remap correction?

    curlb = 0.0_num
    p_visc = 0.0_num
    eta = 0.0_num

    data_dir = 'Data'            ! Output directory - must be 4 characters
    restart = .TRUE.            ! restart from old data?
    restart_snapshot = 34         ! snapsot in data_dir to restart from
    IF (.NOT. restart)  restart_snapshot = 0

  END SUBROUTINE control_variables



  SUBROUTINE equilibrium  


    !Fan stratified atmosphere
    REAL(num) :: z0,ztr,z1 !Position of start of {photosphere,transition region,corona}
    REAL(num) :: T0, T1 !Temperature of photosphere(isothermal), temperature of bottom of corona
    REAL(num) :: h0 !Scale height for problem
    REAL(num) :: r00 ! Density at start of photosphere
    REAL(num) :: p0 !Pressure at start of photosphere
!    REAL(num) :: D0,B00_1,B00_2,B00_3,db,zb0,zb0_1,zb0_3,lamb !Used for tube
    REAL(num) :: D0,B00_1,B00_2,db,zb0,zb0_1,lamb !Used for tube
    REAL(num) :: phBxa,phBxb,bcor !Height at which field starts to ramp up, height at which field reaches TR value, |B| in TR and corona
    REAL(num) :: bcxa,bcxb,bcshr !Field twists by bcshr*pi between heights bcxa and bcxb


    REAL(num) :: rb,scr1,scr2,scr3,scr4,scr5,scr6,A,B,pxa !Work variables currently named to match Klaus' code
    REAL(num) :: pxb,p1,r11,p2,r22,h1
    REAL(num) :: igval = 0.0_num

    REAL(num) :: dztrmin,dztrmin_global !Other work variables using my naming convention
    REAL(num) :: position = 15.0_num !x offset of paired tubes
!    REAL(num) :: sposition = 2.0_num
    INTEGER   :: ix,iy,iz !Others

    INTEGER :: flip=0

    !Shift z axis to run 0->length_z, not -length_z/2->length_z/2
    !xc=xc-length_x/2.0_num
    !xb=xb-length_x/2.0_num

    !yc=yc-length_x/2.0_num
    !yb=yb-length_x/2.0_num

    zc=zc+2.9999999999997584E-002 - 10.0_num + length_z/2.0_num
    zb=zb+2.9999999999997584E-002 - 10.0_num + length_z/2.0_num

    !Initialise variables (may put something to read in files here, but for now, can't be bothered)
    z0  = 22.0_num
    ztr = 32.0_num
    z1  = 42.0_num
    T0  = 1.0_num
    T1  = 150.0_num
    h0  = 1.0_num
    r00 = 1.0_num
    p0  = 1.0_num

    !Tubestuff here
    D0   = 2.5_num
    B00_1  = 5.0_num
    B00_2  = 5.0_num
!    B00_3  = 5.0_num
    db   = 0.4_num
    zb0  = 8.8_num
    zb0_1 = 4.8_num 
!    zb0_1 = 5.8_num
!    zb0_3 = 2.8_num
    lamb = 20.0_num

    phBxa = 27.0_num
    phBxb = 32.0_num
    bcor  = 1.0e-2_num
    bcxa  = 32.0_num
    bcxb  = 92.0_num
    bcshr = 2.0_num

    !Finally, set gravity
    grav=1.0_num

    !Want the start of the transition region to lie on a grid point, so check this
    dztrmin = MINVAL(ABS(zc-ztr))
    CALL MPI_ALLREDUCE(dztrmin, dztrmin_global, 1, mpireal, MPI_MIN, &
         comm, errcode)
    IF (dztrmin_global > 1.0e-6_num) THEN
       !If ztr is too far from grid point then warn
       IF (rank == 0) THEN 
          PRINT *,"WARNING!"
          PRINT *,"--------"
          PRINT *,'Transition region start lies too far from grid point. Please recalibrate'
          PRINT *,'Minimum dz to transition region is ',dztrmin_global
          PRINT *,'This may lead to shocks forming at the Photosphere/Transition region boundary'
          PRINT *,'Code will continue to execute'
       ENDIF
       !    STOP
    ENDIF



    !Setup constants
    A = pi / (phBxb - phBxa)
    B = pi * (-phBxa - phBxb) / (2.0_num *(phBxb - phBxa)) 
    pxa = p0 * EXP(-(phBxa - z0) / h0)
    pxb =  -1.0_num / 8.0_num * A * bcor**2 * (-2.0_num * A * COS(2.0_num * A * phBxb + 2.0_num * B) / (1.0_num / h0**2 + 4.0_num * A**2)&
         +SIN(2.0_num * A * phBxb + 2.0_num * B) / (h0 * (1.0_num / h0**2 + 4.0_num * A**2))) &
         -1.0_num/4.0_num * A * bcor**2 * (COS(A * phBxb + B) / (h0 * (1.0_num / h0**2 + A**2))&
         +A * SIN(A * phBxb + B) / (1.0_num / h0**2 + A**2))&
         +(pxa + 1.0_num / 8.0_num * A * bcor**2 * (2.0_num * A / (1.0_num / h0**2 + 4.0_num * A**2)&
         - 2.0_num * A / (1.0_num / h0**2 + A**2))) * EXP(-(phBxb - phBxa) / h0)
    p1 = pxb * EXP(-(ztr - phBxb) / h0)
    r11 = p1 / p0 * r00
    p2  = EXP((z1 - ztr) / h0 / LOG(T1 / T0) * (T0 / T1 - 1.0_num)) * p1
    r22 = (p2 / p1) / (T1 / T0) * r11
    h1 = h0 * T1 / T0

    DO ix = -1, nx+1
       DO iy = -1, ny+1
          DO iz = 0, nz+1

             !Reset work arrays
             scr1 = 0.0_num
             scr2 = 0.0_num
             scr3 = 0.0_num
             scr4 = 0.0_num
             scr5 = 0.0_num
             scr6 = 0.0_num


             !Density and pressure in convection zone (actually convectively stable, but what the hell)
             IF (zc(iz) .LE. z0) THEN
                rb = (1.0_num - (zc(iz) - z0) / h0 * (gamma - 1.0_num) / gamma)**(1.0_num / (gamma - 1.0_num)) * r00
                scr1 = (1.0_num - (zc(iz) - z0) / h0 * (gamma - 1.0_num) / gamma)**(gamma / (gamma - 1.0_num)) * p0
             ENDIF

             !scr6 is |B|, here ramping up to TR value. Use SIN function since == 0 at finite value and 0 gradient at end
             IF (zc(iz) >= phBxa .AND. zc(iz) <= phBxb) THEN
                scr6 = bcor * (SIN( A * zc(iz) + B) + 1.0_num)/2.0_num
             ENDIF
             !Make sure boundary values are exactly what is expected
             IF (zc(iz) <= z0 .AND. zc(iz+1) > z0) THEN
                scr6 = 0.0_num
             ENDIF
             IF (zc(iz-1) <ztr .AND. zc(iz) >= z0) THEN
                scr6=bcor
             ENDIF

             !If in photosphere, but before magnetic field turns on
             !then have simple exponential pressure stratification
             !with scaleheight h0
             !Isothermal, so same density scaling
             IF (zc(iz) > z0 .AND. zc(iz) <=phBxa) THEN
                scr1 = p0 * EXP(- (zc(iz) - z0) / h0)
                rb = scr1 / p0 * r00
             ENDIF

             !If in photosphere after magnetic field turns on
             !Once again, isothermal assumed
             IF (zc(iz) > phBxa .AND. zc(iz) <= phBxb) THEN
                scr1 = -1.0_num/8.0_num * A * bcor**2 * (-2.0_num * A * COS(2.0_num * A * zc(iz) + 2.0_num * B) / (1.0_num / h0**2 + 4.0_num * A**2)&
                     +SIN(2.0_num * A * zc(iz) + 2.0_num * B) / (h0 *(1.0_num / h0**2 + 4.0_num * A**2)))&
                     -1.0_num / 4.0_num * A * bcor**2 * (COS(A * zc(iz) + B) / (h0 * (1.0_num / h0**2 + A**2))&
                     +A * SIN( A * zc(iz) + B) / (1.0_num / h0**2 + A**2))&
                     + (pxa + 1.0_num / 8.0_num * A * bcor**2 * (2.0_num * A / (1.0_num / h0**2 + 4.0_num * A**2)-2.0_num &
                     * A / (1.0_num / h0**2 + A**2))) * EXP(-(zc(iz) - phBxa) / h0)

                rb = scr1 / p0 * r00
             ENDIF

             !If the field has reached maximum before entering the transition region, then field is force free, go back to simple
             !exponential isothermal stratification until reach transition region
             IF (zc(iz) > phBxb .AND. zc(iz) <= ztr) THEN
                scr1 = pxb * EXP(-(zc(iz) - phBxb) / h0)
                rb = scr1 / p0 * r00
             ENDIF

             !Transition region between ztr and z1.
             IF (zc(iz) > ztr .AND. zc(iz) <= z1) THEN
                scr1 = EXP((z1 - ztr) / h0 / LOG(T1 / T0) * &
                     ((T1 / T0)**(-(zc(iz) - ztr) / (z1-ztr)) - 1.0_num)) * p1
                rb   = (scr1 / p1) / ((T1 / T0)**((zc(iz) - ztr) / (z1 - ztr))) * r11
             ENDIF

             !Finally the corona, again isothermal
             IF (zc(iz) > z1) THEN
                rb   = r22 * exp(-(zc(iz) - z1) / h1)
                scr1 = p2 * exp(-(zc(iz) - z1) / h1)
             ENDIF

             !Calculate ambient magnetic field
             IF (zc(iz) > z0 .AND. zc(iz) <= ztr) THEN
                !     uni-directional coronal field 
                !scr3 = bcxa  * scr6
                !scr4 = bcxb  * scr6
                !scr5 = bcshr * scr6
                !     uniformally sheared coronal field
                scr3 = 0.0_num * scr6
                scr4 = 1.0_num * SIN(pi * bcshr * (bcxa - ztr) / (bcxa - bcxb)) * scr6
                scr5 = 1.0_num * COS(pi * bcshr * (bcxa - ztr) / (bcxa - bcxb)) * scr6
             ENDIF
             !             IF (zc(iz) > ztr) THEN
             !     uni-directional coronal field 
             !scr3 = bcxa  * bcor
             !scr4 = bcxb  * bcor
             !scr5 = bcshr * bcor
             !     uniformally sheared coronal field
             scr3 = 0.00_num !* SQRT(2.0_num)
             scr4 = -0.0_num * SQRT(2.0_num)!0.0_num * SIN(pi * bcshr * (bcxa - zc(iz)) / (bcxa - bcxb)) * bcor
             scr5 = 0.0_num * SQRT(2.0_num)!0.0_num * COS(pi * bcshr * (bcxa - zc(iz)) / (bcxa - bcxb)) * bcor
             !             ENDIF

             !Now actually assign these to primitive variables
             rho(ix,iy,iz) = rb
             energy(ix,iy,iz) = scr1 / (rho(ix,iy,iz) * (gamma - 1.0_num))
             bz(ix,iy,iz) = scr5
             bx(ix,iy,iz) = scr4
             by(ix,iy,iz) = scr3


             !Finally, set the velocity mask
!!$             IF (zc(iz) >= 100.0_num) THEN
!!$                alpha(ix,iy,iz)= (alpha(ix,iy,iz)+(1.0_num-SIN((ABS(zc(iz))-100.0_num)*pi/2.0_num/16.0_num)**2))/2.0_num
!!$             ENDIF
!!$
!!$             IF (ABS(xc(ix)) >= 90.0_num) THEN
!!$                alpha(ix,iy,iz) = (alpha(ix,iy,iz)+(1.0_num-SIN((ABS(xc(ix))-90.0_num)*pi/2.0_num/40.0_num)**2))/2.0_num
!!$             ENDIF
!!$
!!$             IF (ABS(yc(iy)) >= 90.0_num) THEN
!!$                alpha(ix,iy,iz) = (alpha(ix,iy,iz)+(1.0_num-SIN((ABS(yc(iy))-90.0_num)*pi/2.0_num/40.0_num)**2))/2.0_num
!!$             ENDIF

          ENDDO
       ENDDO
    ENDDO

    !Now introduce tube

    p0 = (1.0_num - (zb0 - z0) / h0 * (gamma - 1.0_num) / gamma)**(gamma / (gamma - 1.0_num)) * p0

    bx=0
    by=0
    bz=0

    DO ix = -1, nx+1
       DO iy = -1, ny+1
          DO iz = -1, nz+1
             !Magnetic perturbation
!! second flux tube
             scr2 = (zc(iz) - zb0_1)**2 + (xc(ix)-position)**2
             scr3 = B00_1 * EXP(-scr2 / D0**2) 
             by(ix,iy,iz) = scr3 + by(ix,iy,iz)

             scr2 = (zc(iz) - zb0_1)**2 + (xc(ix)-position)**2
             scr3 =B00_1 * EXP(-scr2 / D0**2)  * db
             bz(ix,iy,iz) = -(xc(ix)-position) * scr3 + bz(ix,iy,iz)

             scr2 = (zc(iz) - zb0_1)**2 + (xc(ix)-position)**2
             scr3 =B00_1 * EXP(-scr2 / D0**2)  * db
             bx(ix,iy,iz) = (zc(iz) - zb0_1) * scr3 + bx(ix,iy,iz)

             !Pressure perturbation
             scr2 = (zc(iz) - zb0_1)**2 + (xc(ix)-position)**2
             scr3 = B00_1**2.0_num / 2.0_num * EXP(-2.0_num * scr2 / D0**2) * &
                  (1.0_num + scr2 * db**2 - db**2 * D0**2 / 2.0_num)
             scr1=energy(ix,iy,iz)*rho(ix,iy,iz)*(gamma-1.0_num)
             rho(ix,iy,iz) = rho(ix,iy,iz) - rho(ix,iy,iz) * scr3 / scr1 * EXP(-yc(iy)**2 / lamb**2)
             scr1 = scr1 - scr3

!!$
!! first flux tube
             scr2 = (zc(iz) - zb0)**2 + (xc(ix)+position)**2
             scr3 = B00_2 * EXP(-scr2 / D0**2)
             by(ix,iy,iz) = scr3 + by(ix,iy,iz)

             scr2 = (zc(iz) - zb0)**2 + (xc(ix)+position)**2
             scr3 =B00_2 * EXP(-scr2 / D0**2)  * db
             bz(ix,iy,iz) = -(xc(ix)+position) * scr3 + bz(ix,iy,iz)

             scr2 = (zc(iz) - zb0)**2 + (xc(ix)+position)**2
             scr3 =B00_2 * EXP(-scr2 / D0**2)  * db
             bx(ix,iy,iz) = (zc(iz) - zb0) * scr3 + bx(ix,iy,iz)

             !Pressure perturbation
             scr2 = (zc(iz) - zb0)**2 + (xc(ix)+position)**2
             scr3 = B00_2**2.0_num / 2.0_num * EXP(-2.0_num * scr2 / D0**2) * &
                  (1.0_num + scr2 * db**2 - db**2 * D0**2 / 2.0_num)
             !scr1=energy(ix,iy,iz)*rho(ix,iy,iz)*(gamma-1.0_num)
             rho(ix,iy,iz) = rho(ix,iy,iz) - rho(ix,iy,iz) * scr3 / scr1 * EXP(-yc(iy)**2 / lamb**2)
             scr1 = scr1 - scr3

!! third flux tube
!             scr2 = (zc(iz) - zb0_3)**2 + (xc(ix)+sposition)**2
!             scr3 = B00_3 * EXP(-scr2 / D0**2)
!             by(ix,iy,iz) = scr3 + by(ix,iy,iz)
!
!             scr2 = (zc(iz) - zb0_3)**2 + (xc(ix)+sposition)**2
!             scr3 =B00_3 * EXP(-scr2 / D0**2)  * db
!             bz(ix,iy,iz) = -(xc(ix)+sposition) * scr3 + bz(ix,iy,iz)
!
!             scr2 = (zc(iz) - zb0_3)**2 + (xc(ix)+sposition)**2
!             scr3 =B00_3 * EXP(-scr2 / D0**2)  * db
!             bx(ix,iy,iz) = (zc(iz) - zb0_3) * scr3 + bx(ix,iy,iz)

             !Pressure perturbation
!             scr2 = (zc(iz) - zb0_3)**2 + (xc(ix)+sposition)**2
!             scr3 = B00_3**2.0_num / 2.0_num * EXP(-2.0_num * scr2 / D0**2) * &
!                  (1.0_num + scr2 * db**2 - db**2 * D0**2 / 2.0_num)
!             rho(ix,iy,iz) = rho(ix,iy,iz) - rho(ix,iy,iz) * scr3 / scr1 * EXP(-yc(iy)**2 / lamb**2)
!             scr1 = scr1 - scr3

             energy(ix,iy,iz) = scr1 / (rho(ix,iy,iz) * (gamma - 1.0_num))
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE equilibrium



  SUBROUTINE grid                 ! stretched and staggered grid

    REAL(num) :: dx, dy, dz, dxmin, dzmin, dymin, xcstar, ycstar, zcstar
    INTEGER :: ix, iy

    ALLOCATE(xb_global(-2:nx_global+2), dxnew(-2:nx_global+2))
    ALLOCATE(yb_global(-2:ny_global+2), dynew(-2:ny_global+2))
    ALLOCATE(zb_global(-2:nz_global+2), dznew(-2:nz_global+2))

    dx = 1.0_num / REAL(nx_global)       ! initially assume uniform grid
    dy = 1.0_num / REAL(ny_global)
    dz = 1.0_num / REAL(nz_global)

    xb_global(0) = -0.5_num               !grid cell boundary for x coordinates
    DO ix = -2, nx_global+2
       xb_global(ix) =  xb_global(0) + REAL(ix, num)*dx
    END DO
    xb_global = xb_global * length_x
    IF (x_stretch) CALL stretch_x     ! stretch grid ?
    !define position of ghost cells using sizes of adjacent cells
    xb_global(nx_global+1) = xb_global(nx_global) + &
         (xb_global(nx_global) - xb_global(nx_global-1)) 
    ! needed for ghost cell
    xb_global(nx_global+2) = xb_global(nx_global+1) &
         + (xb_global(nx_global+1) - xb_global(nx_global))
    xb_global(-1) = xb_global(0) &
         - (xb_global(1) - xb_global(0))
    xb_global(-2) = xb_global(-1) &
         - (xb_global(0) - xb_global(-1))
    xb = xb_global(coordinates(3)*nx-2:coordinates(3)*nx+nx+2)

    DO ix = -1, nx+2
       ixm = ix - 1
       xc(ix) = 0.5_num*(xb(ixm) + xb(ix))     ! cell centre
    END DO
    DO ix = -1, nx+1
       ixp = ix + 1
       dxc(ix) = xc(ixp) - xc(ix)    ! distance between centres
    END DO
    IF (coordinates(3)==nprocx-1) THEN
       dxc(nx+2) = dxc(nx+1)
    ELSE
       xcstar = 0.5_num*(xb(nx+2) + xb_global(coordinates(3)*nx+nx+3))
       dxc(nx+2) = xcstar - xc(nx+2)
    END IF
    DO ix = -1, nx+2
       ixm = ix - 1
       dxb(ix) = xb(ix) - xb(ixm)    ! cell width
    END DO

    yb_global(0) = -0.5_num               !grid cell boundary for y coordinates
    DO iy = -2, ny_global+2
       yb_global(iy) =  yb_global(0) + REAL(iy, num)*dy
    END DO
    yb_global = yb_global * length_y
    IF (y_stretch) CALL stretch_y     ! stretch grid ?
    !define position of ghost cells using sizes of adjacent cells
    yb_global(ny_global+1) = yb_global(ny_global) + &
         (yb_global(ny_global) - yb_global(ny_global-1)) 
    ! needed for ghost cell
    yb_global(ny_global+2) = yb_global(ny_global+1) &
         + (yb_global(ny_global+1) - yb_global(ny_global))
    yb_global(-1) = yb_global(0) &
         - (yb_global(1) - yb_global(0))
    yb_global(-2) = yb_global(-1) &
         - (yb_global(0) - yb_global(-1))
    yb = yb_global(coordinates(2)*ny-2:coordinates(2)*ny+ny+2)
    DO iy = -1, ny+2
       iym = iy - 1
       yc(iy) = 0.5_num*(yb(iym) + yb(iy))     ! cell centre
    END DO
    DO iy = -1, ny+1
       iyp = iy + 1
       dyc(iy) = yc(iyp) - yc(iy)    ! distance between centres
    END DO
    IF (coordinates(2)==nprocy-1) THEN
       dyc(ny+2) = dyc(ny+1)
    ELSE
       ycstar = 0.5_num*(yb(ny+2) + yb_global(coordinates(2)*ny+ny+3))
       dyc(ny+2) = ycstar - yc(ny+2)
    END IF
    DO iy = -1, ny+2
       iym = iy - 1
       dyb(iy) = yb(iy) - yb(iym)    ! cell width
    END DO

    zb_global(0) = -0.5_num               !grid cell boundary for z coordinates
    DO iz = -2, nz_global+2
       zb_global(iz) =  zb_global(0) + REAL(iz, num)*dz
    END DO
    zb_global = zb_global * length_z
    IF (z_stretch) CALL stretch_z     ! stretch grid ?
    !define position of ghost cells using sizes of adjacent cells
    zb_global(nz_global+1) = zb_global(nz_global) + &
         (zb_global(nz_global) - zb_global(nz_global-1)) 
    ! needed for ghost cell
    zb_global(nz_global+2) = zb_global(nz_global+1) &
         + (zb_global(nz_global+1) - zb_global(nz_global))
    zb_global(-1) = zb_global(0) &
         - (zb_global(1) - zb_global(0))
    zb_global(-2) = zb_global(-1) &
         - (zb_global(0) - zb_global(-1))
    zb = zb_global(coordinates(1)*nz-2:coordinates(1)*nz+nz+2)
    DO iz = -1, nz+2
       izm = iz - 1
       zc(iz) = 0.5_num*(zb(izm) + zb(iz))     ! cell centre
    END DO
    DO iz = -1, nz+1
       izp = iz + 1
       dzc(iz) = zc(izp) - zc(iz)    ! distance between centres
    END DO
    IF (coordinates(1)==nprocz-1) THEN
       dzc(nz+2) = dzc(nz+1)
    ELSE
       zcstar = 0.5_num*(zb(nz+2) + zb_global(coordinates(1)*nz+nz+3))
       dzc(nz+2) = zcstar - zc(nz+2)
    END IF
    DO iz = -1, nz+2
       izm = iz - 1
       dzb(iz) = zb(iz) - zb(izm)    ! cell width
    END DO

    DO ix = -1, nx+2
       DO iy = -1, ny+2
          DO iz = -1, nz+2
             cv(ix,iy,iz) = dxb(ix) * dyb(iy) * dzb(iz)    ! define the cell area
          END DO
       END DO
    END DO

    CALL MPI_ALLREDUCE(MINVAL(dxb(1:nx)), dxmin, 1, mpireal, MPI_MIN, &
         comm, errcode)
    CALL MPI_ALLREDUCE(MINVAL(dyb(1:ny)), dymin, 1, mpireal, MPI_MIN, &
         comm, errcode)
    CALL MPI_ALLREDUCE(MINVAL(dzb(1:nz)), dzmin, 1, mpireal, MPI_MIN, &
         comm, errcode)
    min_grid_spacing = MIN(dxmin, dymin, dzmin)

    DEALLOCATE(xb_global, dxnew)
    DEALLOCATE(yb_global, dynew)
    DEALLOCATE(zb_global, dznew)

  END SUBROUTINE grid



  SUBROUTINE stretch_x   ! replace with any stretching algorithm as needed

    REAL(num) :: width, dx, L, f, lx_new

    lx_new = 200.0_num                ! new tolal length
    L = length_x / 2.0_num       ! centre of tanh stretching in unstretched coordinates
    width = length_x / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (lx_new - length_x)/(length_x - L)/2.0_num

    dx = length_x / REAL(nx_global,num)  
    dxnew = dx + f*(1.0_num+TANH((ABS(xb_global)-L)/width))*dx
    
!!$    DO ix = nx_global/2+1, nx_global+2
!!$       xb_global(ix) = xb_global(ix-1) + dxnew(ix)
!!$    ENDDO
!!$    DO ix = nx_global/2-1, -2, -1
!!$       xb_global(ix) = xb_global(ix+1) - dxnew(ix)
!!$    ENDDO

    DO ix = 1, nx_global+2
       xb_global(ix) = xb_global(ix-1) + dxnew(ix)
    ENDDO
    
  END SUBROUTINE stretch_x



  SUBROUTINE stretch_y 

    REAL(num) :: width, dy, L, f, ly_new

    ly_new = 33.0_num                ! new tolal length
    L = 2.0_num * length_y / 3.0_num       ! centre of tanh stretching in unstretched coordinates
    width = length_y / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (ly_new - length_y)/(length_y - L)/2.0_num

    dy = length_y / REAL(ny_global,num)  
    dynew = dy + f*(1.0_num+TANH((ABS(yb_global)-L)/width))*dy

!!$    DO iy = ny_global/2+1, ny_global+2
!!$       yb_global(iy) = yb_global(iy-1) + dynew(iy)
!!$    ENDDO
!!$    DO iy = ny_global/2-1, -2, -1
!!$       yb_global(iy) = yb_global(iy+1) - dynew(iy)
!!$    ENDDO

    DO iy = 1, ny_global+2
       yb_global(iy) = yb_global(iy-1) + dynew(iy)
    ENDDO
    
  END SUBROUTINE stretch_y




  SUBROUTINE stretch_z 

    REAL(num) :: width, dz, L, f, lz_new

    lz_new = 33.0_num                ! new tolal length
    L = 2.0_num * length_z / 3.0_num       ! centre of tanh stretching in unstretched coordinates
    width = length_z / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (lz_new - length_z)/(length_z - L)/2.0_num

    dz = length_z / REAL(nz_global,num)  
    dznew = dz + f*(1.0_num+TANH((ABS(zb_global)-L)/width))*dz

!!$    DO iz = nz_global/2+1, nz_global+2
!!$       zb_global(iz) = zb_global(iz-1) + dznew(iz)
!!$    ENDDO
!!$    DO iz = nz_global/2-1, -2, -1
!!$       zb_global(iz) = zb_global(iz+1) - dznew(iz)
!!$    ENDDO

    DO iz = 1, nz_global+2
       zb_global(iz) = zb_global(iz-1) + dznew(iz)
    ENDDO
    
  END SUBROUTINE stretch_z


  SUBROUTINE open_files

    CHARACTER(LEN=15) :: file2
    CHARACTER(LEN=11) :: file3

    IF (rank == 0) THEN
       WRITE(file2, '(a,"/lare3d.dat")') data_dir
       OPEN(unit=20, STATUS = 'REPLACE',FILE = file2)
       WRITE(file3, '(a,"/en.dat")') data_dir
       OPEN(unit=30, STATUS = 'REPLACE',FILE = file3)
    END IF

  END SUBROUTINE open_files



  SUBROUTINE close_files

    IF (rank == 0) THEN
       CLOSE(unit=20)
       CLOSE(unit=30)
    END IF

  END SUBROUTINE close_files



  SUBROUTINE restart_data

    CHARACTER(LEN=16) :: filename

    WRITE(filename, '(a,"/",i3.3,i4.4,".dat")') data_dir, rank, restart_snapshot

    OPEN(unit=50, FORM = 'UNFORMATTED',STATUS = 'OLD', &
         ACTION='READ',FILE = filename)
    READ(50) time, time, time
    READ(50) time, time, time
    READ(50) time
    READ(50) rho(1:nx,1:ny,1:nz), energy(1:nx,1:ny,1:nz)
    READ(50) vx(1:nx,1:ny,1:nz), vy(1:nx,1:ny,1:nz), vz(1:nx,1:ny,1:nz)
    READ(50) bx(1:nx,1:ny,1:nz), by(1:nx,1:ny,1:nz), bz(1:nx,1:ny,1:nz)
    READ(50) xc(1:nx), yc(1:ny), zc(1:nz)
    READ(50) xb(1:nx+1), yb(1:ny+1), zb(1:nz+1)


    CLOSE(50)

  END SUBROUTINE restart_data


END MODULE initialise

