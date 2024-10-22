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
    REAL(num) :: r,RC,R0 !minor radius, cylindrical radius, major radius
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
    D0   = 2.5_num     ! tube radius
    B00_1  = 5.0_num   !B strength on axis
    B00_2  = 5.0_num
!    B00_3  = 5.0_num
    db   = 0.4_num     ! twist parameter
    zb0  = 8.8_num
    zb0_1 = 4.8_num 
!    zb0_1 = 5.8_num
!    zb0_3 = 2.8_num
    lamb = 20.0_num    ! density pertubation length

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
    R0=50.0_num 
    Rmin = 0.00001_num

    DO ix = -1, nx+1
       DO iy = -1, ny+1
          DO iz = -1, nz+1
             !Magnetic perturbation - Bx, By and Bz are located on faces of cells
!!toroidal flux tube
! x component at xb,yc and zc
             RC = SQRT(yc(iy)**2 + (zc(iz)+10.0_num)**2)
             r = SQRT(xb(ix)**2 + (RC - R0)**2)
             Bphi= B00_1*EXP(-r**2/D0**2)
             Btheta = db*r*Bphi
             IF (RC < Rmin)THEN
               bx(ix,iy,iz) = 0.0_num
             ELSE
               bx(ix,iy,iz)=db*Bphi*(RC-R0)
             ENDIF
! y component at xc, yb, zc
             RC = SQRT(yb(iy)**2 + (zc(iz)+10.0_num)**2)
             r = SQRT(xc(ix)**2 + (RC - R0)**2)
             Bphi= B00_1*EXP(-r**2/D0**2)
             Btheta = db*r*Bphi
             IF(RC < Rmin)THEN
               by(ix,iy,iz) = 0.0_num
             ELSE
               by(ix,iy,iz)=-Bphi*(zc(iz)+10.0_num)/RC-deb*Bphi*(xc(ix)*yb(iy))/(RC)
             ENDIF
! z component at xc, yc, zb
             RC = SQRT(yc(iy)**2 + (zb(iz)+10.0_num)**2)
             r = SQRT(xc(ix)**2 + (RC - R0)**2)
             Bphi= B00_1*EXP(-r**2/D0**2)
             Btheta = db*r*Bphi
             IF(RC < Rmin)THEN
               bz(ix,iy,iz) = 0.0_num
             ELSE
               bz(ix,iy,iz)=Bphi*yc(iy)/RC-db*Bphi*(xc(ix)*(zb(iz)+10.0_num))/(RC)
             ENDIF
!Pressure perturbation at xc, yc, zc
             RC = SQRT(yc(iy)**2 + (zc(iz)+10.0_num)**2)
             r = SQRT(xc(ix)**2 + (RC - R0)**2)
             Bphi= B00_1*EXP(-r**2/D0**2)
             Btheta = db*r*Bphi
             scr2 = B00_1**2*EXP(-2*r**2/D0**2)*((db*D0)**2-2.0_num-2.0_num*(db*r)**2)/4.0_num 
! scr2 is pressure change
             scr1 = energy(ix,iy,iz)*rho(ix,iy,iz)*(gamma-1.0_num)
!             rho(ix,iy,iz) = rho(ix,iy,iz) - rho(ix,iy,iz) * scr2 / scr1 * EXP(-yc(iy)**2 / lamb**2)

! scr1 is background pressure and now calculate new pressure
             scr1 = scr1 + scr2
! assume temperature is held constant so calculate new density
             rho(ix,iy,iz) = scr1/(energy(ix,iy,iz)*(gamma-1.0_num))
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE equilibrium



 