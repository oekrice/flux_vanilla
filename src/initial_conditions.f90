MODULE initial_conditions

  USE shared_data
  USE strings

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: HandleCustomBlock,Equilibrium,flux_rope


CONTAINS

  !-----------------------------------------------------------------------------
  !This function contains the equilibrium
  !-----------------------------------------------------------------------------

  SUBROUTINE Equilibrium
    INTEGER :: ix, iy, iz
    REAL(num) :: photo, trans, cor, base, dg, a, b0, alpha, bphi, btheta, RC, r, pext, pb, Rmin, R0, rho_ph, delta, x0, z0, lam, b1, r2, p, d1, d2, b00, w00, pi, backfield_strength

! Set velocities and magnetic fields to zero

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num


    grav = 1.0_num

! Parameters

    photo = 0.0_num          ! base of photosphere
    trans = 10.0_num         ! base of transition region
    cor = 20.0_num           ! base of corona
    base = -25.0_num         ! base of numerical box
    a = 3.0_num              ! tube radius (minor radius)
    R0 = 15.0_num            ! major radius
    b0 = 6.0_num             ! field strength at axis
    alpha = 0.4_num          ! twist parameter
    Rmin = 0.00001_num       ! minimum radius (to stop dividing by zero)
    delta = 1.0_num          ! super-adiabatic parameter
    
    pi = 3.14159265_num

    DO ix = -1, nx+1
       DO iy = -1, ny+1
          DO iz = -1, nz+1
             
             ! Energy profile

             IF (zc(iz) .LE. photo) THEN
                energy(ix,iy,iz) = 1.0_num - delta*zc(iz)*(gamma-1.0_num)/gamma 
             ELSEIF (zc(iz) .LE. trans) THEN
                energy(ix,iy,iz) = 1.0_num
             ELSEIF (zc(iz) .LE. cor) THEN
                energy(ix,iy,iz) = (150.0_num)**((zc(iz)-trans)/(cor-trans))
             ELSE
                energy(ix,iy,iz) = 150.0_num
             ENDIF

          ENDDO
       ENDDO
    ENDDO

    energy = energy/(gamma-1.0_num)

!!$! Solve the hydrostatic equation to get rho
!!$
!!$    ! Get density at the base of the box
!!$
    rho = (1.0_num - (base-photo)*delta*(gamma-1.0_num)/gamma)**((gamma*(1.0_num-delta)+delta)/(delta*(gamma-1.0_num)))
!!$

    IF (coordinates(1)/=0) THEN
      CALL MPI_RECV(rho(:,:,-1),(nx+4)*(ny+4),mpireal,front,tag,comm,status,errcode)
      CALL MPI_RECV(energy(:,:,-1),(nx+4)*(ny+4),mpireal,front,tag,comm,status,errcode)
    ENDIF

    DO ix = -1,nx+2
       DO iy = -1,ny+2
          DO iz = 0,nz+2
             dg = 1.0_num/(dzb(iz)+dzb(iz-1))
             rho(ix,iy,iz) = rho(ix,iy,iz-1)*(energy(ix,iy,iz-1)/dzc(iz-1)*(gamma-1.0_num) &
                  - grav(iz-1)*dzb(iz-1)*dg)
             rho(ix,iy,iz) = rho(ix,iy,iz)/(energy(ix,iy,iz)/dzc(iz-1)*(gamma-1.0_num) &
                  + grav(iz-1)*dzb(iz)*dg)
          ENDDO
       ENDDO
    ENDDO

    IF (coordinates(1)/=nprocz-1) THEN
      CALL MPI_SEND(rho(:,:,nz-1),(nx+4)*(ny+4),mpireal,back,tag,comm,errcode)
      CALL MPI_SEND(energy(:,:,nz-1),(nx+4)*(ny+4),mpireal,back,tag,comm,errcode)
    ENDIF

  !Add flux rope
  x0 = 0.0_num       ! Tube positions
  z0 = -10.0_num
  b0 = 5.0_num       ! Field strengths
  alpha = 0.4_num    ! Twists
  a = 2.5_num        ! Radius
  lam = 20.0_num     ! Middle curvature

  DO ix = -1, nx+1
     DO iy = -1, ny+1
        DO iz = -1, nz+1

           ! First tube
           r2 = (xc(ix)-x0)**2 + (zc(iz)-z0)**2
           b1 = b0*EXP(-r2/a**2)
           by(ix,iy,iz) = b1 + by(ix,iy,iz)
          
           
           r2 = (xb(ix)-x0)**2 + (zc(iz)-z0)**2
           b1 = b0*EXP(-r2/a**2)
           bx(ix,iy,iz) = -b1*alpha*(zc(iz)-z0) + bx(ix,iy,iz)

           r2 = (xc(ix)-x0)**2 + (zb(iz)-z0)**2
           b1 = b0*EXP(-r2/a**2)
           bz(ix,iy,iz) = b1*alpha*(xc(ix)-x0) + bz(ix,iy,iz)

           r2 = (xc(ix)-x0)**2 + (zc(iz)-z0)**2

           b1 = b0**2/2.0*EXP(-2.0_num*r2/a**2) * &    ! Not field but pexc.
                (1.0_num + r2*alpha**2 - alpha**2*a**2/2.0_num)
           p = energy(ix,iy,iz)*rho(ix,iy,iz)*(gamma-1.0_num)
           rho(ix,iy,iz) = rho(ix,iy,iz) - rho(ix,iy,iz) * &
                b1/p*EXP(-yc(iy)**2/lam**2)    
           p = p - b1
           energy(ix,iy,iz) = p/(rho(ix,iy,iz)*(gamma-1.0_num))
         END DO
     END DO
  END DO  

  !Add background magnetic field.
  !Zero at the top of the photosphere, rising linearly to backfield_strength above


  backfield_strength = 0.005_num
  DO ix = -1, nx+1
     DO iy = -1, ny+1
        DO iz = -1, nz+1
          bx(ix, iy,iz) = bx(ix,iy,iz) + backfield_strength
        END DO
     END DO
  END DO

  energy(:,:,nz+1) = energy(:,:,nz)
  rho(:,:,nz+1) = rho(:,:,nz+1)

  END SUBROUTINE Equilibrium


  SUBROUTINE flux_rope

    INTEGER:: ix, iy, iz
    REAL :: z0, b0, alpha, a, lam, delta, p, b1, r1, r2, RC, bphi, btheta, Rmin, base, pb, d1, d2, R0, R,pext 

    base = -25.0_num         ! base of numerical box
    a = 2.5_num              ! tube radius (minor radius)
    R0 = 15.0_num            ! major radius
    b0 = -5.0_num             ! field strength at axis
    alpha = 0.4_num          ! twist parameter
    Rmin = 0.00001_num       ! minimum radius (to stop dividing by zero)
    delta = 1.0_num          ! super-adiabatic parameter
    d1 = 0.0_num            ! distance of first tube
    d2 = -30.0_num           ! position of second tube

    DO ix = -1,nx+1
       DO iy = -1,ny+1
          DO iz = -1,nz+1
             
             ! x-component
             RC = SQRT(yc(iy)**2 + (zc(iz)-base)**2)
             r = SQRT((xb(ix)-d2)**2 + (RC-R0)**2)
             bphi = b0*exp(-r**2/a**2)
             btheta = alpha*r*bphi
             IF (RC < Rmin) THEN
                bx(ix,iy,iz) = 0.0_num + bx(ix,iy,iz)
             ELSE
                bx(ix,iy,iz) = alpha*bphi*(RC-R0) + bx(ix,iy,iz)
             ENDIF

             ! y-component
             RC = SQRT(yb(iy)**2 + (zc(iz)-base)**2)
             r = SQRT((xc(ix)-d2)**2 + (RC-R0)**2)
             bphi = b0*exp(-r**2/a**2)
             btheta = alpha*r*bphi
             IF (RC < Rmin) THEN
                by(ix,iy,iz) = 0.0_num + by(ix,iy,iz)
             ELSE
                by(ix,iy,iz) = by(ix,iy,iz)-bphi*(zc(iz)-base)/RC-alpha*bphi*((xc(ix)-d2)*yb(iy))/RC
             ENDIF

             ! z-component
             RC = SQRT(yc(iy)**2 + (zb(iz)-base)**2)
             r = SQRT((xc(ix)-d2)**2 + (RC-R0)**2)
             bphi = b0*exp(-r**2/a**2)
             btheta = alpha*r*bphi
             IF (RC < Rmin) THEN
                bz(ix,iy,iz) = 0.0_num + bz(ix,iy,iz)
             ELSE
                bz(ix,iy,iz) = bz(ix,iy,iz) + bphi*yc(iy)/RC-alpha*bphi*((xc(ix)-d2)*(zb(iz)-base))/RC
             ENDIF

             ! pressure perturbation
             RC = SQRT(yc(iy)**2 + (zc(iz)-base)**2)
             r = SQRT((xc(ix)-d2)**2 + (RC-R0)**2)
             bphi = b0*exp(-r**2/a**2)
             btheta = alpha*r*bphi
             pext = b0**2*exp(-2.0_num*r**2/a**2)*((alpha*a)**2-2.0_num-2.0_num*(alpha*r)**2)/4.0_num
             pb = energy(ix,iy,iz)*rho(ix,iy,iz)*(gamma-1.0_num)
             pb = pb + pext
             rho(ix,iy,iz) = pb/(energy(ix,iy,iz)*(gamma-1.0_num))
             energy(ix,iy,iz) = pb/(rho(ix,iy,iz)*(gamma-1.0_num))
          ENDDO
       ENDDO
    ENDDO

             
    energy(:,:,nz+1) = energy(:,:,nz)
    energy(:,:,nz+2) = energy(:,:,nz)
    rho(:,:,nz+1) = rho(:,:,nz)
    rho(:,:,nz+2) = rho(:,:,nz)

   END SUBROUTINE flux_rope

  !-----------------------------------------------------------------------------
  !These functions contain the user input deck elements
  !-----------------------------------------------------------------------------

  FUNCTION HandleCustomBlock(blockname,Element,Value)

    CHARACTER(len=30),INTENT(IN)::blockname,Element,Value
    INTEGER :: HandleCustomBlock
    LOGICAL :: Result

    !The following line must always be present
    HandleCustomBlock=ERR_UNKNOWN_BLOCK
  
  END FUNCTION HandleCustomBlock


END MODULE initial_conditions
