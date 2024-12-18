!Kink unstable loop from Arber et al, 1999

  SUBROUTINE equilibrium  

    INTEGER :: ix, iy, iz
    REAL(num) :: rc, x1, y1, b_theta, amp, k, r0, a, r1, mu
    REAL(num):: b0, bz0

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

    bx = 1.0_num
    by = 0.0_num
    bz = 0.0_num
    
    r0 = 1.0_num / SQRT(6.0_num)
    a = 5.0_num / (6.0_num * SQRT(6.0_num))

    bz0 = 1.0_num
    b0 = 4.3_num

    DO ix = -1, nx+2             ! setup static equilibrium values
       DO iy = -1, ny+2
          DO iz = -1, nz+2

             rc = SQRT(xc(ix)**2 + yc(iy)**2)
             IF (rc >= 1.0_num) rc = 1.0_num
             bz(ix,iy,iz) = (1.0/2.0)*rc**2 - (3.0/8.0)*(rc**4/r0**2)    &
                  + (7.0*a/25.0)*(rc**5/r0**3) + (1.0/12.0)*(rc**6/r0**4)     &
                  - (9.0*a/70.0)*(rc**7/r0**5) + (1.0*a**2/20.0)*(rc**8/r0**6)
             bz(ix,iy,iz) = SQRT(bz0**2 - b0**2*bz(ix,iy,iz))
             rho(ix,iy,iz) = 0.45_num*(1.0_num+COS(pi*rc))+0.1_num

             x1 = xb(ix) 
             y1 = yc(iy) 
             rc = SQRT(x1**2 + y1**2)
             IF (rc >= 1.0_num) rc = 1.0_num
             b_theta = rc/2.0 - rc**3/(4.0*r0**2) + a*rc**4/(5.0*r0**3)
             bx(ix,iy,iz) = - b0 * b_theta * y1 / rc  

             x1 = xc(ix) 
             y1 = yb(iy) 
             rc = SQRT(x1**2 + y1**2)
             IF (rc >= 1.0_num) rc = 1.0_num
             b_theta = rc/2.0 - rc**3/(4.0*r0**2) + a*rc**4/(5.0*r0**3)
             by(ix,iy,iz) = b0 * b_theta * x1 / rc

          END DO
       END DO
    END DO
    bx(-2,:,:) = bx(-1,:,:)
    by(:,-2,:) = by(:,-1,:)
    bz(:,:,-2) = bz(:,:,-1)


    WHERE (rho < 0.1_num) rho = 0.1_num
    energy = 0.01_num / (rho * (gamma-1.0_num))

    k = 2.0_num * pi / length_z     ! apply velocity perturbation 
    amp = 1.e-2_num
    r1 = 0.95_num
    mu = 0.2_num

    DO ix = -1, nx+2
       DO iy = -1, ny+2
          DO iz = -1, nz+2
             rc = SQRT(xb(ix)**2 + yb(iy)**2)
             IF (rc < r1) THEN
                vx(ix,iy,iz) = amp*COS(2.5_num*k*zc(iz))   &
                     * (1.0_num+COS(k*zc(iz)))*(1.0_num - (rc/r1)**2)**mu
                vy(ix,iy,iz) = amp*SIN(2.5_num*k*zc(iz))   &
                     * (1.0_num+COS(k*zc(iz)))*(1.0_num - (rc/r1)**2)**mu
             ELSE
                vx(ix,iy,iz) = 0.0_num
                vy(ix,iy,iz) = 0.0_num
             END IF
          END DO
       END DO
    END DO


  END SUBROUTINE equilibrium
