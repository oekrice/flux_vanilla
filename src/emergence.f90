! All the subroutines in this module are for the partially ionised flux emergence
! simulations; see Leake & Arber, 2006

MODULE emergence
  USE shared_data
  USE boundary

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: perpendicular_resistivity, newton_relax, neutral_fraction

CONTAINS

  
  SUBROUTINE perpendicular_resistivity

    !normalising values are L0 = 150km, v0 = 6.5km/s, rho0 = 2.7e-4 kg/m3
    ! t0 = 23s, T0 = 6420K, P0 = 1.2e4 Pa, B0 = 1200G (0.12T)

    REAL(num), PARAMETER :: eta_bar = 6.56548e-11_num !D=6.4074843321e-11_num
    REAL(num) :: f, b, r, xi_v, bxv, byv, bzv, bfieldsq, rho_v, T_v
    INTEGER :: ixp, iyp, izp

!    CALL neutral_fraction

    DO iz = 0, nz 
       DO iy = 0, ny 
          DO ix = 0, nx 
             ixp = ix + 1
             iyp = iy + 1
             izp = iz + 1

             rho_v = rho(ix,iy,iz)*cv(ix,iy,iz) + rho(ixp,iy,iz)*cv(ixp,iy,iz) &
                  + rho(ix,iyp,iz)*cv(ix,iyp,iz) +rho(ixp,iyp,iz)*cv(ixp,iyp,iz) &
                  + rho(ix,iy,izp)*cv(ix,iy,izp) + rho(ixp,iy,izp)*cv(ixp,iy,izp) &
                  + rho(ix,iyp,izp)*cv(ix,iyp,izp) +rho(ixp,iyp,izp)*cv(ixp,iyp,izp)
             rho_v = rho_v / (cv(ix,iy,iz) + cv(ixp,iy,iz) + cv(ix,iyp,iz) + cv(ixp,iyp,iz)  &
                  + cv(ix,iy,izp) + cv(ixp,iy,izp) + cv(ix,iyp,izp) + cv(ixp,iyp,izp)) 

             bxv = (bx(ix,iy,iz) + bx(ix,iyp,iz) + bx(ix,iy,izp) &
                  + bx(ix,iyp,izp)) / 4.0_num
             byv = (by(ix,iy,iz) + by(ixp,iy,iz) + by(ix,iy,izp) &
                  + by(ixp,iy,izp)) / 4.0_num
             bzv = (bz(ix,iy,iz) + bz(ixp,iy,iz) + bz(ix,iyp,iz) &
                  + bz(ixp,iyp,iz)) / 4.0_num
             bfieldsq = bxv**2 + byv**2 + bzv**2

             T_v = ion_mass * (gamma - 1.0_num) * (energy(ix,iy,iz) &
                      - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot) / (2.0_num - xi_n(ix,iy,iz))  
             T_v = T_v + ion_mass * (gamma - 1.0_num) * (energy(ixp,iy,iz) &
                      - (1.0_num - xi_n(ixp,iy,iz)) * ionise_pot) / (2.0_num - xi_n(ixp,iy,iz))  
             T_v = T_v + ion_mass * (gamma - 1.0_num) * (energy(ix,iyp,iz) &
                      - (1.0_num - xi_n(ix,iyp,iz)) * ionise_pot) / (2.0_num - xi_n(ix,iyp,iz))  
             T_v = T_v + ion_mass * (gamma - 1.0_num) * (energy(ixp,iyp,iz) &
                      - (1.0_num - xi_n(ixp,iyp,iz)) * ionise_pot) / (2.0_num - xi_n(ixp,iyp,iz))  
             T_v = T_v + ion_mass * (gamma - 1.0_num) * (energy(ix,iy,izp) &
                      - (1.0_num - xi_n(ix,iy,izp)) * ionise_pot) / (2.0_num - xi_n(ix,iy,izp))  
             T_v = T_v + ion_mass * (gamma - 1.0_num) * (energy(ixp,iy,izp) &
                      - (1.0_num - xi_n(ixp,iy,izp)) * ionise_pot) / (2.0_num - xi_n(ixp,iy,izp))  
             T_v = T_v + ion_mass * (gamma - 1.0_num) * (energy(ix,iyp,izp) &
                      - (1.0_num - xi_n(ix,iyp,izp)) * ionise_pot) / (2.0_num - xi_n(ix,iyp,izp))  
             T_v = T_v + ion_mass * (gamma - 1.0_num) * (energy(ixp,iyp,izp) &
                      - (1.0_num - xi_n(ixp,iyp,izp)) * ionise_pot) / (2.0_num - xi_n(ixp,iyp,izp))  
             T_v = T_v / 8.0_num

             f = f_bar * SQRT(T_v) * EXP(- MIN(T_bar / T_v,15.0_num)) ! T from b has been cancelled
             b = Tr_bar * EXP(0.25_num * MIN(T_bar / T_v,15.0_num) * (Tr_bar * T_v - 1.0_num))
             r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho_v * b / f))
             xi_v = r / (1.0_num + r) 
             
             f = MAX(1.0_num - xi_v, none_zero)
             eta_perp(ix,iy,iz) = eta_bar * xi_v / f * bfieldsq / rho_v**2 / SQRT(T_v)

             eta_perp(ix,iy,iz) = MIN(eta_perp(ix,iy,iz), 10.0_num)

! Simple eta_perp model

 !            eta_perp(ix,iy,iz) = 400.0_num*bfieldsq*exp(-(zc(iz)-10.0_num)*(zc(iz)-10.0_num)/5.0_num);


          END DO
       END DO
    END DO
    

  END SUBROUTINE perpendicular_resistivity
  


  SUBROUTINE neutral_fraction

    REAL(num) :: bof, r, T, rho0, e0, dx, x
    REAL(num), DIMENSION(2) :: Ta, fa, ba, ra, xi_a
    INTEGER :: loop

    !Variable bof is b/f in the original version

    DO iz = -1, nz+2 
       DO iy = -1, ny+2 
          DO ix = -1, nx+2 
             rho0 = rho(ix,iy,iz)
             e0 = energy(ix,iy,iz)
             Ta = ion_mass * (gamma - 1.0_num) * (/  MAX((e0 - ionise_pot) / 2.0_num, none_zero), e0 /)
             IF (Ta(1) > Ta(2)) THEN 
                PRINT*, "Temperature bounds problem", Ta
                STOP
             ENDIF
             dx = Ta(2) - Ta(1)
             T = Ta(1)

             DO loop = 1, 100
                dx = dx / 2.0_num
                x = T  + dx
                bof=Tr_bar/(f_bar * SQRT(x)) * EXP((0.25_num * (Tr_bar * x - 1.0_num) + 1.0_num) * T_bar / x)
                r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar*rho0 * bof))
                xi_a(1) = r / (1.0_num + r)
                fa(1) = x - ion_mass * (gamma - 1.0_num) * (e0 &
                     - (1.0_num - xi_a(1)) * ionise_pot) / (2.0_num - xi_a(1))  
                IF (fa(1) <= 0.0_num) T = x
                IF (ABS(dx) < 1.e-8_num .OR. fa(1) == 0.0_num) EXIT
             END DO
             bof=Tr_bar/(f_bar * SQRT(T)) * EXP((0.25_num * (Tr_bar * T - 1.0_num) + 1.0_num) * T_bar / T)
             r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho0 * bof))
             xi_n(ix,iy,iz) = r / (1.0_num + r)             
          END DO
       END DO
    END DO


  END SUBROUTINE neutral_fraction
  


  SUBROUTINE newton_relax
    
    INTEGER, DIMENSION(1) :: ref_index, z0(1) = 1
    LOGICAL :: first_call = .TRUE., run_loop = .TRUE.

    ! This should only be run above the photosphere so first call sets up
    ! the lowest value of iz to use if at all, the -2 is due to zc starting at -1
    IF (first_call) THEN
       z0 = MINLOC(ABS(zc - 0.0_num)) - 2      
       IF (z0(1) > nz) run_loop = .FALSE. ! This process doesn't have any cells in the corona
       IF (z0(1) < 1) z0(1) = 1 ! only need to run over the internal domain
       first_call = .FALSE.
    ENDIF
       
    ! For every point need to find the reference density value and hence
    ! the tau and temperature
    IF (run_loop) THEN
       DO iz = z0(1), nz
          DO iy = 1, ny
             DO ix = 1, nx
                ! the 2 is subtracted due to rho_ref starting at -1
                ref_index = MINLOC(ABS(rho(ix,iy,iz) - rho_ref)) - 2
                energy(ix,iy,iz) = (energy(ix,iy,iz) + dt / tau_ref(ref_index(1)) * &
                     T_ref(ref_index(1)) / (gamma - 1.0_num)) &
                     / (1.0_num + dt / tau_ref(ref_index(1)))
             END DO
          END DO
       END DO
  !  IF(rank.eq. 0) print*, z0(1)
    ENDIF
        
!    CALL energy_bcs

  END SUBROUTINE newton_relax
  

END MODULE emergence
