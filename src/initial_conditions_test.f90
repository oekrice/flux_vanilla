MODULE initial_conditions

  USE shared_data
  USE strings

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: HandleCustomBlock,Equilibrium


CONTAINS

  !-----------------------------------------------------------------------------
  !This function contains the equilibrium
  !-----------------------------------------------------------------------------

  SUBROUTINE Equilibrium


INTEGER:: ix, iy

vx = 0.0_num
vy = 0.0_num
vz = 0.0_num
bx = 0.0_num
by = 0.0_num
bz = 0.0_num

gamma = 5.0_num / 3.0_num

! Energy profile

   DO iy = -1, ny+2
      IF(yc(iy) < 0.0_num) THEN
         energy(:,iy) = (1.0_num - yc(iy)*(gamma - 1.0_num)/gamma)
      ELSEIF((yc(iy) >= 0.0_num) .AND. (yc(iy) <= 10.0_num)) THEN
         energy(:,iy) = 1.0
      ELSEIF((yc(iy) > 10.0_num) .AND. (yc(iy) <= 20.0_num)) THEN
         energy(:,iy) = (150.0_num)**((yc(iy) - 10.0_num)/10.0_num)
      ELSE
         energy(:,iy) = 150.0_num
      END IF
   END DO


energy = energy/(gamma - 1.0_num)


! Density profile


   DO iy = -1, ny+2
      IF(yc(iy) < 0.0_num) THEN
         rho(:,iy) = (1.0_num - yc(iy)*(gamma - 1.0_num)/gamma)**(1.0_num/(gamma - 1.0_num))
      ELSEIF((yc(iy) >= 0.0_num) .AND. (yc(iy) <= 10.0_num)) THEN
         rho(:,iy) = exp(-yc(iy))
      ELSEIF((yc(iy) > 10.0_num) .AND. (yc(iy) <= 20.0_num)) THEN
      rho(:,iy) = exp(10.0_num/log(150.0_num)*(1.0_num/((150.0_num)**((yc(iy) - 10.0_num)/10.0_num)) - 1.0_num) - 10.0_num)
         rho(:,iy) = rho(:,iy)*(150.0_num)**((10.0_num - yc(iy))/10.0_num)
      ELSE
   rho(:,iy) = exp(-10.0_num)*exp(10.0_num/log(150.0_num)*(1.0_num/150.0_num - 1.0_num))*exp(-(yc(iy) - 20.0_num)/150.0_num)
         rho(:,iy) = rho(:,iy)/150.0_num
      END IF
   END DO





  END SUBROUTINE Equilibrium

  !-----------------------------------------------------------------------------
  !These functions contain the user input deck elements
  !It doesn't work properly yet, don't touch it
  !-----------------------------------------------------------------------------

  FUNCTION HandleCustomBlock(blockname,Element,Value)
    CHARACTER(len=30),INTENT(IN)::blockname,Element,Value
    INTEGER :: HandleCustomBlock

  !This line must always be present
     HandleCustomBlock=ERR_UNKNOWN_ELEMENT
  
  END FUNCTION HandleCustomBlock


END MODULE initial_conditions
