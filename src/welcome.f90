MODULE welcome

  USE shared_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message

CONTAINS

  SUBROUTINE welcome_message

    CHARACTER, DIMENSION(4) :: clrstr = (/ ' ', '[', '2', 'J' /)

    IF (rank /= 0) RETURN

    clrstr(1) = CHAR(27)
    WRITE(*,'(1x, 4a1)') clrstr

    PRINT*,'    @@          @@@@    @@@@@@@@  @@@@@@@@    @@@@@@@@  @@@@@@   '
    PRINT*,'    @@          @@@@    @@@@@@@@  @@@@@@@@    @@@@@@@@  @@@@@@   '
    PRINT*,'    @@        @@@@@@@@  @@    @@  @@                @@  @@  @@   '
    PRINT*,'    @@        @@@@@@@@  @@    @@  @@                @@  @@  @@   '
    PRINT*,'    @@        @@    @@  @@    @@  @@                @@  @@    @@ '
    PRINT*,'    @@        @@    @@  @@    @@  @@                @@  @@    @@ '
    PRINT*,'    @@        @@@@@@@@  @@@@@@@@  @@@@@@@@    @@@@@@@@  @@    @@ '
    PRINT*,'    @@        @@@@@@@@  @@@@@@@@  @@@@@@@@    @@@@@@@@  @@    @@ '
    PRINT*,'    @@        @@    @@  @@@@      @@                @@  @@    @@ '
    PRINT*,'    @@        @@    @@  @@@@      @@                @@  @@    @@ '
    PRINT*,'    @@        @@    @@  @@  @@    @@                @@  @@  @@   '
    PRINT*,'    @@        @@    @@  @@  @@    @@                @@  @@  @@   '
    PRINT*,'    @@@@@@@@  @@    @@  @@    @@  @@@@@@@@    @@@@@@@@  @@@@@@   '
    PRINT*,'    @@@@@@@@  @@    @@  @@    @@  @@@@@@@@    @@@@@@@@  @@@@@@   '
    PRINT*
    PRINT*

    PRINT *,""
    WRITE(*,'("Welcome to Lare3D Version ",I1,".",I1)') Version,Revision
    PRINT *,""
    IF (.NOT. deckfile .AND. rank .EQ. 0) THEN
      PRINT *,"***SERIOUS WARNING*** No input deck found. Lare3D will abort"
      CALL MPI_ABORT(MPI_COMM_WORLD,errcode)
    ENDIF

  END SUBROUTINE welcome_message

END MODULE welcome
