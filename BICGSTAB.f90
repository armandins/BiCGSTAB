! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ARMAN DINDAR SAFA
!   12/ 15 / 2024
! 
! 
! THIS PROGRAM SOLVES FOR AX = B USING BICONDUGATE STABILITZED METHOD
! 
! INPUT FILE FORMAT IN ORDER: 
! 1. ROWS AND COLS OF MATRIX A
! 2. A MATRIX 
! 3. B MATRIX
! 4. ANALYTICAL SOLUTION IF AVAILABLE ( OPTIONAL )
! --------------
! EXAMPLE INPUT:
! 4 4                    ! 4 x 4 A
! 1.0 1.0 1.0 1.0        !       A
! 2.0 3.0 1.0 1.0
! 3.0 4.0 1.0 1.0
! 3.0 4.0 1.0 2.0
! 50.0 96.0 135.0 143.0  ! B
! 5.0 8.0 10.0 12.0      ! ANALYTICAL
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE BICGSTAB
    IMPLICIT NONE

    INTEGER (KIND=2), PARAMETER    :: SP = KIND(1.000)
    INTEGER (KIND=2), PARAMETER    :: DP = KIND(1.0D0)
    INTEGER (KIND=2), PARAMETER    :: WP = SP

CONTAINS
    FUNCTION BICGSTAB(A, B) RESULT(X)
        IMPLICIT NONE

        REAL (KIND=WP), INTENT(IN )                   :: A (:,:)
        REAL (KIND=WP), INTENT(IN )                   :: B ( : )
        REAL (KIND=WP), DIMENSION(1:SIZE(B, DIM=1))   :: X

        REAL (KIND=WP), DIMENSION(1:SIZE(B, DIM=1))   :: R, RS, V, P, S, T
        REAL (KIND=WP), PARAMETER                     :: E = 1D-33
        REAL (KIND=WP)                                :: RHO, RHO_PREV
        REAL (KIND=WP)                                :: ALPHA, OMEGA, BETA
        REAL (KIND=WP)                                :: NORM_R, NORM_B       
        INTEGER (KIND=SP)                             :: IT=0

        IF (SIZE(A, DIM=1) /= SIZE(A, DIM=2)) STOP "ERROR: IMPROPER DIMENSIONS FOR MATRIX A IN BICGSTAB."

        X      = 0.0_WP
        R      = B - MATMUL(A, X)
        RS     = R
        RHO    = 1.0_WP; ALPHA = 1.0_WP; OMEGA = 1.0_WP
        V      = 0.0_WP; P  = 0.0_WP
        NORM_R = SQRT(DOT_PRODUCT(R, R))
        NORM_B = SQRT(DOT_PRODUCT(B, B))

        DO WHILE (NORM_R .GT. E*NORM_B)
            RHO_PREV = RHO
            RHO      = DOT_PRODUCT(RS, R)
            BETA     = (RHO/RHO_PREV) * (ALPHA/OMEGA)
            P        = R + BETA * (P - OMEGA*V)
            V        = MATMUL(A, P)
            ALPHA    = RHO / DOT_PRODUCT(RS, V)
            S        = R - ALPHA*V
            T        = MATMUL(A, S)
            OMEGA    = DOT_PRODUCT(T, S) / DOT_PRODUCT(T, T)
            X        = X + ALPHA*P + OMEGA*S
            R        = S - OMEGA*T
            NORM_R   = SQRT(DOT_PRODUCT(R, R))
            NORM_B   = SQRT(DOT_PRODUCT(B, B))
            IT       = IT + 1
        END DO

        WRITE(*,*) "REQ. ITERATIONS: ", IT

        RETURN
    END FUNCTION BICGSTAB     
END MODULE BICGSTAB

PROGRAM MAIN
    USE BICGSTABMOD
    IMPLICIT NONE

    INTEGER (KIND=SP), PARAMETER             :: FILE_UNIT_BCG = 10 
    INTEGER (KIND=SP), PARAMETER             :: FILE_UNIT_IN  = 20
    INTEGER (KIND=SP), PARAMETER             :: FILE_UNIT_ERR = 30
    REAL    (KIND=WP), ALLOCATABLE           :: A(:,:), B(:), X_CALCULATED(:), X_ACTUAL(:)
    INTEGER                                  :: M, N, I

    ! Open the input file
    OPEN (UNIT=FILE_UNIT_IN, FILE="BICGSTAB_INPUT.TXT", STATUS="OLD", ACTION="READ")
    READ (FILE_UNIT_IN,*) M, N
    ALLOCATE (A(1:M,1:N), B(1:M), X_CALCULATED(1:M), X_ACTUAL(1:M))
    DO I = 1, M
        READ (FILE_UNIT_IN,*) A(I, :)
    END DO
    READ  (FILE_UNIT_IN,*) B(:)
    READ  (FILE_UNIT_IN,*) X_ACTUAL(:)
    CLOSE (FILE_UNIT_IN)

    X_CALCULATED = BICGSTAB(A, B)

    OPEN  (UNIT=FILE_UNIT_BCG, FILE="BICGSTAB_SOLUTION.TXT", STATUS="REPLACE", ACTION="WRITE")
    WRITE (FILE_UNIT_BCG,*) X_CALCULATED
    CLOSE (FILE_UNIT_BCG)

    OPEN  (UNIT=FILE_UNIT_ERR, FILE="BICGSTAB_ERROR.TXT", STATUS="REPLACE", ACTION="WRITE")
    WRITE (FILE_UNIT_ERR,*) ABS(X_ACTUAL - X_CALCULATED)/X_ACTUAL
    CLOSE (FILE_UNIT_ERR)

    WRITE(*,*) "//"
    WRITE(*,*) "RESULTS HAVE BEEN WRITTEN TO 'BICGSTAB_SOLUTION.TXT' AND 'BICGSTAB_ERROR.TXT'."

    DEALLOCATE (A, B, X_CALCULATED, X_ACTUAL)
END PROGRAM MAIN
