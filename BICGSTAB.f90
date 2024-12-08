MODULE BICGSTABMOD
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

        IF (SIZE(A, DIM=1) /= SIZE(A, DIM=2)) STOP "ERROR: IMPROPER DIMENSION OF MATRIX A IN BICGSTAB."

        X  = 0.0_WP
        R  = B - MATMUL(A, X)
        RS = R
        RHO   = 1.0_WP; ALPHA = 1.0_WP; OMEGA = 1.0_WP
        V  = 0.0_WP; P  = 0.0_WP
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
            IT = IT + 1
        END DO

        WRITE(*,*) "ITERATIONS REQUIRED :", IT

        RETURN
    END FUNCTION BICGSTAB     
END MODULE BICGSTABMOD

PROGRAM MAIN
    USE BICGSTABMOD
    IMPLICIT NONE

    INTEGER (KIND=SP), PARAMETER             :: M=4, N=4
    REAL    (KIND=WP), DIMENSION(1:M,1:N)    :: A
    REAL    (KIND=WP), DIMENSION(1:M    )    :: X_CALCULATED, X_ACTUAL, B
    INTEGER                                   :: FILE_UNIT

    A(1,:) = [ 1.0_WP, 1.0_WP, 1.0_WP, 1.0_WP]
    A(2,:) = [ 2.0_WP, 3.0_WP, 1.0_WP, 1.0_WP]
    A(3,:) = [ 3.0_WP, 4.0_WP, 1.0_WP, 1.0_WP]
    A(4,:) = [ 3.0_WP, 4.0_WP, 1.0_WP, 2.0_WP]

    X_ACTUAL(:) = [ 5.0_WP, 8.0_WP,10.0_WP,12.0_WP]
    B = MATMUL(A, X_ACTUAL)

    X_CALCULATED = BICGSTAB(A, B)

    FILE_UNIT = 10
    OPEN (UNIT=FILE_UNIT, FILE="BICGSTAB_RESULTS.TXT", STATUS="REPLACE", ACTION="WRITE")
    WRITE(FILE_UNIT,*) "ANALYTICAL =", X_ACTUAL
    WRITE(FILE_UNIT,*) "BICGSTAB   =", X_CALCULATED
    WRITE(FILE_UNIT,*) "ERROR      =", ABS(X_ACTUAL - X_CALCULATED)/X_ACTUAL
    CLOSE(FILE_UNIT  )

    WRITE(*,*) "//"
    WRITE(*,*) "RESULTS HAVE BEEN WRITTEN TO 'BICGSTAB_RESULTS.TXT'."

END PROGRAM MAIN