PROGRAM x_acf
    ! Variables
    IMPLICIT NONE
    REAL, PARAMETER   :: Conv = 1e-8
    INTEGER           :: ios, NMat, N, M
    INTEGER           :: i, j, Ind0, Indi
    REAL              :: CorrTime, frameTime, dt
    REAL              :: xdata, Var, Mean
    REAL              :: characTime, Diffusion
    REAL, ALLOCATABLE :: Mat(:), Blocks(:,:), A(:), t(:)
    CHARACTER(len=32) :: arg, DataFile, SaveFile

    ! Trapezoidal subroutine interface
    INTERFACE
        SUBROUTINE trapezoidal(x, y, area)
            REAL, DIMENSION(:) :: x, y
        END SUBROUTINE trapezoidal
    END INTERFACE

    ! Command line arguments
    CALL GETARG(1, arg)       ! 
    READ(arg, *) CorrTime     ! Autocorrelation cutoff (ps)
    CALL GETARG(2, arg)       !
    READ(arg, *) frameTime    ! Time step, dt, between frames (ps)
    CALL GETARG(3, DataFile)  ! Trajectory File
    CALL GETARG(4, SaveFile)  ! Ouptut File

    ! Determine number of lines in trajectory file
    NMat = 0
    OPEN(UNIT=1, FILE=DataFile)
    DO
        READ(1, *, IOSTAT=ios) xdata
        IF (ios /= 0) EXIT
        NMat = NMat + 1
    END DO
    
    ! Determine block sizes
    N    = CEILING(CorrTime/frameTime)
    NMat = NMat - MOD(NMat,N)
    M    = INT(NMat/N)

    ! Allocate Matrix
    ALLOCATE(Mat(NMat))
    ALLOCATE(Blocks(N,M))
    ALLOCATE(A(N))
    ALLOCATE(t(N))

    ! Initialize time array
    dt   = CorrTime/DBLE(N)
    t    = (/((i*dt), i=1,N)/)

    ! Store data into Matrix
    REWIND(1)
    DO i=1,NMat
        READ(1, *) xdata, Mat(i)
    END DO
    CLOSE(1)

    ! Calculate mean and variance
    Mean = SUM(Mat)/NMat
    Mat  = Mat - Mean
    Var  = SUM(Mat**2)/NMat

    ! Split data into blocks of ensemble
    Blocks = RESHAPE(Mat, (/N, M/))
    DEALLOCATE(Mat)
    
    ! Multiply all point with 1st element - Q(0).Q(t)
    DO i=1,M
        Blocks(:,i) = Blocks(:,i)*Blocks(1,i)
    END DO

    ! Ensemble Average - <Q(0).Q(t)>
    DO i=1,N
        A(i) = SUM(Blocks(i,:))/M
    END DO
    DEALLOCATE(Blocks)

    ! Normalization factor
    A = A / Var

    ! Diffusivity
    CALL trapezoidal(t, A, characTime)
    Diffusion = Var/characTime*Conv

    ! Save ACF
    OPEN(UNIT=2, FILE=SaveFile)
    WRITE(2,'(A, F8.5)') '# Frame time (ps)          : ', frameTime
    WRITE(2,'(A, F8.5)') '# Correlation time (ps)    : ', CorrTime
    WRITE(2,'(A, F8.5)') '# Characteristic time (ps) : ', characTime
    WRITE(2,'(A, F8.5)') '# Diffusion (x10^-9 m^2/s) : ', Diffusion*1e9
    DO i=1,N
        WRITE(2, '(F10.5, F10.5)') t(i), A(i)
    END DO
    CLOSE(2)
    DEALLOCATE(A, t)
END

! Numerical integration
SUBROUTINE trapezoidal(x, y, area)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(in) :: x, y
    REAL, INTENT(OUT)              :: area
    REAL                           :: dx
    INTEGER                        :: i, N

    ! Initialize
    N    = SIZE(x)
    area = 0.0

    ! Calculate area of trapezoid
    DO i=1,N-1
        dx   = x(i+1) - x(i)
        area = area + (dx*(y(i+1) + y(i))/2)
    END DO
END SUBROUTINE trapezoidal
