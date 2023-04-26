!! filename = CH.f90
!! Vahid Attari
!! Created: 30 Feb. 2016
!! Modified: ....
!! Arroyave Research Group, Department of Materials Science & Engineering, Texas A&M University
!!
!! Acknowledgements:  Based on Cahn-Hilliard 1965 paper
!!
!! Purpose:
!!   - Phase Field Modeling with dynamic coupling to thermodyanmic and kinetic databases
!!     to self consistantly model the Spinodal Composition Phenomenon
!!   
!! General Algorithm function:
!!
!!   1. Retrieve parameter data from file "parameters.dat"
!!   2. Assess thermodynamics of the associated system 
!!   3. Reads initial phase distribution from "phase.dat" file
!!   4. Calculate Phase Evolution with time integration
!!      -  Nucleate Phases
!!      -  Resolve boundary conditions 
!!         -- Periodic boundaries in all directions
!!      -  Solve differential equations via 9-stencil finite difference
!!      -  Update phase information and concentration data
!!
!! Compilation instructions: >> make
!!    - Manual: >>  ifort -o a.out CH.f90
!!
!! Execution: >> ./a.out 
!!                                     
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!!====================================================================================

!   This code simulates the early stage of Spinodal Decomposition...



!   Cahn- Hilliard solver 
!   with periodic boundary conditions.
!   Nonlinear term: f(u) = u - u**3

module INITIALIZATION

    integer, parameter::nx=128, ny=128

    real(8),parameter :: Lx = 2e-6
    real(8),parameter :: Ly = Lx

    real(8),parameter :: dx = Lx/(real(nx)-1)
    real(8),parameter :: dy = dx
    real(8),parameter :: h  = dx
    
    real(8),parameter :: dt = 1e-12 !   dx)**4/32    !0.0001


    real, dimension(nx,ny) :: phi,lap_phi,f,df,d2f
    real, dimension(nx,ny) :: lap_chem_pot,chem_pot
    real,save :: phiold, pi
    real :: phitot
    real :: Betta_C,Betta_Max,Landa_C,Landa_Max,R_Max,betta,R
    integer, save :: l1, m1
   
end module INITIALIZATION

!==============================================================

module mts_pars

	! Energy parameters
    real :: aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk
	
	! Kinetic parameters
    real,parameter  :: mobility = 1.0e-10

	! Alloy parameters
    real :: C0 
	real, parameter :: noise_amplitude = 0.02

    REAL, PARAMETER :: W       = 1e6
    REAL, PARAMETER :: EPSSL   = sqrt(1e-11)
    REAL, PARAMETER :: kappa   = EPSSL**2
    REAL, PARAMETER :: sigma   = sqrt(EPSSL**2*W)/(3.0*sqrt(2.0))
    REAL, PARAMETER :: delta   = sqrt(EPSSL**2/W)   ! interface width
    
end module mts_pars

!==============================================================

MODULE wrt_opts

    CHARACTER ch*9
    CHARACTER(LEN=100) :: VAL
    CHARACTER(len=255) :: cwd,fileplace

    INTEGER, PARAMETER :: wrt_cycle  = 5000           !100000
    INTEGER, PARAMETER :: file_cycle = 10*wrt_cycle   !20000
    INTEGER, PARAMETER :: stop_iter  = 1*file_cycle   !100000
    INTEGER :: NNN3      = 1000
    CHARACTER(100) :: CMSG
	INTEGER :: CSTAT, ESTAT

    integer            :: IPC
    integer            :: itimes
    real               :: Time

CONTAINS

    SUBROUTINE mk_dir

        !!***** MAKE DIRs *****
        call system('rm -r ../parameters')
		call system('rm -r ../microstructure')
        call system('rm -r ../results')
        call system('rm -r ../images')
        !call system('rm -r heat')
        call system('rm *.mod')
        call system('rm *.dat')
        call system('rm *.plt')

        call system('mkdir ../images')
        call system('mkdir ../results')
        call system('mkdir ../microstructure')
        !call system('mkdir parameters')
        !call system('mkdir heat')

    END SUBROUTINE

    SUBROUTINE wrt_disk

        !!***** MAKE DIRs *****
        CALL getcwd(cwd)
        WRITE(fileplace,*) ADJUSTL(TRIM(cwd))

    END SUBROUTINE

END MODULE wrt_opts

!==============================================================

! -----------------------------------------------
MODULE String_Functions  ! by David Frank  dave_frank@hotmail.com
IMPLICIT NONE            ! http://home.earthlink.net/~dave_gemini/strings.f90

! Copy (generic) char array to string or string to char array
! Clen           returns same as LEN      unless last non-blank char = null
! Clen_trim      returns same as LEN_TRIM    "              "
! Ctrim          returns same as TRIM        "              "
! Count_Items    in string that are blank or comma separated
! Reduce_Blanks  in string to 1 blank between items, last char not blank
! Replace_Text   in all occurances in string with replacement string
! Spack          pack string's chars == extract string's chars
! Tally          occurances in string of text arg
! Translate      text arg via indexed code table
! Upper/Lower    case the text arg

INTERFACE Copy    ! generic
   MODULE PROCEDURE copy_a2s, copy_s2a
END INTERFACE Copy

CONTAINS
! ------------------------
PURE FUNCTION Copy_a2s(a)  RESULT (s)    ! copy char array to string
CHARACTER,INTENT(IN) :: a(:)
CHARACTER(SIZE(a)) :: s
INTEGER :: i
DO i = 1,SIZE(a)
   s(i:i) = a(i)
END DO
END FUNCTION Copy_a2s

! ------------------------
PURE FUNCTION Copy_s2a(s)  RESULT (a)   ! copy s(1:Clen(s)) to char array
CHARACTER(*),INTENT(IN) :: s
CHARACTER :: a(LEN(s))
INTEGER :: i
DO i = 1,LEN(s)
   a(i) = s(i:i)
END DO
END FUNCTION Copy_s2a

! ------------------------
PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
INTEGER :: i
Clen = LEN(s)
i = LEN_TRIM(s)
IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
END FUNCTION Clen

! ------------------------
PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
INTEGER :: i                       ! then len of C string is returned, note:
                                   ! Ctrim is only user of this function
i = LEN_TRIM(s) ; Clen_trim = i
IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
END FUNCTION Clen_trim

! ----------------
FUNCTION Ctrim(s1)  RESULT(s2)     ! returns same result as TRIM unless:
CHARACTER(*),INTENT(IN)  :: s1     ! last non-blank char is null in which
CHARACTER(Clen_trim(s1)) :: s2     ! case trailing blanks prior to null
s2 = s1                            ! are output
END FUNCTION Ctrim

! --------------------
INTEGER FUNCTION Count_Items(s1)  ! in string or C string that are blank or comma separated
CHARACTER(*) :: s1
CHARACTER(Clen(s1)) :: s
INTEGER :: i, k

s = s1                            ! remove possible last char null
k = 0  ; IF (s /= ' ') k = 1      ! string has at least 1 item
DO i = 1,LEN_TRIM(s)-1
   IF (s(i:i) /= ' '.AND.s(i:i) /= ',' &
                    .AND.s(i+1:i+1) == ' '.OR.s(i+1:i+1) == ',') k = k+1
END DO
Count_Items = k
END FUNCTION Count_Items

! --------------------
FUNCTION Reduce_Blanks(s)  RESULT (outs)
CHARACTER(*)      :: s
CHARACTER(LEN_TRIM(s)) :: outs
INTEGER           :: i, k, n

n = 0  ; k = LEN_TRIM(s)          ! k=index last non-blank (may be null)
DO i = 1,k-1                      ! dont process last char yet
   n = n+1 ; outs(n:n) = s(i:i)
   IF (s(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
END DO
n = n+1  ; outs(n:n)  = s(k:k)    ! last non-blank char output (may be null)
IF (n < k) outs(n+1:) = ' '       ! pad trailing blanks
END FUNCTION Reduce_Blanks

! ------------------
FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
CHARACTER(*)        :: s,text,rep
CHARACTER(LEN(s)+100) :: outs     ! provide outs with extra 100 char len
INTEGER             :: i, nt, nr

outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
DO
   i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
   outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
END DO
END FUNCTION Replace_Text

! ---------------------------------
FUNCTION Spack (s,ex)  RESULT (outs)
CHARACTER(*) :: s,ex
CHARACTER(LEN(s)) :: outs
CHARACTER :: aex(LEN(ex))   ! array of ex chars to extract
INTEGER   :: i, n

n = 0  ;  aex = Copy(ex)
DO i = 1,LEN(s)
   IF (.NOT.ANY(s(i:i) == aex)) CYCLE   ! dont pack char
   n = n+1 ; outs(n:n) = s(i:i)
END DO
outs(n+1:) = ' '     ! pad with trailing blanks
END FUNCTION Spack

! --------------------
INTEGER FUNCTION Tally (s,text)
CHARACTER(*) :: s, text
INTEGER :: i, nt

Tally = 0 ; nt = LEN_TRIM(text)
DO i = 1,LEN(s)-nt+1
   IF (s(i:i+nt-1) == text(:nt)) Tally = Tally+1
END DO
END FUNCTION Tally

! ---------------------------------
FUNCTION Translate(s1,codes)  RESULT (s2)
CHARACTER(*)       :: s1, codes(2)
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER            :: i, j

DO i = 1,LEN(s1)
   ch = s1(i:i)
   j = INDEX(codes(1),ch) ; IF (j > 0) ch = codes(2)(j:j)
   s2(i:i) = ch
END DO
END FUNCTION Translate

! ---------------------------------
FUNCTION Upper(s1)  RESULT (s2)
CHARACTER(*)       :: s1
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
INTEGER            :: i

DO i = 1,LEN(s1)
   ch = s1(i:i)
   IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
   s2(i:i) = ch
END DO
END FUNCTION Upper

! ---------------------------------
FUNCTION Lower(s1)  RESULT (s2)
CHARACTER(*)       :: s1
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
INTEGER            :: i

DO i = 1,LEN(s1)
   ch = s1(i:i)
   IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
   s2(i:i) = ch
END DO
END FUNCTION Lower

END MODULE String_Functions

!==============================================================

program main
    use INITIALIZATION
    use wrt_opts
    use mts_pars
    implicit none
    integer :: i,j
    real::DT1
			
	! Read energy and chem_pot coefs:		
	open(10,file='../input/poten900.txt')
    read(10,*) aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk
    close(10)

	! Read energy and chem_pot coefs:
	open(10,file='../input/mts_pars900.txt')
    read(10,*) C0
    close(10)

	
    DT1 = 0.2*DX*DX/(mobility*EPSSL**2)
    Write(*,*) 'DT1       =',DT1
			
			
    WRITE(*,*) 'Nx        =',nx
    WRITE(*,*) 'Lx        =',Lx
    WRITE(*,*) 'DX        =',DX
    WRITE(*,*) 'DT        =',DT
    WRITE(*,*)
    WRITE(*,*) 'Alloy Comp=',C0
    WRITE(*,*) 'Noise Ampl=',noise_amplitude
    WRITE(*,*)
    WRITE(*,*) 'EPSILON_SL=',EPSSL
    WRITE(*,*) 'sigma_SL  =',sigma
    WRITE(*,*) 'Delta_SL  =',delta
	WRITE(*,*) 
			
    IPC    = 0
    itimes = 0
	Time   = 0		    
    !! ************* INITIALIZATION
	call mk_dir		
    call init

    !! *********************************************
	OPEN(unit=3,file='../results/itimes.dat')
	WRITE(3,*) IPC,itimes,phitot,Betta_C,Betta_Max
	!! *********************************************
	    
    !! ************* CALCULATING
    do itimes=1,stop_iter
		
        IF ( itimes.EQ. 1 .or. MOD(itimes,wrt_cycle).EQ.0 ) THEN
            IPC = IPC + 1
            !CALL analyze
            CALL printdata(itimes)
            !! *********************************************
			WRITE(3,*) IPC,itimes,phitot,Betta_C,Betta_Max
			!! *********************************************
        ENDIF
        
        call gradandlaplace
        call evolution

        Time = Time + dt
        
    end do

    print *, 'SIMULATION FINISHED... GO HOME!!!'

    stop

end program main

!=================================================================

subroutine init
    use INITIALIZATION
    use mts_pars, only : C0
    implicit none
    integer::i,j
    real::ranum


    l1=nx
    m1=ny
            
    !setting the domain
    do i=1,nx
        do j=1,ny

            call random_number(ranum)
            ranum = (1.0D0*ranum - 0.5D0)
            phi(i,j)= C0 + 0.02*ranum
            phitot = phitot + phi(i,j)
 
        end do
    end do
    
    
end subroutine init

!=======================================================================================================================

subroutine freeenergy(i,j,FE,dfdc,d2fdc2)
    
    use INITIALIZATION, only : nx,ny,phi
    use mts_pars, only : aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk

    implicit none
    integer :: i,j
    real :: x
    real, dimension(nx,ny) :: FE,dfdc,d2fdc2
    
    
    x = phi(i,j)
	
    FE(i,j)  = aa*x**10 + bb*x**9 + cc*x**8 + dd*x**7 + ee*x**6 + &
                       ff*x**5 + gg*x**4 + hh*x**3 + ii*x**2 + jj*x + kk
                
	dfdc(i,j)= 10*aa*x**9 + 9*bb*x**8 + 8*cc*x**7 + 7*dd*x**6 + 6*ee*x**5 + &
                        5*ff*x**4 + 4*gg*x**3 + 3*hh*x**2 + 2*ii*x + jj
	
    d2fdc2(i,j) = 0.0 ! we don't need this.
    
end subroutine

!=======================================================================================================================

subroutine analyze
    use INITIALIZATION
    use mts_pars
    implicit none
    integer :: i,j,nnn 
        
    PI=4.D0*DATAN(1.D0)
    
    !---------------------------
    !! System Total Mass : For Controlling the mass conservation
    phitot = phitot/DFLOAT((nx)*(ny))
    !---------------------------
        
    !---------------------------
    do i=1,nx; do j=1,ny;
        
        call freeenergy(i,j,f,df,d2f)
        ! wavenumber
        Betta_Max  = 0.5*sqrt(-d2f(i,j)/kappa)
        Betta_C    = Betta_Max*sqrt(2.0D0)

        ! wavelength
        Landa_Max  = 2*PI/Betta_Max                 ! The fastest growing wavelength
        Landa_C    = 2*PI/Betta_C 
        
        ! amplification factor
        R_Max      = mobility*((d2f(i,j))**2/8*kappa)

    end do; end do;
    !---------------------------
    
    !---------------------------
    betta = 0
!    WRITE(4,*) 'ZONE'
    do i=0,100,1
        R          = -mobility*(d2f(10,10))*betta**2-2*mobility*kappa*betta**4
        betta      = betta + 0.02 
        !WRITE(4,*) betta_C/betta,R/R_Max
    enddo
        
            
end subroutine analyze   
    
!================================================================================

subroutine gradandlaplace
    use INITIALIZATION
    use mts_pars
    implicit none
    integer :: i,j,ip,im,jp,jm
    integer :: ipp,imm,jpp,jmm
    real :: result

    do i=1,nx
        do j=1,ny


            !finite difference step
            ip = i+1  ; im = i-1  ; jp = j+1  ; jm = j-1  ; 
            ipp = i+2 ; imm = i-2 ; jpp = j+2 ; jmm = j-2 ;
            
            !periodic boundary condition
            if (im==0) im=l1
            if (jm==0) jm=m1
            if (ip==l1+1) ip=1
            if (jp==m1+1) jp=1
            
            if (imm==-1) imm=l1-1
            if (jmm==-1) jmm=m1-1
            if (imm==0)  imm=l1
            if (jmm==0)  jmm=m1
            
            if (ipp==l1+1) ipp=1
            if (jpp==m1+1) jpp=1
            if (ipp==l1+2) ipp=2
            if (jpp==m1+2) jpp=2

            !---------------------------

            !laplacians (9-stencil)
            lap_phi(i,j) = (2.0*(phi(ip,j)+phi(im,j)+phi(i,jp)+ phi(i,jm)) + phi(ip,jp)+phi(im,jm)+phi(im,jp)+ &
            phi(ip,jm) - 12.0*phi(i,j))/(3.0*dx*dx)
            
            !laplacians (13-stencil)
!            lap_phi(i,j) = (20*phi(i,j) - 8*(phi(ip,j)+phi(im,j)+phi(i,jp)+phi(i,jm)) &
!            + 2* (phi(im,jm) + phi(im,jp) + phi(ip,jm) + phi(ip,jp)) &
!            + phi(imm,j) + phi(ipp,j) + phi(i,jmm) + phi(i,jpp))/(dx*dx*dx*dx)
            
            result = lap_phi(i,j)
            IF(ISNAN(result)) THEN
            WRITE(*,*) i,j
            WRITE(*,*) lap_phi(i,j)
            WRITE(*,*) phi(i,j),phi(ip,j),phi(im,j),phi(i,jp),phi(i,jm), &
            phi(im,jm), phi(im,jp), phi(ip,jm) , phi(ip,jp), &
            phi(imm,j) ,phi(ipp,j) , phi(i,jmm) , phi(i,jpp),dx
            WRITE(*,*) (20*phi(i,j) - 8*(phi(ip,j)+phi(im,j)+phi(i,jp)+phi(i,jm)) &
            + 2*(phi(im,jm) + phi(im,jp) + phi(ip,jm) + phi(ip,jp)) &
            + phi(imm,j) + phi(ipp,j) + phi(i,jmm) + phi(i,jpp))
            WRITE(*,*) (1/dx**4)
            stop
            ENDIF 
            
            call freeenergy(i,j,f,df,d2f)
                                    
            chem_pot(i,j) = df(i,j) - 2*kappa*lap_phi(i,j)
                        
        end do
    end do

end subroutine gradandlaplace
!===============================================================
subroutine evolution
    use INITIALIZATION
    use mts_pars
    implicit none
    integer :: i,j,ip,im,jp,jm
    integer :: ipp,imm,jpp,jmm
    phitot = 0.0D0
     
    do i=1,nx
        do j=1,ny
        
            !finite difference step
            ip = i+1  ; im = i-1  ; jp = j+1  ; jm = j-1  ; 
            ipp = i+2 ; imm = i-2 ; jpp = j+2 ; jmm = j-2 ;
            
            !periodic boundary condition
            if (im==0) im=l1
            if (jm==0) jm=m1
            if (ip==l1+1) ip=1
            if (jp==m1+1) jp=1
            
            if (imm==-1) imm=l1-1
            if (jmm==-1) jmm=m1-1
            if (imm==0)  imm=l1
            if (jmm==0)  jmm=m1
            
            if (ipp==l1+1) ipp=1
            if (jpp==m1+1) jpp=1
            if (ipp==l1+2) ipp=2
            if (jpp==m1+2) jpp=2
            
            ! Laplace of chemical potential (9-stencil) 
            lap_chem_pot(i,j) = (2.0*(chem_pot(ip,j)+chem_pot(im,j)+chem_pot(i,jp)+ &
            chem_pot(i,jm)) + chem_pot(ip,jp)+chem_pot(im,jm)+chem_pot(im,jp)+ &
            chem_pot(ip,jm) - 12.0*chem_pot(i,j))/(3.0*dx*dx)
            
            ! Laplace of chemical potential (13-stencil)
!            lap_chem_pot(i,j) = (20*chem_pot(i,j) - 8*(chem_pot(ip,j)+chem_pot(im,j)+chem_pot(i,jp)+chem_pot(i,jm)) &
!            + 2* (chem_pot(im,jm) + chem_pot(im,jp) + chem_pot(ip,jm) + chem_pot(ip,jp)) &
!            + chem_pot(imm,j) + chem_pot(ipp,j) + chem_pot(i,jmm) + chem_pot(i,jpp))/(dx*dx*dx*dx)

            phi(i,j) = phi(i,j) + mobility*lap_chem_pot(i,j)*dt
      
            phitot = phitot + phi(i,j)

        end do
    end do

end subroutine evolution
!==================================================================

subroutine printdata(k)
    use INITIALIZATION
	use wrt_opts
	use mts_pars
	use String_Functions
    
    IMPLICIT NONE
    INTEGER :: k,i,j
    CHARACTER(len=80) :: FileName2
    
    !! *********** WRITTING OUTPUT *************** 
    WRITE(FileName2,FMT='(A22,I6.6,A4)') trim("../microstructure/phi_"),k,".plt"
    !FileName2 = replace(FileName2,' ','', every=.TRUE.)
    FileName2 = trim(Reduce_Blanks(FileName2))
    write(*,*) FileName2

103 FORMAT('variables =',2X,'"IG"',2X,'"JG"',2X,'"PHI"')

    OPEN(unit=52,file=FileName2)
    WRITE(52,103)
    WRITE(52,*) 'ZONE ','I=',nx,'J= ',ny

!    WRITE(52,770) nx, ny
!770 FORMAT('ZONE',2X,'I=',I5,2X,'J=',I5,2X)

    DO j=1,ny
        DO i=1,nx
            WRITE(52,771) i,j, phi(i,j)
        ENDDO
    END DO

    CLOSE(52)
771 Format(I3,1X,i3,1X,1(1X,F8.4))
    !! *********************************************

!            !status = system('python -m plot_thickness.py &')
!            !if ( status .ne. 0 ) stop 'system: error'
!            call execute_command_line ('python -m plot_thickness ../results/pfm_pde_parts'//ch//'.plt', exitstat=i, wait=.false., CMDSTAT=CSTAT, CMDMSG=CMSG)
!			IF (CSTAT > 0) THEN
!			  PRINT *, 'Command execution failed with error ', TRIM(CMSG)
!			ELSE IF (CSTAT < 0) THEN
!			  PRINT *, 'Command execution not supported'
!			ELSE
!			  PRINT *, 'Command completed with status ', ESTAT
!			END IF
!            !print *, "Exit status of external_prog.exe was ", i
    
    
    !! ************* DISPLAY OUTPUT ****************
    WRITE(*,*) '******Calculating...******'
    WRITE(*,*) 'IPC=',IPC,'Time=',Time,'itimes=',itimes
    WRITE(*,*) 'phi_min=',minval(phi(1:nx,1:ny) ),'phi_max=',maxval(phi(1:nx,1:ny) )
    WRITE(*,*) 'PHITOT=',phitot
    WRITE(*,*) 'chem_pot_min=',minval(chem_pot(1:nx,1:ny) ),'chem_pot_max=',maxval(chem_pot(1:nx,1:ny) )
    WRITE(*,*) 'lap_phi_min =',minval(lap_phi(1:nx,1:ny) ) ,'lap_phi_max =',maxval(lap_phi(1:nx,1:ny) )
    WRITE(*,*) 'Max. W.#.=',Betta_Max,'C. W.#.=',Betta_C
    WRITE(*,*) 'Max. W.L.=',Landa_Max,'C. W.L.=',Landa_C
    WRITE(*,*) 'R_Max=',R_Max
    WRITE(*,*) 'kappa=',kappa,'M=',Mobility
    WRITE(*,*)   
    !! *********************************************
    
    RETURN

end subroutine printdata  
