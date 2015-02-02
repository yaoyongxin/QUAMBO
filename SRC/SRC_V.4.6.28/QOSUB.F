#include "symbol.inc"
!**********************************************************************
! QO_READER!
! this subroutine reads the required  parameters from the INCAR file 
! LQO     = 0: no QO
!           1: after scf.
!**********************************************************************
      SUBROUTINE QO_READER(IU5,IU0)
      USE vaspxml; USE QO_MODULE; USE prec
      IMPLICIT NONE

      INTEGER IU5,IU0
      INTEGER IDUM, N, IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LOPEN,LDUM
      CHARACTER*1 :: CHARAC

      LOPEN=.FALSE.
      OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')
      QO%KEY=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOKEY','=','#',';','I', &
     &   QO%KEY,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) THEN
            WRITE(IU0,*)'Error reading item ''QOKEY'' from file INCAR.'
         ENDIF
      ENDIF
      CALL XML_INCAR_V('QOKEY','I',QO%KEY,RDUM,CDUM,LDUM,CHARAC,N)
      IF(QO%KEY==0)RETURN

      QO%LSKIPDAV = .FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOLSKIPDAV','=','#',';','L', &
     &            IDUM,RDUM,CDUM,QO%LSKIPDAV,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOLSKIPDAV'' from file INCAR.'
         QO%LSKIPDAV = .FALSE.
      ENDIF
      CALL XML_INCAR('QOLSKIPDAV','L',IDUM,RDUM,CDUM,QO%LSKIPDAV,CHARAC,N)

      QO%DFERMIUP=.5_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'QODFERMIUP','=','#',';','F', &
     &            IDUM,QO%DFERMIUP,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QODFERMIUP'' from file INCAR.'
         QO%DFERMIUP=0
      ENDIF
      CALL XML_INCAR('QODFERMIUP','F',IDUM,QO%DFERMIUP,CDUM,LDUM,CHARAC,N)

      QO%DFERMIDN=-100_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'QODFERMIDN','=','#',';','F', &
     &            IDUM,QO%DFERMIDN,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QODFERMIDN'' from file INCAR.'
         QO%DFERMIDN=0
      ENDIF
      CALL XML_INCAR('QODFERMIDN','F',IDUM,QO%DFERMIDN,CDUM,LDUM,CHARAC,N)

      QO%RCUT=1.E2_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'QORCUT','=','#',';','F', &
     &            IDUM,QO%RCUT,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QORCUT'' from file INCAR.'
         QO%RCUT=1.e2_q
      ENDIF
      CALL XML_INCAR('QORCUT','F',IDUM,QO%RCUT,CDUM,LDUM,CHARAC,N)

      QO%CCUT=-.5_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOCCUT','=','#',';','F', &
     &            IDUM,QO%CCUT,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOCCUT'' from file INCAR.'
         QO%CCUT=.5_q
      ENDIF
      CALL XML_INCAR('QOCCUT','F',IDUM,QO%CCUT,CDUM,LDUM,CHARAC,N)

      QO%RBV=100._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'QORBV','=','#',';','F', &
     &            IDUM,QO%RBV,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QORBV'' from file INCAR.'
         QO%RBV=20._q
      ENDIF
      CALL XML_INCAR('QORBV','F',IDUM,QO%RBV,CDUM,LDUM,CHARAC,N)

      QO%LORTHLOC = .FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOLORTHLOC','=','#',';','L', &
     &            IDUM,RDUM,CDUM,QO%LORTHLOC,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOLORTHLOC'' from file INCAR.'
         QO%LORTHLOC = .FALSE.
      ENDIF
      CALL XML_INCAR('QOLORTHLOC','L',IDUM,RDUM,CDUM,QO%LORTHLOC,CHARAC,N)

      QO%LORTHALL = .FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOLORTHALL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,QO%LORTHALL,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOLORTHALL'' from file INCAR.'
         QO%LORTHALL = .FALSE.
      ENDIF
      CALL XML_INCAR('QOLORTHALL','L',IDUM,RDUM,CDUM,QO%LORTHALL,CHARAC,N)

      QO%LWRTQO = .FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOLWRTQO','=','#',';','L', &
     &            IDUM,RDUM,CDUM,QO%LWRTQO,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOLWRTQO'' from file INCAR.'
         QO%LWRTQO = .FALSE.
      ENDIF
      CALL XML_INCAR('QOLWRTQO','L',IDUM,RDUM,CDUM,QO%LWRTQO,CHARAC,N)

      QO%LWSCEL = .FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOLWSCEL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,QO%LWSCEL,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOLWSCEL'' from file INCAR.'
         QO%LWSCEL = .FALSE.
      ENDIF
      CALL XML_INCAR('QOLWSCEL','L',IDUM,RDUM,CDUM,QO%LWSCEL,CHARAC,N)

      QO%NBMIN=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'QONBMIN','=','#',';','I', &
     &   QO%NBMIN,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) THEN
            WRITE(IU0,*)'Error reading item ''QONBMIN'' from file INCAR.'
         ENDIF
      ENDIF
      CALL XML_INCAR_V('QONBMIN','I',QO%NBMIN,RDUM,CDUM,LDUM,CHARAC,N)

      QO%NBMAX=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'QONBMAX','=','#',';','I', &
     &   QO%NBMAX,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) THEN
            WRITE(IU0,*)'Error reading item ''QONBMAX'' from file INCAR.'
         ENDIF
      ENDIF
      CALL XML_INCAR_V('QONBMAX','I',QO%NBMAX,RDUM,CDUM,LDUM,CHARAC,N)

      QO%EWINDN=-100._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOEWINDN','=','#',';','F', &
     &            IDUM,QO%EWINDN,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOEWINDN'' from file INCAR.'
         QO%EWINDN=-100._q
      ENDIF
      CALL XML_INCAR('QOEWINDN','F',IDUM,QO%EWINDN,CDUM,LDUM,CHARAC,N)

      QO%EWINUP= 100._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOEWINUP','=','#',';','F', &
     &            IDUM,QO%EWINUP,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOEWINUP'' from file INCAR.'
         QO%EWINUP= 100._q
      ENDIF
      CALL XML_INCAR('QOEWINUP','F',IDUM,QO%EWINUP,CDUM,LDUM,CHARAC,N)

      SELECT CASE (QO%KEY)
      CASE (0)
        WRITE(QOIO,*)'No QO analysis!'
      CASE (-1)
        WRITE(QOIO,*)'Original QUAMBO with XFQIAN method.'
      CASE (1)
        WRITE(QOIO,*)'XFQIAN''s optimised virtual method.'
      CASE (-2)
        WRITE(QOIO,*)'Original QUAMBO with YXYAO method.'
      CASE (2)
        WRITE(QOIO,*)'YXYAO''s optimised virtual method.'
      CASE (3)
        WRITE(QOIO,*)'PURE AO LCAO TB PARAMETERS.'
      CASE (5)
        WRITE(QOIO,*)'SIMPLE PROJECTED WANNIER FUNCTION.'
      CASE DEFAULT
        WRITE(QOIO,*)'WARNING: ILLEGAL QO%KEY!'; STOP
      END SELECT

      END SUBROUTINE QO_READER

!=======================================================================
! INITIALIZE QO
!=======================================================================
      SUBROUTINE INI_QO(T_INFO,P)
      USE pseudo; USE poscar; USE QO_MODULE; USE prec
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
! Local
      INTEGER NT,LMAX,NMAX,IR,IL,LL,LMBASE
      REAL(q) RTMP
      CHARACTER CSEL*1,LABEL*2,UNITS*1

      IF(QO%KEY.EQ.0)RETURN
      QO%ITYP=>T_INFO%ITYP
      QO%NIONS=T_INFO%NIONS
      WRITE(80,'("# PSW,AEW,RW:")')
      QO%NTYP=T_INFO%NTYP
      QO%NITYP=>T_INFO%NITYP
      ALLOCATE(QOYS(T_INFO%NTYP))
      QO%NPSNL=NPSNL
      QO%NQBND=0; QO%LDIM=0; QO%LMDIM=0
      OPEN(QOIU3,FILE='QOCAR',STATUS='OLD')
      DO NT=1,T_INFO%NTYP
        LMAX=1+MAXVAL(P(NT)%LPS(1:P(NT)%LMAX))
        IF(LMAX*2.NE.P(NT)%LMAX)THEN
        WRITE(*,'(A4,I2,1x,A14,I4,1x,I4)')' NT=',NT,'LMAX(QO.vs.P):',LMAX*2,P(NT)%LMAX
        WRITE(*,*)'Warning: NOT two eps for each l in POTCAR!'
        ENDIF
        READ(QOIU3,*)CSEL,LABEL,LL,NMAX,UNITS
        IF(UNITS.NE.'a'.AND.UNITS.NE.'A'.AND.UNITS.NE.'b'.AND.UNITS.NE.'B')THEN
        WRITE(0,*)'Illegal unit:',UNITS; STOP
        ENDIF
        IF(LL.NE.LMAX)THEN
        WRITE(*,'(" WARNING: NT=",I2," LMAX=",I2," LL=",I2," SET LMAX=LL")')NT,LMAX,LL
        LMAX=LL 
        ENDIF
        IF(INDEX(LABEL,P(NT)%ELEMENT,.TRUE.).EQ.0)THEN
          WRITE(*,'(A9,2(2x,A2))')' ELEMENT=',LABEL,P(NT)%ELEMENT
          WRITE(*,'(A6,2I6  )')' LMAX=',LMAX,QOYS(NT)%LMAX
          WRITE(*,*)'QOCAR and POTCAR mis-match!'
          STOP
        ENDIF
        ALLOCATE(QOYS(NT)%LIN(LMAX),QOYS(NT)%LPS(LMAX))
        ALLOCATE(QOYS(NT)%ALPHA(LMAX))
        QOYS(NT)%ALPHA=0; QOYS(NT)%LIN=.TRUE.
        READ(QOIU3,*)CSEL,QOYS(NT)%LIN
        READ(QOIU3,*)CSEL,QOYS(NT)%ALPHA
        QOYS(NT)%R%NMAX=NMAX
! READ ATOMIC ORBITALS, UNIT is Bohr for r and 1/sqrt(Bohr) for r*Phi !
        ALLOCATE(QOYS(NT)%R%R(NMAX),QOYS(NT)%RW(NMAX,LMAX),QOYS(NT)%QW(0:NPSNL,LMAX))
        QOYS(NT)%R%R=0; QOYS(NT)%RW=0; QOYS(NT)%QW=0
        DO IL=1,LMAX
        READ(QOIU3,*)CSEL,LL
        QOYS(NT)%LPS(IL)=LL
        DO IR=1,NMAX
        READ(QOIU3,*)QOYS(NT)%R%R(IR),QOYS(NT)%RW(IR,IL)
        ENDDO
        ENDDO  ! IL
        READ(QOIU3,*)
        QOYS(NT)%LMAX =LMAX
        QOYS(NT)%LMMAX=SUM(2*QOYS(NT)%LPS+1)
        QO%LDIM =MAX(QOYS(NT)%LMAX ,QO%LDIM )
        QO%LMDIM=MAX(QOYS(NT)%LMMAX,QO%LMDIM)
        QOYS(NT)%LMQD=QOYS(NT)%LMMAX
        DO IL=1,LMAX
          IF(QOYS(NT)%LIN(IL)) CYCLE
          QOYS(NT)%LMQD=QOYS(NT)%LMQD-(2*QOYS(NT)%LPS(IL)+1)
          WRITE(QOIO,'( " NT=",I2," L=",I2," skipped for QO!")')NT,QOYS(NT)%LPS(IL)
        ENDDO
        ALLOCATE(QOYS(NT)%LMCHAR(QOYS(NT)%LMQD))
        LMBASE=0
        DO IL=1,LMAX
         IF(.NOT.QOYS(NT)%LIN(IL)) CYCLE
         SELECT CASE (QOYS(NT)%LPS(IL))
          CASE (0)
           QOYS(NT)%LMCHAR(LMBASE+1)='  s'
          CASE (1)
           QOYS(NT)%LMCHAR(LMBASE+1)=' py'
           QOYS(NT)%LMCHAR(LMBASE+2)=' pz'
           QOYS(NT)%LMCHAR(LMBASE+3)=' px'
          CASE (2)
           QOYS(NT)%LMCHAR(LMBASE+1)='dxy'
           QOYS(NT)%LMCHAR(LMBASE+2)='dyz'
           QOYS(NT)%LMCHAR(LMBASE+3)='dz2'
           QOYS(NT)%LMCHAR(LMBASE+4)='dxz'
           QOYS(NT)%LMCHAR(LMBASE+5)='dx2'
          CASE (3)
           QOYS(NT)%LMCHAR(LMBASE+1)='f-3'
           QOYS(NT)%LMCHAR(LMBASE+2)='f-2'
           QOYS(NT)%LMCHAR(LMBASE+3)='f-1'
           QOYS(NT)%LMCHAR(LMBASE+4)='f-0'
           QOYS(NT)%LMCHAR(LMBASE+5)='f+1'
           QOYS(NT)%LMCHAR(LMBASE+6)='f+2'
           QOYS(NT)%LMCHAR(LMBASE+7)='f+3'
          CASE DEFAULT
           STOP 'LMCHAR: LM OUT OF RANGE!'
         END SELECT
         LMBASE=LMBASE+2*QOYS(NT)%LPS(IL)+1
        ENDDO
        QO%NQBND=QO%NQBND +QOYS(NT)%LMQD*T_INFO%NITYP(NT)
        QOYS(NT)%PSMAXN=P(NT)%PSMAXN
        CALL QONORM(NT,P(NT),UNITS)
        CALL QOGREP(P(NT),NT)
      ENDDO ! NT
      CLOSE(QOIU3)
      WRITE(QOIO,'(A10,I4)')' QO%NQBND=',QO%NQBND
      RETURN
   
      END SUBROUTINE INI_QO

!======================================================================= 
! Set up Wigner_Seitz supercell coordinates.
!=======================================================================
      SUBROUTINE SET_WS_VEC()
      USE QO_MODULE; USE prec
      IMPLICIT NONE
! Local
      INTEGER RDIM(2,3),ISITE,IX,IY,IZ,I

      IF(QO%WS%NGRID(1).LE.0.OR.QO%WS%NGRID(1).GT.30)THEN
        WRITE(0,'(" FETAL ERROR: QO%WS%NGRID=",3I)')QO%WS%NGRID
        STOP 'CHECK KPOINTS FILE!'
      ENDIF

      CALL hamiltonian_wigner_seitz(.TRUE.)
      ALLOCATE(QO%WS%IR(3,QO%WS%DIM),QO%WS%NDEG(QO%WS%DIM))
      QO%WS%IR=0; QO%WS%NDEG=0
      CALL hamiltonian_wigner_seitz(.FALSE.)
      WRITE(QOIO,*)'SET_WS_VEC DONE!'
      RETURN

      END SUBROUTINE SET_WS_VEC

!=======================================================================
! Set up supercell coordinates.
!=======================================================================
      SUBROUTINE SET_CELL()
      USE QO_MODULE; USE prec
      IMPLICIT NONE
! Local
      INTEGER RDIM(2,3),ISITE,IX,IY,IZ,I

      RDIM(1,:) =-FLOOR(ABS((QO%WS%NGRID-0.5_q)/2))
      DO I=1,3
        IF(MOD(QO%WS%NGRID(I),2).EQ.0)THEN
        RDIM(2,I)=-RDIM(1,I)+1
        ELSE
        RDIM(2,I)=-RDIM(1,I)
        ENDIF
      ENDDO
      QO%WS%DIM=(RDIM(2,1)-RDIM(1,1)+1)*(RDIM(2,2)-RDIM(1,2)+1)*(RDIM(2,3)-RDIM(1,3)+1)
      ALLOCATE(QO%WS%IR(3,QO%WS%DIM),QO%WS%NDEG(QO%WS%DIM))
      QO%WS%IR=0; QO%WS%NDEG=1

      WRITE(QOIO,'(1x,A5,6I4)')'RDIM:',RDIM
      ISITE=0
      DO IX=RDIM(1,1),RDIM(2,1); DO IY=RDIM(1,2),RDIM(2,2); DO IZ=RDIM(1,3),RDIM(2,3)
      ISITE=ISITE+1
      QO%WS%IR(1,ISITE)=IX; QO%WS%IR(2,ISITE)=IY; QO%WS%IR(3,ISITE)=IZ
      ENDDO; ENDDO; ENDDO
      WRITE(QOIO,*)'SET_CELL DONE!'
      RETURN

      END SUBROUTINE SET_CELL

!======================================================================= 
! Contract Plane Wave representation to QO rep'n.
! Now just for serial version!
!=======================================================================
      SUBROUTINE PW2QO(GRID,P,T_INFO,SYMM,INFO,LATT_CUR,NONL_S,W,WDES,LMDIM,CDIJ,CQIJ,SV,EFERMI)
      USE pseudo; USE poscar; USE wave; USE dfast; USE base; USE lattice; USE nonl
      USE mgrid ; USE hamil; USE QO_MODULE; USE prec
      IMPLICIT NONE
      TYPE (grid_3d)     GRID          ! descriptor for FFT grids
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (symmetry)    SYMM
      TYPE (info_struct) INFO          ! INFO structure of VASP
      TYPE (latt)        LATT_CUR
      TYPE (nonl_struct) NONL_S        ! descriptor for non local part of PP (reciprocal space)
      TYPE (wavespin)    W             ! array for wavefunction
      TYPE (wavedes)     WDES          ! descriptor for wavefunction
      INTEGER            LMDIM
      RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      REAL(q)            EFERMI
! Local
      INTEGER ISP,NK,NB,NBKP,NBP,HOMO,NBIN,IQ,IQP,NQBND,NPL,LWORK,IERR,NBVIR,NBSUB,NSKIP
      INTEGER NI,NIS,NT,LL,MM,LM,LLP,MMP,LMP,LMQ,LMQP,NBVIRT,LWORKM
      INTEGER,PARAMETER::NKSHOW=1
      GDEF,ALLOCATABLE       :: CPROW(:,:)  ! stores the projected coefficients of wave functions * Q
      TYPE (wavedes1)    WDES1
      TYPE (wavefun1),ALLOCATABLE::W1(:)
      REAL(q) NEL,NELS,RES,ETOP,EBOT,RES2,FERWE
      COMPLEX(q) CES,CONE,CZERO
      COMPLEX(q),ALLOCATABLE :: A(:,:),B(:,:),ZWORK(:),CPJBUF(:,:),CWBUF(:,:),WOA(:,:)
      REAL(q),ALLOCATABLE :: WEIGM(:),RWORKM(:)
      COMPLEX(q),ALLOCATABLE :: AM(:,:),ZWORKM(:)
      REAL(q),ALLOCATABLE :: WEIG(:),RWORK(:),WQQW(:,:,:,:),LOCC(:,:,:,:)
      INTEGER             :: W_NBANDS
      REAL(q),ALLOCATABLE :: W_FERWE(:)
      COMPLEX(q),ALLOCATABLE :: W_CELEN(:),W_CPTWFP(:,:),W_CPROJ(:,:)
      COMPLEX(q),ALLOCATABLE :: W_CPTWFP_(:),W_CPROJ_(:)

! test
!      logical lbin(20)
!      lbin=.true.; lbin(1:5)=.false.

      WRITE(*,'(" PW2QO...")')
#ifdef MPI
      WRITE(*,*)'QO IS ONLY FOR SERIAL VERSION CURRENTLY!'; STOP 
#endif
      IF(INFO%LREAL)THEN
      WRITE(*,*)'SET LREAL=F IN INCAR FOR QO!'; RETURN
      ENDIF

      NQBND=QO%NQBND; CONE=1._q; CZERO=0
      QO%RSPIN=WDES%RSPIN;   QO%NBANDS=WDES%NBANDS
      QO%ISPIN=WDES%ISPIN;   QO%NKPTS =WDES%NKPTS
      QO%NRPLWV=WDES%NRPLWV; QO%NPROD =WDES%NPROD
      WRITE(QOIO,'(" QO:",f12.1," MBYTES are necessary.")')16*QO%NRPLWV*NQBND*QO%NKPTS*2/1.E6_q

      QO%NBSUB=QO%NBMAX-QO%NBMIN+1
      IF(QO%NBSUB.LT.NQBND.AND.QO%KEY.EQ.5)STOP" ERROR: NBSUB<NQBND FOR KEY=5!"
      IF(QO%KEY.NE.5)QO%NBSUB=NQBND
      NBSUB=QO%NBSUB

      ALLOCATE(QO%CPTWFP (QO%NRPLWV,NBSUB                  ), &
              &QO%HOMO   (                QO%NKPTS,QO%ISPIN), &
              &QO%CPROJ  (QO%NPROD ,NBSUB                  ), &
              &QO%CELEN  (          NBSUB,QO%NKPTS,QO%ISPIN), &
              &QO%WATOM  (QO%NRPLWV,NQBND                  ), &
              &QO%WAPRJ  (QO%NPROD ,NQBND                  ), &
              &QO%PSICOEF(NBSUB    ,NQBND,QO%NKPTS,QO%ISPIN), &
              &QO%IONIQ  (QO%LMDIM ,WDES%NIONS)             , &
              &QO%IQION  (NQBND)                            , &
              &QO%IQLMCHAR(NQBND)                             )
      QO%WATOM=0; QO%WAPRJ=0; QO%IONIQ=0; QO%IQION=0
      QO%HOMO=0
      QO%CELEN=0; QO%PSICOEF=0
      QO%EFERMI=EFERMI
      ETOP=QO%DFERMIUP+EFERMI; EBOT=QO%DFERMIDN+EFERMI
      QO%EWINDN=QO%EWINDN+EFERMI; QO%EWINUP=QO%EWINUP+EFERMI
      ALLOCATE(A(NQBND,NQBND),B(NQBND,NQBND),WEIG(NQBND))
      A=0; B=0; WEIG=0
      LWORK=64*NQBND
      ALLOCATE(ZWORK(LWORK),RWORK(3*NQBND-2)); ZWORK=0; RWORK=0
      ALLOCATE(CPROW(QO%NPROD,NBSUB)); CPROW=0
      ALLOCATE(WQQW(QO%LMDIM,QO%LMDIM,WDES%NIONS,QO%ISPIN)); WQQW=0
      ALLOCATE(LOCC(QO%LMDIM,QO%LMDIM,WDES%NIONS,QO%ISPIN)); LOCC=0
      ALLOCATE(CPJBUF(QO%NPROD,NQBND),CWBUF(QO%NRPLWV,NQBND))
      ALLOCATE(W1(NQBND))
      DO IQ=1,NQBND; ALLOCATE(W1(IQ)%CR(GRID%MPLWV)); ENDDO

      IF(QO%KEY.EQ.0)RETURN
      IF(QO%NBMAX.EQ.0)QO%NBMAX=WDES%NBANDS
! Quick check total number of electrons
      NELS=0
      DO ISP=1,QO%ISPIN; DO NK =1,QO%NKPTS
      NELS=NELS+SUM(W%FERWE(QO%NBMIN:QO%NBMAX,NK,ISP))*WDES%WTKPT(NK)*QO%RSPIN
      ENDDO; ENDDO
      WRITE(QOIO,'(" EXACT NUMBER OF ELECTRONS in SUB SPACE:",f7.1)')NELS
      NELS=0

      W_NBANDS=QO%NBMAX-QO%NBMIN+1
      ALLOCATE(W_FERWE(W_NBANDS),W_CELEN(W_NBANDS),W_CPTWFP(QO%NRPLWV,W_NBANDS),W_CPROJ(QO%NPROD,W_NBANDS))
      ALLOCATE(W_CPTWFP_(QO%NRPLWV),W_CPROJ_(QO%NPROD))

      IF(QO%LWRTQO)THEN
        OPEN (QOIU6,FILE='LCAO.QO',STATUS='REPLACE')
        WRITE(QOIU6,'(4I5,"   ! ISPIN,NRPLWV,NK,NQ")')QO%ISPIN,QO%NRPLWV,QO%NKPTS,QO%NQBND
      ENDIF

!=======================================================================
      spin:    DO ISP=1,QO%ISPIN
      kpoints: DO NK =1,QO%NKPTS
!=======================================================================
      W_FERWE =W%FERWE (  QO%NBMIN:QO%NBMAX,NK,ISP)
      W_CELEN =W%CELEN (  QO%NBMIN:QO%NBMAX,NK,ISP)
      W_CPTWFP=W%CPTWFP(:,QO%NBMIN:QO%NBMAX,NK,ISP)
      W_CPROJ =W%CPROJ (:,QO%NBMIN:QO%NBMAX,NK,ISP)
      QO%CPTWFP=0; QO%CPROJ=0
      CALL GEN_LCAO(LATT_CUR,WDES,NONL_S,NK)
! Check states to be preserved!
      HOMO=0; NSKIP=0
      DO NB=1,W_NBANDS
        IF(ABS(W_FERWE(NB)).GT.1.E-4_q)HOMO=NB ! TETRAHEDRAN may have negative occ.
        IF(REAL(W_CELEN(NB),q).GT.ETOP.AND.ABS(W_FERWE(NB)).LT.1.E-4_q)EXIT
        IF(REAL(W_CELEN(NB),q).LT.EBOT)NSKIP=NSKIP+1
      ENDDO
      NBKP=MAX(NB-1,HOMO)
      
! Quick test
      IF(ETOP.LT.EBOT)THEN
        NBKP=0; NSKIP=0
      ENDIF

! Reorder
      DO NB=1,NSKIP
        CES=W_CELEN(1); W_CPTWFP_=W_CPTWFP(:,1); W_CPROJ_=W_CPROJ(:,1); FERWE=W_FERWE(1)
        DO NBP=1,W_NBANDS-1
          W_CELEN (NBP)=W_CELEN(NBP+1)
          W_FERWE (NBP)=W_FERWE(NBP+1)
          W_CPTWFP(:,NBP)=W_CPTWFP(:,NBP+1)
          W_CPROJ (:,NBP)=W_CPROJ (:,NBP+1)
        ENDDO
        W_CELEN (W_NBANDS)=CES
        W_FERWE (W_NBANDS)=FERWE
        W_CPTWFP(:,W_NBANDS)=W_CPTWFP_
        W_CPROJ (:,W_NBANDS)=W_CPROJ_
      ENDDO
      NBKP=NBKP-NSKIP; HOMO=HOMO-NSKIP

! test
!      HOMO=6; NBKP=MAX(NB-1,HOMO)

      DO NB=1,W_NBANDS
        IF(REAL(W_CELEN(NB),q).GT.QO%EWINDN.AND.REAL(W_CELEN(NB),q).LT.QO%EWINUP)CYCLE
        W_CPTWFP(:,NB)=0; W_CPROJ(:,NB)=0
      ENDDO

      CALL SETWDES(WDES,WDES1,NK)
      WDES1%NBANDS=1    ! used this only here not quite clean
      NPL=WDES1%NGVECTOR

      IF(QO%KEY.EQ.3)THEN  ! ATOMIC LCAO
      NBKP=0; CPROW=0
      DO IQ=1,NQBND
        CALL OVERL(WDES1,INFO%LOVERL,LMDIM,CQIJ,QO%WAPRJ(1,IQ),CPROW(1,IQ)) ! Q|LCAO>
      ENDDO
      A=0
      DO IQ=1,NQBND; DO IQP=IQ,NQBND
        CALL WPSW(QO%WATOM(:,IQP),QO%WATOM(:,IQ), &
                 &QO%WAPRJ(:,IQP),CPROW   (:,IQ),NPL,QO%NPROD,A(IQP,IQ))
      IF(IQP.GT.IQ)A(IQ,IQP)=CONJG(A(IQP,IQ))
      ENDDO; ENDDO
      QO%CPTWFP=QO%WATOM
      QO%CPROJ =QO%WAPRJ
      CALL ORTHOK2(A,QO%CPTWFP,QO%CPROJ)
      GOTO 3109
      ENDIF ! QO%KEY.EQ.3

      IF(QO%KEY.EQ.5)THEN ! Simple projected WF
      DO NB=1,NBSUB
        QO%CPTWFP(:,NB)=W_CPTWFP(:,NB)
        QO%CPROJ (:,NB)=W_CPROJ (:,NB)
        QO%CELEN(NB,NK,ISP)=REAL(W_CELEN(NB),q)
      ENDDO
      GOTO 3110
      ENDIF ! QO%KEY.EQ.5

! Calc QIJ|LCAO>
      CPROW=0
      DO IQ=1,NQBND
        CALL OVERL(WDES1,INFO%LOVERL,LMDIM,CQIJ,QO%WAPRJ(1,IQ),CPROW(1,IQ)) ! Q|LCAO>
      ENDDO
! Calc WOA=<Bloch|S|LCAO> for preserved Bloch states
      IF(NBKP.GT.0)THEN
        ALLOCATE(WOA(NBKP,NQBND)); WOA=0
      ENDIF
      DO IQ=1,NQBND; DO NB=1,NBKP

! test
!        if(NB.LT.6.AND.IQ.LT.11)CYCLE
!        if(NB.GT.7.AND.IQ.GE.11)CYCLE

        CALL WPSW(W_CPTWFP(:,NB),QO%WATOM(:,IQ), &
                 &W_CPROJ (:,NB),CPROW   (:,IQ),NPL,QO%NPROD,WOA(NB,IQ))
      ENDDO; ENDDO

! MATERIAL SPECIFIC! OCCUPIED BANDS FILTER
      NBIN=0

! test, for BaFe2As2 
!      RES =SUM(ABS(WOA(1,:))**2)
!      RES2=SUM(ABS(WOA(2,:))**2)
!      IF(RES.GT.RES2)THEN
!        WOA(2,:)=0
!      ELSE
!        WOA(1,:)=0
!      ENDIF

!      RES =SUM(ABS(WOA(6,1:10))**2)
!      RES2=SUM(ABS(WOA(7,1:10))**2)
!      IF(RES.LT.RES2)THEN
!        WOA(6,1:10)=0; WOA(7,11:)=0
!        lbin(6)=.false.; lbin(7)=.true.
!      ELSE
!        WOA(7,1:10)=0; WOA(6,11:)=0
!        lbin(7)=.false.; lbin(6)=.true.
!      ENDIF

      DO NB=1,NBKP
        IF(QO%CCUT.GT.0)THEN
!          DO IQ=1,2 ! GRAPHENE
!            RES=RES+ABS(WOA(NBIN+1,4*IQ-1))**2 ! Pz components
!          ENDDO
!          RES=RES+ABS(WOA(NBIN+1,129))**2 ! Nd-Graphene
          RES=SUM(ABS(WOA(NBIN+1,:))**2)
          RES=SQRT(RES)
        ELSE
          RES=1._q
        ENDIF

        IF(RES.GT.QO%CCUT)THEN ! CRITERIA
          NBIN=NBIN+1
          QO%CPTWFP(:,NBIN)=W_CPTWFP(:,NBIN)
          QO%CPROJ (:,NBIN)=W_CPROJ (:,NBIN)
          QO%CELEN (  NBIN,NK,ISP)=REAL(W_CELEN(NBIN),q)
          WRITE(QOIO,'(" NB=",I3," WEIGHT=",f8.3," ACCEPTED!")')NB,RES
          IF(NBIN.EQ.NQBND) EXIT
          GOTO 110
        ENDIF

        IF(NB.LE.HOMO)THEN
          WRITE(0,'(" WARNING: OCCUPIED BAND BUT SMALL WEIGHT!")')
          WRITE(0,'(" NB=",I3," WEIGHT=",f8.3)')NB,RES
        ENDIF
        CES=W_CELEN(1+NBIN); W_CPTWFP_=W_CPTWFP(:,1+NBIN); W_CPROJ_=W_CPROJ(:,1+NBIN); FERWE=W_FERWE(1+NBIN)
        DO NBP=1+NBIN,NBKP-1
          WOA(NBP,:)=WOA(NBP+1,:)
          W_FERWE(NBP)=W_FERWE(NBP+1)
          W_CELEN(NBP)=W_CELEN(NBP+1)
          W_CPTWFP(:,NBP)=W_CPTWFP(:,NBP+1)
          W_CPROJ (:,NBP)=W_CPROJ (:,NBP+1)
        ENDDO
        W_FERWE(NBKP)=FERWE
        W_CELEN(NBKP)=CES
        W_CPTWFP(:,NBKP)=W_CPTWFP_
        W_CPROJ (:,NBKP)=W_CPROJ_
        WRITE(QOIO,'(" NB=",I3," WEIGHT=",f8.3," REJECTED!")')NB,RES
110     CONTINUE
        IF(NB.EQ.HOMO)QO%HOMO(NK,ISP)=NBIN
      ENDDO ! NB=1,NBKP
      NBKP=NBIN
      HOMO=QO%HOMO(NK,ISP)
      IF(MOD(NK,NKSHOW)==0)THEN
      WRITE(QOIO,'(" NK=",I4," ISP=",I2," HOMO=",I3," GROSS NBKP=",I3)')NK,ISP,HOMO,NBKP
      ENDIF

      IF(NBKP.GT.0)THEN
        B(1:NBKP,:)=WOA(1:NBKP,:); DEALLOCATE(WOA)
      ENDIF

      IF(NBKP.GT.NQBND)THEN
        WRITE(0,'(" FETAL ERROR: NK=",I4," ISP=",I2," NBKP=",I3," NQBND=",I3)')
        WRITE(0,'(" NBKP>NQBND! INCREASE the dimension of the basis set orbitals!")')
        STOP
      ENDIF
      IF(NBKP.EQ.NQBND)GOTO 3110

!*********************************************************************
      ov: IF(QO%KEY.EQ.1)THEN ! XFQ's method
! PROJECT OUT: I-sum_{|Psi_occ><Psi_occ|S}|AO>
      DO IQ=1,NQBND
      CWBUF (:,IQ)=QO%WATOM(:,IQ)
      CPJBUF(:,IQ)=QO%WAPRJ(:,IQ)
      DO NB=1,NBKP
        CWBUF (:,IQ)=CWBUF (:,IQ)-W_CPTWFP(:,NB)*B(NB,IQ)
        CPJBUF(:,IQ)=CPJBUF(:,IQ)-W_CPROJ (:,NB)*B(NB,IQ)
      ENDDO ! NB
      ENDDO ! IQ, Still a bit different from XF's version.
! Calc B=<AO(S)|I-sum_{|Psi_occ><Psi_occ|S}|AO>
      DO IQ=1,NQBND; DO NB=1,NQBND
        CALL WPSW(CWBUF (:,NB),QO%WATOM(:,IQ), &
                 &CPJBUF(:,NB),CPROW   (:,IQ),NPL,QO%NPROD,B(NB,IQ))
      ENDDO; ENDDO
! A=sum_LCAO {<Psi|LCAO><LCAO|Psi'>}
      A=MATMUL(B,CONJG(TRANSPOSE(B)))
! Calc B=<Bloch_buf|S|Bloch_buf'>
      CPROW=0
      DO IQ=1,NQBND
        CALL OVERL(WDES1,INFO%LOVERL,LMDIM,CQIJ,CPJBUF(1,IQ),CPROW(1,IQ)) ! Q|Bloch_buf>
      ENDDO
      B=0
      DO IQ=1,NQBND; DO IQP=1,NQBND
        CALL WPSW(CWBUF (:,IQP),CWBUF(:,IQ), &
                 &CPJBUF(:,IQP),CPROW(:,IQ),NPL,QO%NPROD,B(IQP,IQ))
      ENDDO; ENDDO
      CALL ZHEGV(1,'V','U',NQBND,A,NQBND,B,NQBND,WEIG,ZWORK,LWORK,RWORK,IERR)
      IF(IERR.NE.0)THEN
        WRITE(*,*)'IERR=',IERR; STOP 'Error: ZHEGV!'
      ENDIF
      IF(MOD(NK,NKSHOW)==0)THEN
        WRITE(QOIO,'(10f7.4)')WEIG
      ENDIF
! antibonding parts for Bloch basis.
      IF(ABS(WEIG(NBKP+1)).LT..01_q)STOP'Check WEIG1!'
      DO NB=NBKP+1,NQBND; DO IQ=1,NQBND
        QO%CPTWFP(:,NB)=QO%CPTWFP(:,NB)+CWBUF (:,IQ)*A(IQ,NB)
        QO%CPROJ (:,NB)=QO%CPROJ (:,NB)+CPJBUF(:,IQ)*A(IQ,NB)
      ENDDO; ENDDO

! YYX's method
      ELSEIF(ABS(QO%KEY).EQ.2)THEN ov 
! A=sum_psi {<LCAO|Psi><Psi|LCAO'>}
      A=MATMUL(CONJG(TRANSPOSE(B(1:NBKP,:))),B(1:NBKP,:))
! Calc B=<LCAO'|S|LCAO>
      B=0
      DO IQ=1,NQBND; DO IQP=1,NQBND
        CALL WPSW(QO%WATOM(:,IQP),QO%WATOM(:,IQ), &
                 &QO%WAPRJ(:,IQP),CPROW   (:,IQ),NPL,QO%NPROD,B(IQP,IQ)) ! IQP,IQ order is important!
      ENDDO; ENDDO
      CALL ZHEGV(1,'V','U',NQBND,A,NQBND,B,NQBND,WEIG,ZWORK,LWORK,RWORK,IERR)
      IF(IERR.NE.0)THEN
        WRITE(*,*)'IERR=',IERR; STOP 'Error: ZHEGVI!'
      ENDIF
      IF(MOD(NK,NKSHOW)==0)THEN
        WRITE(QOIO,'(10F8.3)')WEIG
      ENDIF
      NBVIR=NQBND-NBKP
      IF(NBVIR+1.LE.NQBND)THEN
      IF(WEIG(NBVIR+1).LT..1_q)THEN
        WRITE(QOIO,'(" WARNING: TOO SMALL WEIGI!")')
        WRITE(0,'(" WARNING: TOO SMALL WEIGI: NBKP=",I5)')NBKP
        WRITE(0,'(10F8.3)')WEIG; READ(*,*)
      ENDIF; ENDIF
! antibonding parts for Bloch basis.
      IF(ABS(WEIG(NBVIR)).GT.1.E-6_q)STOP'Check WEIGI!'

! YYX's QO scheme
      IF(QO%KEY.EQ.2)THEN
      DO NB=NBKP+1,NQBND; DO IQ=1,NQBND
        QO%CPTWFP(:,NB)=QO%CPTWFP(:,NB)+QO%WATOM(:,IQ)*A(IQ,NB-NBKP)
        QO%CPROJ (:,NB)=QO%CPROJ (:,NB)+QO%WAPRJ(:,IQ)*A(IQ,NB-NBKP)
      ENDDO; ENDDO
! YYX's original quambo scheme: need virtual bands as input
      ELSE
      NBVIRT=W_NBANDS-NBKP
      ALLOCATE(AM(MAX(NBVIRT,NQBND),MAX(NBVIRT,NQBND))); AM=0
      DO IQ=1,NQBND; DO NB=NBKP+1,W_NBANDS
      CALL WPSW(W_CPTWFP(:,NB),QO%WATOM(:,IQ), &
               &W_CPROJ (:,NB),CPROW   (:,IQ),NPL,QO%NPROD,AM(NB-NBKP,IQ))
      ENDDO; ENDDO
      AM(1:NBVIRT,1:NQBND)=MATMUL(AM(1:NBVIRT,1:NQBND),A)
      AM(1:NBVIRT,1:NBVIRT)=-MATMUL(AM(1:NBVIRT,1:NBVIR),CONJG(TRANSPOSE(AM(1:NBVIRT,1:NBVIR))))

      LWORKM=64*NBVIRT
      ALLOCATE(WEIGM(NBVIRT),ZWORKM(LWORK),RWORKM(3*NBVIRT-2))
      CALL ZHEEV('V','U',NBVIRT,AM(1:NBVIRT,1:NBVIRT),NBVIRT,WEIGM,ZWORKM,LWORK,RWORKM,IERR)
      IF(IERR.NE.0)THEN
        WRITE(0,'(" IERR=",I3)')IERR
        STOP 'FETAL ERROR IN ZHEEV-IV!'
      ENDIF
      WRITE(QOIO,'(" WEIGII:")')
      WRITE(QOIO,'(10F8.3)')WEIGM
      IF(ABS(WEIGM(NBVIR)).LT..1_q)THEN
        WRITE(QOIO,'(" WARNING: WEIGII FOR QUAMBO PROCEDURE TOO SMALL!")')
        WRITE(0,'(" WARNING--WEIGII FOR QUAMBO PROCEDURE:")')
        WRITE(0,'(10F8.3)')WEIGM; READ(*,*)
      ENDIF
      DO IQ=1,NBVIR; DO NB=1,NBVIRT
        QO%CPTWFP(:,IQ+NBKP)=QO%CPTWFP(:,IQ+NBKP)+W_CPTWFP(:,NB+NBKP)*AM(NB,IQ)
        QO%CPROJ (:,IQ+NBKP)=QO%CPROJ (:,IQ+NBKP)+W_CPROJ (:,NB+NBKP)*AM(NB,IQ)
      ENDDO; ENDDO
      DEALLOCATE(AM,WEIGM,ZWORKM,RWORKM)
      ENDIF
      ENDIF ov

3109  CONTINUE
! I actually need the local eigen function
! Construct local sub-Hamiltonian matrix
      NBVIR=NQBND-NBKP
      A=0
      DO IQ=NBKP+1,NQBND
      W1(IQ)%CPTWFP=>QO%CPTWFP(:,IQ); W1(IQ)%CPROJ =>QO%CPROJ (:,IQ)
      W1(IQ)%CR=0
      CALL FFTWAV(NPL,WDES%NINDPW(1,NK),W1(IQ)%CR(1),W1(IQ)%CPTWFP(1),GRID)
      DO IQP=NBKP+1,IQ
      CALL ECCP(WDES1,W1(IQ),W1(IQP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),A(IQP,IQ))
      IF(IQP.LT.IQ)A(IQ,IQP)=CONJG(A(IQP,IQ))
      ENDDO; ENDDO
      CALL ZHEEV('V','U',NBVIR,A(NBKP+1:NQBND,NBKP+1:NQBND),NBVIR,WEIG(NBKP+1:NQBND),ZWORK,LWORK,RWORK,IERR)
      IF(IERR.NE.0)THEN
        WRITE(0,'(" IERR=",I3)')IERR
        STOP 'FETAL ERROR IN ZHEEVII!'
      ENDIF

      CWBUF(:,1:NBVIR) =QO%CPTWFP(:,NBKP+1:NQBND)
      CALL ZGEMM('N','N',QO%NRPLWV,NBVIR,NBVIR,CONE,CWBUF (:,1:NBVIR),QO%NRPLWV,A(NBKP+1:NQBND,NBKP+1:NQBND), &
          &NBVIR,CZERO,QO%CPTWFP(:,NBKP+1:NQBND),QO%NRPLWV)
      CPJBUF(:,1:NBVIR)=QO%CPROJ (:,NBKP+1:NQBND)
      CALL ZGEMM('N','N',QO%NPROD ,NBVIR,NBVIR,CONE,CPJBUF(:,1:NBVIR),QO%NPROD ,A(NBKP+1:NQBND,NBKP+1:NQBND), &
          &NBVIR,CZERO,QO%CPROJ (:,NBKP+1:NQBND),QO%NPROD )
      QO%CELEN(NBKP+1:NQBND,NK,ISP)=WEIG(NBKP+1:NQBND)
      IF(QO%KEY.EQ.1.OR.QO%KEY.EQ.2)THEN
      DO IQ=NBKP+1,NQBND
        IF(QO%CELEN(IQ,NK,ISP).LT.QO%CELEN(HOMO,NK,ISP))THEN
        WRITE(QOIO,*)
        WRITE(QOIO,'(" WARNING: POSSIBLE GHOST BAND ENERGY (HOMO/IQ):",2F8.3)') &
              &QO%CELEN(HOMO,NK,ISP),QO%CELEN(IQ,NK,ISP)
        ENDIF
      ENDDO
      ENDIF ! QO%KEY.EQ.1.OR.QO%KEY.EQ.2

3110  CONTINUE
! Calc Q|PSI>
      CPROW=0
      DO IQ=1,NBSUB
        CALL OVERL(WDES1,INFO%LOVERL,LMDIM,CQIJ,QO%CPROJ(1,IQ),CPROW(1,IQ))
      ENDDO
! Check normalization of QO%CPTWFP
      IF(MOD(NK,NKSHOW)==0)THEN
        WRITE(QOIO,*)'Check <PSI|PSIP>'
        DO NB=1,NBSUB
        CALL WPSW(QO%CPTWFP(:,1),QO%CPTWFP(:,NB), &
                 &QO%CPROJ (:,1),CPROW    (:,NB),NPL,QO%NPROD,CES)
        WRITE(QOIO,'(2(f9.5,1x))')CES
        ENDDO
        WRITE(QOIO,*)'Check <PSI|PSI>' 
        DO NB=1,NBSUB
        CALL WPSW(QO%CPTWFP(:,NB),QO%CPTWFP(:,NB), & 
                 &QO%CPROJ (:,NB),CPROW    (:,NB),NPL,QO%NPROD,CES) 
        WRITE(QOIO,'(2(f9.5,1x))')CES
        ENDDO 
      ENDIF
! Record the coefficients of QO
      DO IQ=1,NQBND; DO NB=1,NBSUB

! test
!        if(.NOT.lbin(NB).AND.IQ.LT.11)CYCLE
!        if(lbin(NB).and.IQ.GE.11)cycle

        CALL WPSW(QO%CPTWFP(:,NB),QO%WATOM(:,IQ), &
                 &   CPROW (:,NB),QO%WAPRJ(:,IQ),NPL,QO%NPROD,QO%PSICOEF(NB,IQ,NK,ISP))
      ENDDO; ENDDO
! Orthogonal QUAMBO
      IF(QO%LORTHALL)THEN
        CALL ORTHOK_ALL(NK,ISP)
      ELSEIF(QO%LORTHLOC)THEN
        CALL ORTHOK_ONSITE(NK,ISP)
      ENDIF

      IF(QO%LWRTQO)THEN
        QO%CPTWFP(:,1:NQBND)=MATMUL(QO%CPTWFP,QO%PSICOEF(:,:,NK,ISP))
        DO IQ=1,NQBND; WRITE(QOIU6,*)QO%CPTWFP(:,IQ); ENDDO
      ENDIF

      NELS=NELS+SUM(W_FERWE(1:QO%HOMO(NK,ISP)))*WDES%WTKPT(NK)*QO%RSPIN
! Normalization for QUAMBOs, symmetrization is needed.
      IQ=0; NIS=0
      DO NT =1,NONL_S%NTYP;   DO NI=1,NONL_S%NITYP(NT); NIS=NIS+1
      LMQ=0
      DO LL =1,QOYS(NT)%LMAX; IF(.NOT.QOYS(NT)%LIN(LL ))CYCLE
      DO MM =1,2*QOYS(NT)%LPS(LL)+1
      LM=MM
      IF(LL.GT.1)LM=LM+SUM(2*QOYS(NT)%LPS(1:LL-1)+1)
      LMQ=LMQ+1
      LMQP=0
      DO LLP=1,QOYS(NT)%LMAX; IF(.NOT.QOYS(NT)%LIN(LLP))CYCLE
      DO MMP=1,2*QOYS(NT)%LPS(LLP)+1
      LMP=MMP
      IF(LLP.GT.1)LMP=LMP+SUM(2*QOYS(NT)%LPS(1:LLP-1)+1)
      LMQP=LMQP+1
      WQQW(LMP,LM,NIS,ISP)=WQQW(LMP,LM,NIS,ISP)+SUM(CONJG(QO%PSICOEF(:,IQ+LMQP,NK,ISP))* &
                      &QO%PSICOEF(:,IQ+LMQ,NK,ISP))*WDES%WTKPT(NK)
      LOCC(LMP,LM,NIS,ISP)=LOCC(LMP,LM,NIS,ISP)+SUM(CONJG(QO%PSICOEF(:,IQ+LMQP,NK,ISP))* &
               &W_FERWE(:)*QO%PSICOEF(:,IQ+LMQ,NK,ISP))* &
               &WDES%WTKPT(NK)*QO%RSPIN
      ENDDO; ENDDO ! LLP,MMP
      ENDDO; ENDDO ! LL ,MM
      IQ=IQ+QOYS(NT)%LMQD
      ENDDO; ENDDO ! NT,NI
      ENDDO kpoints
      ENDDO spin 

      IF(QO%LWRTQO) CLOSE(QOIU6)
      DEALLOCATE(W_FERWE,W_CELEN,W_CPTWFP,W_CPROJ)
      QO%NET=NELS
      WRITE(QOIO,'(" EXACT NUMBER OF ELECTRONS in FILTERED SUB SPACE:",f7.1)')NELS

      NELS=0
      DO ISP=1,QO%ISPIN
      IQ=0; NIS=0
      DO NT=1,NONL_S%NTYP
      DO NI=1,NONL_S%NITYP(NT); NIS=NIS+1
      LMQ=0
      DO LL =1,QOYS(NT)%LMAX; IF(.NOT.QOYS(NT)%LIN(LL ))CYCLE
      DO MM =1,2*QOYS(NT)%LPS(LL)+1
      LM=MM
      IF(LL.GT.1)LM=LM+SUM(2*QOYS(NT)%LPS(1:LL-1)+1)
      LMQ=LMQ+1
      QO%PSICOEF(:,IQ+LMQ,:,ISP)=QO%PSICOEF(:,IQ+LMQ,:,ISP)/SQRT(WQQW(LM,LM,NIS,ISP))
      DO LLP=1,QOYS(NT)%LMAX; IF(.NOT.QOYS(NT)%LIN(LLP))CYCLE
      DO MMP=1,2*QOYS(NT)%LPS(LLP)+1
      LMP=MMP
      IF(LLP.GT.1)LMP=LMP+SUM(2*QOYS(NT)%LPS(1:LLP-1)+1)
      LOCC(LM,LMP,NIS,ISP)=LOCC(LM,LMP,NIS,ISP)/SQRT(WQQW(LM,LM,NIS,ISP))/SQRT(WQQW(LMP,LMP,NIS,ISP))
      ENDDO; ENDDO; ENDDO; ENDDO
      IQ=IQ+QOYS(NT)%LMQD
      ENDDO; ENDDO ! NT,NI

      WRITE(QOIO,*)
      WRITE(QOIO,*) 'DEVIATION:'
      DO NIS=1,NONL_S%NIONS
        WRITE(QOIO,'(I4,1x,16f7.3)') NIS,(             WQQW(LM,LM,NIS,ISP)  ,LM=1,QO%LMDIM)
        WRITE(QOIO,'(I4,1x,16f7.3)') NIS,(SQRT(2*ABS(1-WQQW(LM,LM,NIS,ISP))),LM=1,QO%LMDIM)
      ENDDO
      WRITE(QOIO,*)
      WRITE(QOIO,*) 'LOCC: (MEANINGFUL FOR LORTHALL=.T.)'
      DO NIS=1,NONL_S%NIONS 
        NEL=0
        DO MM=1,QO%LMDIM; NEL=NEL+LOCC(MM,MM,NIS,ISP); ENDDO
        NELS=NELS+NEL
        WRITE(QOIO,'(I4,1x,f12.3)') NIS,NEL
        DO LM=1,QO%LMDIM
        WRITE(QOIO,'(16f7.3)')LOCC(LM,1:QO%LMDIM,NIS,ISP)
        ENDDO
      ENDDO; ENDDO ! ISP
      WRITE(QOIO,'(A14,f12.3)')' TOTAL CHARGE:',NELS
      DEALLOCATE(WQQW,LOCC)

      IF(SYMM%ISYM.NE.0)THEN
        WRITE(0,*)'WARNING: SET ISYM=0 INCAR TO GET TB PARAMETERS!'
        WRITE(0,*)'YOU MIGHT TRY A NON-SELFCONSISTENT LOOP!'
      ENDIF
      CALL READ_POSREF(LATT_CUR,NONL_S)
      IF(QO%LWSCEL)THEN
        CALL SET_WS_VEC()
      ELSE
        CALL SET_CELL()
      ENDIF
      ALLOCATE(QO%TBH(NQBND,NQBND,QO%WS%DIM,QO%ISPIN), &
              &QO%TBS(NQBND,NQBND,QO%WS%DIM,QO%ISPIN))
      QO%TBH=0; QO%TBS=0
      CALL TBPAR(WDES,NONL_S,SYMM,LATT_CUR)
      CALL OUT_TBHS(0)
      CALL OUT_TBHS(1)
      DO IQ=1,NQBND; IF(ASSOCIATED(W1(IQ)%CR))DEALLOCATE(W1(IQ)%CR); ENDDO
      DEALLOCATE(W1)
      CALL TBTEST(WDES,W)
      CALL TBBAND(LATT_CUR)
      RETURN

      END SUBROUTINE PW2QO

!=======================================================================
! Assuming standard Cartesian format
!=======================================================================
      SUBROUTINE READ_POSREF(LATT_CUR,NONL_S)
      USE nonl; USE lattice; USE QO_MODULE; USE prec
      IMPLICIT NONE
      TYPE (latt)        LATT_CUR
      TYPE (nonl_struct) NONL_S        ! descriptor for non local part of PP (reciprocal space)
! LOCAL
      INTEGER I,J,NSCALE
      REAL(q)        SCALEX,SCALEY,SCALEZ
      CHARACTER*255  INPLIN,INPWRK
      CHARACTER*1 STMP
      INTEGER, EXTERNAL :: NITEMS

      ALLOCATE(QO%POSREF(3,QO%NIONS))
      OPEN(QOIU6,FILE='POSREF',STATUS='OLD',ERR=100)
      READ(QOIU6,*)
! one scaling parameter or one for x, y and z
      READ(QOIU6,'(A)') INPLIN
! how many words/data items? --> number of ion types on file POSCAR!
      NSCALE=NITEMS(INPLIN,INPWRK,.TRUE.,'F')
      IF (NSCALE==1) THEN
         READ(INPLIN,*) QO%LAT%SCALE
         SCALEX=1
         SCALEY=1
         SCALEZ=1
      ELSE IF (NSCALE==3) THEN
         QO%LAT%SCALE=1
         READ(INPLIN,*) SCALEX,SCALEY,SCALEZ
      ELSE
         STOP 'ERROR: there must be 1 or 3 items on line 2 of POSREF!'
      ENDIF
      DO I=1,3; READ(QOIU6,*)QO%LAT%A(1:3,I); ENDDO
      IF(QO%LAT%SCALE<0._q)STOP'ERROR: ASSUME SCALE>0 IN POSREF!'
      QO%LAT%A(1,:) =QO%LAT%A(1,:)*SCALEX*QO%LAT%SCALE
      QO%LAT%A(2,:) =QO%LAT%A(2,:)*SCALEY*QO%LAT%SCALE
      QO%LAT%A(3,:) =QO%LAT%A(3,:)*SCALEZ*QO%LAT%SCALE
      CALL LATTIC(QO%LAT)
      READ(QOIU6,*)
      READ(QOIU6,*)STMP
      IF(STMP.EQ.'C'.OR.STMP.EQ.'c'.OR.STMP.EQ.'K'.OR.STMP.EQ.'k')STMP='K'
      DO I=1,QO%NIONS
      READ(QOIU6,*)QO%POSREF(1:3,I)
      ENDDO
      CLOSE(QOIU6)
      IF (STMP=='K') THEN
        QO%POSREF(1,:)=QO%LAT%SCALE*QO%POSREF(1,:)*SCALEX
        QO%POSREF(2,:)=QO%LAT%SCALE*QO%POSREF(2,:)*SCALEY
        QO%POSREF(3,:)=QO%LAT%SCALE*QO%POSREF(3,:)*SCALEZ
        CALL KARDIR(QO%NIONS,QO%POSREF,QO%LAT%B)
      ENDIF
      CALL TOPRIM(QO%NIONS,QO%POSREF)
      DO I=1,3; DO J=1,QO%NIONS
        QO%POSREF(I,J)=QO%POSREF(I,J)+NINT(NONL_S%POSION(I,J)-QO%POSREF(I,J))
      ENDDO; ENDDO
! COMPARE WITH POSCAR
      WRITE(QOIO,'(" POSCAR vs POSREF:")')
      DO I=1,QO%NIONS
        WRITE(QOIO,'(" NI:",I3," POS:",3F7.3)')I,NONL_S%POSION(:,I)
        WRITE(QOIO,'(12X,3F7.3)')QO%POSREF(:,I)
      ENDDO
      DO I=1,3; DO J=1,QO%NIONS
      IF(ABS(QO%POSREF(I,J)-NONL_S%POSION(I,J)).LT..1_q)CYCLE
      WRITE(0,'(" CHECK POS IJ:",2I4)')I,J; STOP
      ENDDO; ENDDO
      RETURN
100   CONTINUE
      WRITE(*   ,'(" NO POSREF, USE POSCAR INSTEAD.")')
      WRITE(QOIO,'(" NO POSREF, USE POSCAR INSTEAD.")')
      QO%LAT%A=LATT_CUR%A; QO%LAT%B=LATT_CUR%B
      QO%LAT%SCALE=LATT_CUR%SCALE
      QO%POSREF=NONL_S%POSION
      RETURN

      END SUBROUTINE READ_POSREF

!================================================================================!
      subroutine hamiltonian_wigner_seitz(count_pts)
    !================================================================================!
    ! Calculates a grid of points that fall inside of (and eventually on the         !
    ! surface of) the Wigner-Seitz supercell centered on the origin of the B         !
    ! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3  !
    !================================================================================!
      USE prec; USE QO_MODULE
    ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
    !                   irvec(1:3,irpt) in the basis of the lattice vectors
    ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
    ! nrpts             number of Wigner-Seitz grid points
    implicit none
    logical, intent(in) :: count_pts
    integer       :: ndiff (3),mp_grid(3),nrpts
    real(q) :: dist(125),tot,dist_min,DPOS(3)
    integer       :: n1,n2,n3,i1,i2,i3,icnt,i,j
    INTEGER,POINTER :: irvec(:,:),ndegen(:)
    ! The Wannier functions live in a supercell of the real space unit cell
    ! this supercell is mp_grid unit cells long in each direction
    ! We loop over grid points r on a unit cell that is 8 times larger than this
    ! primitive supercell.
    ! One of these points is in the W-S cell if it is closer to R=0 than any of the
    ! other points, R (where R are the translation vectors of the supercell)
    ! In the end nrpts contains the total number of grid
    ! points that have been found in the Wigner-Seitz cell
    mp_grid=QO%WS%NGRID
    IF(.NOT.count_pts)THEN
      irvec=>QO%WS%IR; ndegen=>QO%WS%NDEG
    ENDIF
    nrpts = 0
    do n1 = -mp_grid(1) , mp_grid(1)
       do n2 = -mp_grid(2), mp_grid(2)
          do n3 = -mp_grid(3),  mp_grid(3)
             ! Loop over the 125 points R. R=0 corresponds to i1=i2=i3=0, or icnt=63
             icnt = 0
             do i1 = -2, 2
                do i2 = -2, 2
                   do i3 = -2, 2
                      icnt = icnt + 1
                      ! Calculate distance squared |r-R|^2
                      ndiff(1) = n1 - i1 * mp_grid(1)
                      ndiff(2) = n2 - i2 * mp_grid(2)
                      ndiff(3) = n3 - i3 * mp_grid(3)
                      DPOS=MATMUL(QO%LAT%A,ndiff)
                      dist(icnt) = SUM(DPOS**2)
                   enddo
                enddo
             enddo
             dist_min=minval(dist)
             if (abs(dist(63) - dist_min ) .lt. 1.E-7_q ) then
                nrpts = nrpts + 1
                if(.not. count_pts) then
                   ndegen(nrpts)=0
                   do i=1,125
                      if (abs (dist (i) - dist_min) .lt. 1.E-7_q ) ndegen(nrpts)=ndegen(nrpts)+1
                   end do
                   irvec(1, nrpts) = n1
                   irvec(2, nrpts) = n2
                   irvec(3, nrpts) = n3
                endif
             end if
          enddo
       enddo
    enddo
    if(count_pts)THEN
      QO%WS%DIM=nrpts; return
    ENDIF
    write(QOIO,'(1x,i4,a,/)') nrpts,  ' lattice points in Wigner-Seitz supercell:'
    do i=1,nrpts
       write(QOIO,'(4x,a,3(i3,1x),a,i2)') '  vector ', irvec(1,i),irvec(2,i),&
            irvec(3,i),'  degeneracy: ', ndegen(i)
    enddo
    ! Check the "sum rule"
    tot = 0
    do i = 1, nrpts
       tot = tot + 1.0_q/real(ndegen(i),q)
    enddo
    if (abs (tot - real(mp_grid(1) * mp_grid(2) * mp_grid(3),q) ) > 1.E-8_q) then
      STOP 'ERROR in hamiltonian_wigner_seitz: error in finding Wigner-Seitz points!'
    endif
    return

  end subroutine hamiltonian_wigner_seitz

!=======================================================================
      SUBROUTINE CALC_PAIR_DIST(NG,IAG,ICG,NGT)
      USE prec; USE QO_MODULE
      IMPLICIT NONE
      INTEGER NG(2),IAG(MAXVAL(NG),2),NGT,ICG(3,NGT)
! LOCAL
      INTEGER I1,I2,NI1,NI2,ICC(3)
      REAL(q) DPOS(3)
      
      DO I1=1,NG(1); NI1=IAG(I1,1); DO I2=1,NG(2); NI2=IAG(I2,2)
      ICC=ICG(:,I2+NG(1))-ICG(:,I1)
      DPOS=ICC+QO%POSREF(:,NI2)-QO%POSREF(:,NI1)
      DPOS=MATMUL(QO%LAT%A,DPOS)
      WRITE(QOIO,'(" GROUP I ATOM",I3," GROUP II ATOM",I3," DIST=",F8.2)')I1,I2,SQRT(SUM(DPOS**2))
      ENDDO; ENDDO
      RETURN
      END SUBROUTINE CALC_PAIR_DIST

!=======================================================================
      SUBROUTINE ORTHOK_ONSITE(NK,ISP)
      USE QO_MODULE; USE prec
      IMPLICIT NONE
      INTEGER NK,ISP
! LOCAL
      INTEGER NI,NT,LMMAX,LM1,LM2,IQ1,IQ2,IQ3,NB
      INTEGER,PARAMETER::NMAX=35
      COMPLEX(q) SK(NMAX,NMAX),SNHK(NMAX,NMAX)
      INTEGER LWORK,IERR
      REAL(q) RWORK(3*NMAX),W(NMAX)
      COMPLEX(q) WORK(32*NMAX)

      DO NI=1,QO%NIONS; NT=QO%ITYP(NI); LMMAX=QOYS(NT)%LMQD
      IF(LMMAX.LE.1)CYCLE
      IF(LMMAX.GT.NMAX)STOP' NEED TO INCREASE NMAX!'
      SK=0
      DO LM1=1,LMMAX; IQ1=QO%IONIQ(LM1,NI)
      DO LM2=1,LMMAX; IQ2=QO%IONIQ(LM2,NI)
        SK(LM1,LM2)=SUM(CONJG(QO%PSICOEF(:,IQ1,NK,ISP))*QO%PSICOEF(:,IQ2,NK,ISP))
      ENDDO; ENDDO
      LWORK=32*LMMAX
      CALL ZHEEV('V','U',LMMAX,SK(1:LMMAX,1:LMMAX),LMMAX,W(1:LMMAX),WORK,LWORK,RWORK,IERR)
      IF(IERR.NE.0)STOP 'Check ORTHOK_ONSITE!'
      DO IQ3=1,LMMAX
      IF(W(IQ3).LT.1.E-6_q)THEN
      WRITE(0,*)
      WRITE(0,*)'EIGENVALUES:'
      WRITE(0,'(10f8.3)')W(1:LMMAX)
      STOP 'Check ORTHOK_ONSITE!'
      ENDIF
      ENDDO ! IQ3
! Set up S^(-1/2)
      SNHK=0
      DO IQ1=1,LMMAX; DO IQ2=1,LMMAX; DO IQ3=1,LMMAX
        SNHK(IQ2,IQ1)=SNHK(IQ2,IQ1)+SK(IQ2,IQ3)*CONJG(SK(IQ1,IQ3))/SQRT(W(IQ3))
      ENDDO; ENDDO; ENDDO
      DO NB=1,QO%NBSUB
      WORK(1:LMMAX)=QO%PSICOEF(NB,QO%IONIQ(1,NI):QO%IONIQ(LMMAX,NI),NK,ISP)
      QO%PSICOEF(NB,QO%IONIQ(1,NI):QO%IONIQ(LMMAX,NI),NK,ISP)=0
      DO LM1=1,LMMAX; IQ1=QO%IONIQ(LM1,NI)
      DO LM2=1,LMMAX
      QO%PSICOEF(NB,IQ1,NK,ISP)=QO%PSICOEF(NB,IQ1,NK,ISP)+SNHK(LM2,LM1)*WORK(LM2)
      ENDDO; ENDDO; ENDDO
      ENDDO
      RETURN
        
      END SUBROUTINE ORTHOK_ONSITE

!=======================================================================
      SUBROUTINE ORTHOK_ALL(NK,ISP)
      USE QO_MODULE; USE prec
      IMPLICIT NONE
      INTEGER NK,ISP
! LOCAL
      INTEGER NQBND,IQ1,IQ2,IQ3,NB
      COMPLEX(q),ALLOCATABLE :: SK(:,:),SNHK(:,:),WORK(:)
      INTEGER LWORK,IERR
      REAL(q),ALLOCATABLE :: RWORK(:),W(:)

      NQBND=QO%NQBND
      ALLOCATE(SK(NQBND,NQBND),SNHK(NQBND,NQBND),WORK(32*NQBND),RWORK(3*NQBND),W(NQBND))
      SK=0; SNHK=0; WORK=0; RWORK=0; W=0
      DO IQ1=1,NQBND; DO IQ2=1,NQBND
        SK(IQ1,IQ2)=SUM(CONJG(QO%PSICOEF(:,IQ1,NK,ISP))*QO%PSICOEF(:,IQ2,NK,ISP))
      ENDDO; ENDDO
      LWORK=32*NQBND
      CALL ZHEEV('V','U',NQBND,SK,NQBND,W,WORK,LWORK,RWORK,IERR)
      IF(IERR.NE.0)STOP 'Check ORTHOK_ALL!'

      DO IQ3=1,NQBND
      IF(W(IQ3).LT.1.E-6_q)THEN
      WRITE(0,*)
      WRITE(0,*)'EIGENVALUES:'
      WRITE(0,'(10f8.3)')W
      STOP 'Check ORTHOK_ALL!'
      ENDIF
      ENDDO ! IQ3

! Set up S^(-1/2)
      SNHK=0
      DO IQ1=1,NQBND; DO IQ2=1,NQBND; DO IQ3=1,NQBND
        SNHK(IQ2,IQ1)=SNHK(IQ2,IQ1)+SK(IQ2,IQ3)*CONJG(SK(IQ1,IQ3))/SQRT(W(IQ3))
      ENDDO; ENDDO; ENDDO
      DO NB=1,QO%NBSUB
      WORK(1:NQBND)=QO%PSICOEF(NB,:,NK,ISP); QO%PSICOEF(NB,:,NK,ISP)=0
      DO IQ1=1,NQBND; DO IQ2=1,NQBND
      QO%PSICOEF(NB,IQ1,NK,ISP)=QO%PSICOEF(NB,IQ1,NK,ISP)+SNHK(IQ2,IQ1)*WORK(IQ2)
      ENDDO; ENDDO; ENDDO
      DEALLOCATE(SK,SNHK,WORK,RWORK,W)
      RETURN

      END SUBROUTINE ORTHOK_ALL

!=======================================================================
      SUBROUTINE ORTHOK2(SK,WF,WPRJ)
      USE QO_MODULE; USE prec
      IMPLICIT NONE
      COMPLEX(q) SK(QO%NQBND,QO%NQBND),WF(QO%NRPLWV,QO%NQBND)
      COMPLEX(q) WPRJ(QO%NPROD,QO%NQBND)
! Local
      COMPLEX(q) SNHK(QO%NQBND,QO%NQBND)
      INTEGER IQ1,IQ2,IQ3,NB
      INTEGER LWORK,IERR
      REAL(q) RWORK(3*QO%NQBND),W(QO%NQBND)
      COMPLEX(q) WORK(32*QO%NQBND),ALPHA,BETA
      COMPLEX(q),ALLOCATABLE::CBUF(:,:)

      LWORK=32*QO%NQBND
      CALL ZHEEV('V','U',QO%NQBND,SK,QO%NQBND,W,WORK,LWORK,RWORK,IERR)
      IF(IERR.NE.0)STOP 'Check ORTHOK2!'
      DO IQ3=1,QO%NQBND; IF(W(IQ3).GT.1.E-4_q)CYCLE
      WRITE(0,*)
      WRITE(0,*)'EIGENVALUES:'
      WRITE(0,'(10f8.3)')W
      STOP 'Check ORTHOK2!'
      ENDDO
! Set up S^(-1/2)
      SNHK=0
      DO IQ1=1,QO%NQBND; DO IQ2=1,QO%NQBND; DO IQ3=1,QO%NQBND
        SNHK(IQ2,IQ1)=SNHK(IQ2,IQ1)+SK(IQ2,IQ3)*CONJG(SK(IQ1,IQ3))/SQRT(W(IQ3))
      ENDDO; ENDDO; ENDDO
      ALPHA=1._q; BETA=0
      ALLOCATE(CBUF(QO%NRPLWV,QO%NQBND)); CBUF=WF
      CALL ZGEMM('N','N',QO%NRPLWV,QO%NQBND,QO%NQBND,ALPHA,CBUF,QO%NRPLWV,SNHK,QO%NQBND,BETA,WF,QO%NRPLWV)
      DEALLOCATE(CBUF); ALLOCATE(CBUF(QO%NPROD,QO%NQBND))
      CBUF=WPRJ
      CALL ZGEMM('N','N',QO%NPROD,QO%NQBND,QO%NQBND,ALPHA,CBUF,QO%NPROD,SNHK,QO%NQBND,BETA,WPRJ ,QO%NPROD)
      RETURN

      END SUBROUTINE ORTHOK2

!======================================================================= 
! TB parameter generation.
!=======================================================================
      SUBROUTINE TBPAR(WDES,NONL_S,SYMM,LATT_CUR)
      USE constant; USE wave; USE nonl; USE base; USE paw; USE lattice; USE QO_MODULE; USE prec
      IMPLICIT NONE
      TYPE (wavedes)     WDES          ! descriptor for wavefunction
      TYPE (nonl_struct) NONL_S        ! descriptor for non local part of PP (reciprocal space)
      TYPE (symmetry)    SYMM
      TYPE (latt)        LATT_CUR
! Local
      INTEGER IC,NK,NT1,NT2,NI1,NI2,IQ1,IQ2,L1,L2,M1,M2,ISP,NB,IQ1_,IQ2_
      INTEGER NIONS,NIS1,NIS2,I
      REAL(q) POS1(3),POS2(3),DPOS(3)
      COMPLEX(q),ALLOCATABLE::PHASEK(:,:),CBUF(:)
      REAL(q)   ,ALLOCATABLE::RBUF(:)
      
      WRITE(*,'(" TBPAR...")')
! Set up pairdist
      NIONS=NONL_S%NIONS
      ALLOCATE(QO%PAIRDIST(NIONS,NIONS,QO%WS%DIM)); QO%PAIRDIST=0
      ALLOCATE(QO%NINIPR(3,NIONS,NIONS,QO%WS%DIM)); QO%NINIPR  =0
      DO NI1=1,NIONS; POS1=QO%POSREF(:,NI1)
      DO NI2=1,NIONS; POS2=QO%POSREF(:,NI2)
      DO IC=1,QO%WS%DIM
      DPOS=QO%WS%IR(:,IC)+POS2-POS1
      DPOS=MATMUL(QO%LAT%A,DPOS)
      QO%NINIPR(:,NI1,NI2,IC)=DPOS
      QO%PAIRDIST(NI1,NI2,IC)=SQRT(SUM(DPOS**2))
      ENDDO; ENDDO; ENDDO ! IC,NI2,NI1
      CALL CHK_SHELL(QO%PAIRDIST,NIONS,QO%WS%DIM,QOIO)
! CHECK CONTCAR TO SEE UNEXPECTED SHIFT
      DO IC=1,3
        IF(QO%WS%NGRID(IC).GT.1)CYCLE
        NI1=1
        DO NI2=2,NIONS
        IF(ABS(NONL_S%POSION(IC,NI2)-NONL_S%POSION(IC,NI1)).LT..9_q)CYCLE
        WRITE(QOIO,'(" WARNING: MAY NEED TO CHECK CONTCAR AND KPOINTS FOR NIs=",2I3)')NI1,NI2
        ENDDO
      ENDDO ! IC

      ALLOCATE(PHASEK(QO%WS%DIM,QO%NKPTS),RBUF(QO%WS%DIM),CBUF(QO%NKPTS)); PHASEK=0
      DO IC=1,QO%WS%DIM; DO NK=1,QO%NKPTS
        PHASEK(IC,NK)=EXP(-CITPI*SUM(QO%WS%IR(:,IC)*WDES%VKPT(:,NK)) ) ! EXP(iK.R)
      ENDDO; ENDDO ! NK,IC

      DO ISP=1,QO%ISPIN
      IQ1=0; NIS1=0
      DO NT1=1,NONL_S%NTYP;   DO NI1=1,NONL_S%NITYP(NT1); NIS1=NIS1+1
      IQ1_=0
      DO L1=1,QOYS(NT1)%LMAX; IF(.NOT.QOYS(NT1)%LIN(L1 ))CYCLE
      DO M1=1,2*QOYS(NT1)%LPS(L1)+1
      IQ1=IQ1+1; IQ1_=IQ1_+1
      IQ2=0; NIS2=0
      DO NT2=1,NONL_S%NTYP;   DO NI2=1,NONL_S%NITYP(NT2); NIS2=NIS2+1
      IQ2_=0
      DO L2=1,QOYS(NT2)%LMAX; IF(.NOT.QOYS(NT2)%LIN(L2 ))CYCLE
      DO M2=1,2*QOYS(NT2)%LPS(L2)+1
      IQ2=IQ2+1; IQ2_=IQ2_+1
      DO NB=1,QO%NBSUB
      CBUF=WDES%WTKPT(:)*CONJG(QO%PSICOEF(NB,IQ1,:,ISP))*QO%PSICOEF(NB,IQ2,:,ISP)
      RBUF=MATMUL(PHASEK,CBUF)
      QO%TBS(IQ1,IQ2,:,ISP)=QO%TBS(IQ1,IQ2,:,ISP)+RBUF  ! <IQ1|IQ2(R)>
      CBUF=CBUF*QO%CELEN(NB,:,ISP)
      RBUF=MATMUL(PHASEK,CBUF)
      QO%TBH(IQ1,IQ2,:,ISP)=QO%TBH(IQ1,IQ2,:,ISP)+RBUF  ! <IQ1|H|IQ2(R)>
      ENDDO  ! NB
      ENDDO; ENDDO; ENDDO; ENDDO ! M2,L2,NI2,NT2
      ENDDO; ENDDO; ENDDO; ENDDO ! M1,L1,NI1,NT1
      ENDDO  ! ISP
      RETURN

      END SUBROUTINE TBPAR

!=======================================================================
! Out out TB parameters
!=======================================================================
      SUBROUTINE OUT_TBHS(MODE)
      USE QO_MODULE; USE prec
      IMPLICIT NONE
      INTEGER MODE
! Local
      INTEGER IC,NK,NT1,NT2,NI1,NI2,IQ1,IQ2,L1,L2,M1,M2,ISP,NB,IQ1_,IQ2_
      INTEGER NIONS,NIS1,NIS2,I,CELLDIM
      INTEGER,PARAMETER :: NC=20
      REAL(q) RC(NC),EC(NC),DIST,ABSE

      WRITE(*,'(" OUT_TBHS...")')
      EC=0
      IF(QO%LWSCEL)THEN
        CELLDIM=-QO%WS%DIM
      ELSE
        CELLDIM= QO%WS%DIM
      ENDIF
      DO I=1,NC; RC(I)=I; ENDDO
      SELECT CASE(MODE)
       CASE(1) 
        OPEN(QOIU5,FILE='TBHC.INP',STATUS='REPLACE')
        OPEN(QOIU6,FILE='TBSC.INP',STATUS='REPLACE')
        WRITE(QOIU5,'(" DOS%NET,DOS%EFERMI:")')
        WRITE(QOIU5,'(2F8.3)')QO%NET,QO%EFERMI
        WRITE(QOIU5,'(" WS%DIM,WS%IR,WS%NDEG:")')
        WRITE(QOIU5,'(I7)')CELLDIM
        DO IC=1,QO%WS%DIM; WRITE(QOIU5,'(4I4)')QO%WS%IR(:,IC),QO%WS%NDEG(IC); ENDDO
        WRITE(QOIU5,'(" QB%DIM,QB%ISPIN,QB%RCUT:")')
        WRITE(QOIU5,'(I4,I2,F7.2)')QO%NQBND,QO%ISPIN,QO%RCUT
        WRITE(QOIU5,'(" QB%NI,QB%LMCHAR:")')
        WRITE(QOIU5,'(10I4)')QO%IQION
        WRITE(QOIU5,'(10(1X,A3))')QO%IQLMCHAR
        WRITE(QOIU5,'(" QB%HR--ISP,IQ1,IQ2,IC,HR:")')
        WRITE(QOIU6,'(" QB%SR--ISP,IQ1,IQ2,IC,HR:")')
      END SELECT

      DO ISP=1,QO%ISPIN
      IQ1=0; NIS1=0
      DO NT1=1,QO%NTYP;   DO NI1=1,QO%NITYP(NT1); NIS1=NIS1+1
      IQ1_=0
      DO L1=1,QOYS(NT1)%LMAX; IF(.NOT.QOYS(NT1)%LIN(L1 ))CYCLE
      DO M1=1,2*QOYS(NT1)%LPS(L1)+1
      IQ1=IQ1+1; IQ1_=IQ1_+1
      IQ2=0; NIS2=0
      DO NT2=1,QO%NTYP;   DO NI2=1,QO%NITYP(NT2); NIS2=NIS2+1
      IQ2_=0
      DO L2=1,QOYS(NT2)%LMAX; IF(.NOT.QOYS(NT2)%LIN(L2 ))CYCLE
      DO M2=1,2*QOYS(NT2)%LPS(L2)+1
      IQ2=IQ2+1; IQ2_=IQ2_+1
      DO IC=1,QO%WS%DIM
      DIST=QO%PAIRDIST(NIS1,NIS2,IC)
! TO DELETE THE EFFECT OF THE PERIODIC IMAGE! (QO IS ALSO PERIODIC OVER BV SPACE!)
      IF(DIST.GT.QO%RBV)THEN
        QO%TBS(IQ1,IQ2,IC,ISP)=0; QO%TBH(IQ1,IQ2,IC,ISP)=0
      ENDIF
      IF(MODE.GT.0)THEN
        IF(DIST.GT.QO%RCUT)CYCLE
        WRITE(QOIU5,'(3I4,I6,F12.6)')ISP,IQ1,IQ2,IC,QO%TBH(IQ1,IQ2,IC,ISP)
        WRITE(QOIU6,'(3I4,I6,F12.6)')ISP,IQ1,IQ2,IC,QO%TBS(IQ1,IQ2,IC,ISP)
      ENDIF
      ABSE=ABS(QO%TBH(IQ1,IQ2,IC,ISP))
      DO I=1,NC; IF(DIST.GT.RC(I)) EC(I)=MAX(EC(I),ABSE); ENDDO
      ENDDO  ! IC
      ENDDO; ENDDO; ENDDO; ENDDO
      ENDDO; ENDDO; ENDDO; ENDDO ! NT,NI,LL,MM
      ENDDO  ! ISP
      IF(MODE.GT.0)RETURN
      WRITE(*,'(" RCUT          ECUT")')
      WRITE(QOIO,'(" RCUT          ECUT")')
      DO I=1,NC
        WRITE(*,'(2F8.4)')RC(I),EC(I)
        WRITE(QOIO,'(2F8.4)')RC(I),EC(I)
      ENDDO
      WRITE(*,'(" Current rcut=",F8.4)')QO%RCUT
      WRITE(*,'(" Enter the adjusted rcut:(Better if Ecut<0.01eV)")')
!      READ(*,*)QO%RCUT
      WRITE(QOIO,'(" Adjusted rcut=",F8.4)')QO%RCUT  
      RETURN

      END SUBROUTINE OUT_TBHS

!======================================================================= 
! Test TB parameters.
!=======================================================================
      SUBROUTINE TBTEST(WDES,W)
      USE wave; USE QO_MODULE; USE prec
      IMPLICIT NONE
      TYPE (wavedes)     WDES          ! descriptor for wavefunction
      TYPE (wavespin)    W             ! array for wavefunction
! Local 
      INTEGER NK,NB,ISP,IDE(QO%NQBND),NKT
      REAL(q) EVAL(QO%NQBND),ETMP
      INTEGER            :: W_NBANDS
      COMPLEX(q),POINTER :: W_CELEN(:,:,:)

      WRITE(*,'(" TBTEST...")')
      NKT=MIN(7,QO%NKPTS)
      W_NBANDS=WDES%NBANDS-QO%NBMIN+1
      W_CELEN=>W%CELEN(QO%NBMIN:,:,:)
      WRITE(QOIO,*)
      WRITE(QOIO,*)'TB TEST: EQOLDA, EQOTB, ELDA'
      DO ISP=1,QO%ISPIN; WRITE(QOIO,'(A5,I2)')'&ISP=',ISP
      DO NK=1,NKT
        CALL TBEK(ISP,WDES%VKPT(:,NK),EVAL,0)
        WRITE(QOIO,'(A1,I3,3f7.3)')'&',NK,WDES%VKPT(:,NK)
        CALL ORDER(QO%NQBND,QO%CELEN(:,NK,ISP),IDE)
        DO NB=1,QO%NQBND
        IF(NB.LE.W_NBANDS)THEN
          ETMP=REAL(W_CELEN(NB,NK,ISP),q)
        ELSE
          ETMP=0._q
        ENDIF
        WRITE(QOIO,'(I5,1x,3f10.3)')NK,QO%CELEN(IDE(NB),NK,ISP),EVAL(NB),ETMP
        ENDDO
      ENDDO; ENDDO

      END SUBROUTINE TBTEST

!======================================================================= 
! Get TB Band structure 
!=======================================================================
      SUBROUTINE TBBAND(LATT_CUR)
      USE lattice; USE QO_MODULE; USE prec
      IMPLICIT NONE
      TYPE (latt)        LATT_CUR
! Local
      INTEGER       NKPTS,NKP,NINTER,INDEX,N,IERR,ISP,MODE
      CHARACTER*1   CSEL,CLINE
      REAL(q)       VKPT(3,QONKDIM),VKPT2(3,QONKDIM),SHIFT(3)
      REAL(q)       EVAL(QO%NQBND)
      REAL(q)       KLEN,DKLEN

      WRITE(*,'(" TBBABD...")')
      OPEN(QOIU3,FILE='KPOINTS',STATUS='OLD')
      DO 
      READ(QOIU3,'(A1)',ERR=70111,END=70111)CSEL
      IF(CSEL=='Q')EXIT
      ENDDO
      READ(QOIU3,*)
      READ(QOIU3,*)NKPTS
      READ(QOIU3,'(A1)',ERR=70111,END=70111)CSEL
      IF (CSEL/='L'.AND.CSEL/='l')GOTO 70111
      READ(QOIU3,'(A1)',ERR=70111,END=70111)CSEL
      VKPT2=0
      NINTER=NKPTS
      NKP=0 ! counter for the number of k-points already read in
      DO
        NKP=NKP+1
        IF (NKP>QONKDIM) THEN
          WRITE(*,*)'ERROR: increase QONKDIM!'
          STOP
        ENDIF
        READ(QOIU3,*,IOSTAT=IERR) VKPT2(1,NKP),VKPT2(2,NKP),VKPT2(3,NKP)
        IF (IERR/=0) EXIT
      ENDDO
      NKPTS=NKP-1
      CLOSE(QOIU3)
      IF (CSEL=='K'.OR.CSEL=='k'.OR.CSEL=='C'.OR.CSEL=='c') THEN
      VKPT2(:,1:NKPTS)=VKPT2(:,1:NKPTS)/LATT_CUR%SCALE
      CALL KARDIR(NKPTS,VKPT2,LATT_CUR%A)
      ENDIF
      INDEX=0
      ! make NKPTS even
      NKPTS=(NKPTS/2)*2; VKPT=0
      DO NKP=1,NKPTS-1,2
        SHIFT=(VKPT2(:,NKP+1)-VKPT2(:,NKP))/(NINTER-1)
        DO N=0,NINTER-1
        INDEX=INDEX+1
        IF (INDEX>QONKDIM) THEN
          WRITE(*,*)'ERROR: increase QONKDIM'
          STOP
        ENDIF
        VKPT(:,INDEX)=VKPT2(:,NKP)+SHIFT*N
        ENDDO
      ENDDO
      NKPTS=INDEX

      WRITE(QOIO,*)'TB BAND:'
      DO ISP=1,QO%ISPIN
      WRITE(CSEL,'(I1)')ISP
      OPEN(QOIU3,FILE='TBBAND'//CSEL,STATUS='REPLACE')
      KLEN=0
      DO NKP=1,NKPTS
        IF(NKP.EQ.1)THEN; MODE=1; ELSE; MODE=0; ENDIF
        CALL TBEK(ISP,VKPT(:,NKP),EVAL,MODE)
        IF(NKP.EQ.1)THEN
          DKLEN=0
        ELSE
          CALL CALC_DLEN(DKLEN,VKPT(:,NKP-1),VKPT(:,NKP),LATT_CUR%B,LATT_CUR%SCALE)
        ENDIF
        KLEN=KLEN+DKLEN
        WRITE(QOIU3,70110)KLEN,EVAL -QO%EFERMI
      ENDDO
      CLOSE(QOIU3)
      ENDDO
      RETURN
70110 FORMAT(<QO%NQBND+1>f10.3)
70111 CONTINUE
      WRITE(*,*)'CHECK "QO" IN KPOINTS, RESET KPOINTS!'
      RETURN

      END SUBROUTINE TBBAND

!======================================================================= 
! Get TB eigen values @ one K-point.
!=======================================================================
      SUBROUTINE TBEK(ISP,KPT,EVAL,MODE)
      USE constant; USE QO_MODULE; USE prec
      IMPLICIT NONE
      INTEGER ISP,MODE
      REAL(q) KPT(3),EVAL(QO%NQBND)
! Local 
      COMPLEX(q) PHASEK
      INTEGER IC,IQ1,IQ2,LWORK,IERR,NBQ(QO%NQBND),IQR1,IQR2
      COMPLEX(q) TBHK (QO%NQBND,QO%NQBND),TBSK (QO%NQBND,QO%NQBND)
      REAL(q)    RWORK(3*QO%NQBND),WEIGHT
      CHARACTER JOBZ
      COMPLEX(q) WORK (32*QO%NQBND)
      
      PHASEK=0; EVAL=0; TBHK=0; TBSK=0
      DO IC=1,QO%WS%DIM
      PHASEK=EXP(-CITPI*SUM(QO%WS%IR(:,IC)*KPT(:)) ) ! EXP(-iK.R)
      DO IQ1=1,QO%NQBND; DO IQ2=IQ1,QO%NQBND
      IF(QO%PAIRDIST(QO%IQION(IQ1),QO%IQION(IQ2),IC).GT.QO%RCUT)CYCLE
      TBHK(IQ1,IQ2)=TBHK(IQ1,IQ2)+QO%TBH(IQ1,IQ2,IC,ISP)*PHASEK/QO%WS%NDEG(IC)
      TBSK(IQ1,IQ2)=TBSK(IQ1,IQ2)+QO%TBS(IQ1,IQ2,IC,ISP)*PHASEK/QO%WS%NDEG(IC)
      ENDDO; ENDDO
      ENDDO
      LWORK=32*QO%NQBND
      IF(MODE.EQ.1)THEN; JOBZ='V'; ELSE; JOBZ='N'; ENDIF
      CALL ZHEGV(1,JOBZ,'U',QO%NQBND,TBHK,QO%NQBND,TBSK,QO%NQBND,EVAL,WORK, &
                &LWORK,RWORK,IERR)
      IF(IERR.NE.0)THEN
        WRITE(*,*)'IERR:',IERR,QO%NQBND
        STOP ' Check TBEKI!'
      ENDIF
      IF(MODE.EQ.0)RETURN
      WRITE(QOIO,'(" KPT=",3F8.3)')KPT
      NBQ=1
      DO IQ1=1,QO%NQBND
      WEIGHT=ABS(TBHK(1,IQ1))
      DO IQ2=2,QO%NQBND
      IF(WEIGHT.GE.ABS(TBHK(IQ2,IQ1)))CYCLE
      NBQ(IQ1)=IQ2; WEIGHT=ABS(TBHK(IQ2,IQ1))
      ENDDO; ENDDO
      WRITE(QOIO,'(" NB  IQ   WEIGHT")')
      DO IQ1=1,QO%NQBND
      WRITE(QOIO,'(I4,I4,F8.3)')IQ1,NBQ(IQ1),ABS(TBHK(NBQ(IQ1),IQ1))
      ENDDO

      END SUBROUTINE TBEK

!======================================================================= 
! Radial Fourier transformation of RW to get QW. GMAX=P(NT)%PSMAXN
!=======================================================================
      SUBROUTINE QOGREP(P,NT)
      USE pseudo; USE QO_MODULE; USE prec
      IMPLICIT NONE
      INTEGER NT
      TYPE(potcar) P
! Local
      INTEGER LMAX,NMAX,IG,IR,LL,ISGN
      REAL(q) GLEN,RLEN,GSTEP,BJ
      REAL(q),POINTER :: AUX(:)

      LMAX=QOYS(NT)%LMAX
! Follow the recipe for PSPRNL.
      NMAX=QOYS(NT)%R%NMAX
      ALLOCATE(AUX(NMAX)); AUX=0 
      GSTEP=P%PSMAXN/NPSNL
      GLEN=0; ISGN=1

      DO LL=1,LMAX
      DO IG=1,NPSNL
      GLEN=(IG-1)*GSTEP
      DO IR=1,NMAX
        RLEN=QOYS(NT)%R%R(IR)
        CALL SBESSEL(GLEN*RLEN,BJ,QOYS(NT)%LPS(LL))
        AUX(IR)=RLEN*QOYS(NT)%RW(IR,LL)*BJ
      ENDDO ! IR
      QOYS(NT)%QW(IG,LL)=SUM(AUX*QOYS(NT)%R%SI)
      ENDDO ! IG
      QOYS(NT)%QW(0,LL)=QOYS(NT)%QW(2,LL)*ISGN
      ISGN=-ISGN
      ENDDO ! LL

      END SUBROUTINE QOGREP


!=======================================================================
! SETUP REDIAL PART OF QO
!=======================================================================
      SUBROUTINE QONORM(NT,P,UNITS)
      USE pseudo; USE constant; USE QO_MODULE; USE prec
      IMPLICIT NONE
      INTEGER NT
      TYPE (potcar)      P
      CHARACTER UNITS*1 ! B/b: Bohr unit; A: anstrong, AE; a: PS
      REAL(q),EXTERNAL::GET_RNORM
! Local 
      INTEGER LL,LMAX,NMAX,IR
      REAL(q) NORM0

! UNIT conversion.
! Logrithmic grid is a bit different here, R(i)=RSTART * (exp [H (i-1)]-1) 
      NMAX=QOYS(NT)%R%NMAX
      IF(UNITS.EQ.'B'.OR.UNITS.EQ.'b')THEN
      QOYS(NT)%R%R     =QOYS(NT)%R%R*AUTOA
      QOYS(NT)%RW      =QOYS(NT)%RW /SQRT(AUTOA)
      ENDIF
      QOYS(NT)%R%RMAX  =QOYS(NT)%R%R(NMAX)
      QOYS(NT)%R%RSTART=QOYS(NT)%R%R(1)
      QOYS(NT)%R%REND  =QOYS(NT)%R%R(NMAX)
      QOYS(NT)%R%D     =(QOYS(NT)%R%R(NMAX  )-QOYS(NT)%R%R(NMAX-1))/ &
                       &(QOYS(NT)%R%R(NMAX-1)-QOYS(NT)%R%R(NMAX-2))
      QOYS(NT)%R%H     =LOG(QOYS(NT)%R%D)

      IF(ASSOCIATED(P%R%R))THEN
      IF(UNITS.EQ.'A'.OR.UNITS.EQ.'a')THEN
      IF(ABS((QOYS(NT)%R%H-P%R%H)/P%R%H).GT.1.E-4)THEN
      STOP 'MISMATCH RADIAL GRID!'
      ELSE
      QOYS(NT)%R%D=P%R%D; QOYS(NT)%R%H=P%R%H
      ENDIF
      ENDIF
! UPTO NOW, QOYS(NT)%RW IS ACTUALLY AE WAVE, (I USE CUT AND PASTE TO GET THE PSEUDO ONE.)
      IF(UNITS.EQ.'A') CALL AE2PS(NT,P)
      ENDIF

      LMAX=QOYS(NT)%LMAX
! Apply scaling facotr, alpha's and renormalized to the original norm/=1(charge deficiency).
      CALL SET_SIMP(QOYS(NT)%R)
      WRITE(QOIO,'(" NT=",I2)')NT
! Check original norm
      DO LL=1,LMAX
        NORM0=GET_RNORM(QOYS(NT)%LPS(LL),LL,P,NT,1)
        WRITE(QOIO,'(a3,I2,2(2x,a6,f6.3))')' L=',QOYS(NT)%LPS(LL),'NORM0=',NORM0
        QOYS(NT)%RW(:,LL)=QOYS(NT)%RW(:,LL)*EXP(-QOYS(NT)%ALPHA(LL)*QOYS(NT)%R%R)
        NORM0=GET_RNORM(QOYS(NT)%LPS(LL),LL,P,NT,1)
        WRITE(QOIO,'(a3,I2,2(2x,a6,f6.3))')' L=',QOYS(NT)%LPS(LL),'NORM1=',NORM0
        QOYS(NT)%RW(:,LL)=QOYS(NT)%RW(:,LL)*SQRT(ABS(1/NORM0))
! Target orbital might be negative due to output.
        NORM0=SUM(QOYS(NT)%R%SI(:)*QOYS(NT)%RW(:,LL))
        IF(NORM0.LT.0) QOYS(NT)%RW(:,LL)=-QOYS(NT)%RW(:,LL)
      ENDDO ! LL
      RETURN

      END SUBROUTINE QONORM


!=======================================================================
! CUT AND PASTE TO GET PS WAVE
!=======================================================================
      SUBROUTINE AE2PS(NT,P)
      USE pseudo; USE QO_MODULE; USE prec
      IMPLICIT NONE 
      TYPE (potcar)      P
      INTEGER NT
! LOCAL
      INTEGER LL,CH,IND,LL2CH
      INTEGER LMAX,NMAX
      REAL(q) FAKT,FAKT0
 
      LMAX=QOYS(NT)%LMAX; NMAX=P%R%NMAX
      DO LL=1,LMAX; FAKT0=0
      DO CH=1,P%LMAX
      IF(P%LPS(CH).NE.QOYS(NT)%LPS(LL))CYCLE
      FAKT=QOYS(NT)%RW(NMAX,LL)/P%WAE(NMAX,CH)/(QOYS(NT)%RW(NMAX-3,LL)/P%WAE(NMAX-3,CH))
      IF(ABS(ABS(FAKT)-1).LT.ABS(FAKT0-1))THEN
        FAKT0=ABS(FAKT)
        LL2CH=CH
      ENDIF
      ENDDO ! CH
      WRITE(QOIO,'(" WAE MATCHING FAKT0=",F8.2)')FAKT0
      IF(FAKT0.LT..1_q)CYCLE
      FAKT=QOYS(NT)%RW(NMAX,LL)/P%WPS(NMAX,LL2CH)
      QOYS(NT)%RW(1:NMAX,LL)=P%WPS(1:NMAX,LL2CH)*FAKT
      ENDDO

      END SUBROUTINE AE2PS


!=======================================================================
! Get <AO|S|AO> for radial part. L=0,1,2...
!=======================================================================
      FUNCTION GET_RNORM(L,IL,P,NT,MODE)
      USE pseudo; USE constant; USE QO_MODULE; USE prec
      IMPLICIT NONE
      INTEGER L,IL,NT,MODE
      TYPE (potcar)      P 
      REAL(q) GET_RNORM
! Local
      INTEGER CH,IND,I,IS,J,CHS(2)
      REAL(q) ARGSC,DIST,REM,VPS,SNORM(2)

      IS=0; SNORM=0; CHS=1
      ARGSC=NPSRNL/P%PSRMAX
      DO CH=1,P%LMAX
      IF(P%LPS(CH).NE.L)CYCLE
      IS=IS+1; CHS(IS)=CH
      IF(MODE.EQ.1.AND.IS.EQ.1)THEN
      WRITE(80,'("# WPS NT=",I2," L=",I2)')NT,L
      DO IND=1,P%R%NMAX
        WRITE(80,*)P%R%R(IND),P%WPS(IND,CH)
      ENDDO
      WRITE(80,*)
      WRITE(80,'("# WAE NT=",I2," L=",I2)')NT,L
      DO IND=1,P%R%NMAX
        WRITE(80,*)P%R%R(IND),P%WAE(IND,CH) 
      ENDDO 
      WRITE(80,*)
      ENDIF ! MODE.EQ.1.AND.IS.EQ.1
      IF(IS.GT.2)STOP'IS>2!'
      DO IND=1,QOYS(NT)%R%NMAX
      DIST=QOYS(NT)%R%R(IND)
      IF(DIST.GT.P%PSRMAX)EXIT
      I=MIN(INT(DIST*ARGSC)+1,NPSRNL-1)
      REM=DIST-P%PSPRNL(I,1,CH)
      VPS=P%PSPRNL(I,2,CH)+REM*(P%PSPRNL(I,3,CH)+REM*(P%PSPRNL(I,4,CH)+REM*P%PSPRNL(I,5,CH)))
! Ahh: PSPRNL is only the R_l(r) instead of rR_l(r)
      SNORM(IS)=SNORM(IS)+QOYS(NT)%R%SI(IND)*VPS*QOYS(NT)%RW(IND,IL)*DIST
      ENDDO
      ENDDO ! CH

      IF(MODE.EQ.1)THEN
      WRITE(80,'("# UNSCALED RW NT=",I2," L=",I2)')NT,L
      DO IND=1,QOYS(NT)%R%NMAX
        WRITE(80,*)QOYS(NT)%R%R(IND),QOYS(NT)%RW(IND,IL)
      ENDDO
      WRITE(80,*)
      ENDIF ! MODE.EQ.1

      GET_RNORM=SUM(QOYS(NT)%R%SI(:)*QOYS(NT)%RW(:,IL)**2)
      DO I=1,2; DO J=1,2
      GET_RNORM=GET_RNORM+SNORM(I)*SNORM(J)*P%QION(CHS(I),CHS(J))
      ENDDO; ENDDO

      END FUNCTION GET_RNORM

!=======================================================================
! GENERATE LCAO <AO| for one k point and calc their CPROJ
!=======================================================================
      SUBROUTINE GEN_LCAO(LATT_CUR,WDES,NONL_S,NK)
      USE lattice; USE wave; USE nonl; USE constant; USE asa; USE QO_MODULE; USE prec
      IMPLICIT NONE
      TYPE (latt)        LATT_CUR
      TYPE (nonl_struct) NONL_S        ! descriptor for non local part of PP (reciprocal space)
      TYPE (wavedes)     WDES          ! descriptor for wavefunction
      INTEGER NK
! Local
      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W1
      REAL(q) KG(3),KGX,KGY,KGZ,QW1,FAKT
      COMPLEX(q)         CFKT
      INTEGER NT,NI,NIS,IQO,LM,LMQ,NPL,IPL,MM,LL,NITYP,NIONS,LMQD
      REAL(q),ALLOCATABLE:: YLM(:,:),KGNORM(:)
      COMPLEX(q),ALLOCATABLE:: KGREXP(:,:) ! EXP(i(K+G)R) for atoms with same type.

!=======================================================================
! number of G-Vectors on the k-Point
!=======================================================================
      FAKT= 1/SQRT(LATT_CUR%OMEGA)
      NPL= WDES%NGVECTOR(NK)

      ALLOCATE(YLM(NPL,16),KGNORM(NPL)); YLM=0; KGNORM=0 ! maxL=3

! Set up phase and YLM 
      NIONS=NONL_S%NIONS
      ALLOCATE(KGREXP(NPL,NIONS)); KGREXP=0

      DO IPL=1,NPL
      KGX=WDES%IGX(IPL,NK)+WDES%VKPT(1,NK); KG(1)=KGX  
      KGY=WDES%IGY(IPL,NK)+WDES%VKPT(2,NK); KG(2)=KGY  
      KGZ=WDES%IGZ(IPL,NK)+WDES%VKPT(3,NK); KG(3)=KGZ  
      DO NI=1,NIONS
        KGREXP(IPL,NI)=EXP(CITPI*SUM(KG(:)*NONL_S%POSION(:,NI))) 
      ENDDO ! NI
      KG(:)=(KGX*LATT_CUR%B(:,1)+KGY*LATT_CUR%B(:,2)+KGZ*LATT_CUR%B(:,3))*TPI
      KGNORM(IPL)=MAX(SQRT(SUM(KG(:)**2)),1E-10_q)
      KG=KG/KGNORM(IPL)  
      CALL SETYLM(3,1,YLM(IPL:IPL,:),KG(1:1),KG(2:2),KG(3:3)) 
      ENDDO ! IPL 

! Begin generating LCAO and set IONIQ
      NIS=0; IQO=0
      DO NT=1,NONL_S%NTYP
      NITYP=NONL_S%NITYP(NT); LMQD=QOYS(NT)%LMQD

! Set IONIQ
      IF(NK.EQ.1)THEN
      DO NI=1,NITYP; DO MM=1,LMQD
      QO%IONIQ(MM,NIS+NI)=IQO+(NI-1)*LMQD+MM
      QO%IQION(IQO+(NI-1)*LMQD+MM)=NIS+NI
      QO%IQLMCHAR(IQO+(NI-1)*LMQD+MM)=QOYS(NT)%LMCHAR(MM)
      ENDDO; ENDDO ! NI,MM
      ENDIF
      
      LMQ=0
      DO LL=1,QOYS(NT)%LMAX
      IF(.NOT.QOYS(NT)%LIN(LL))CYCLE
      CFKT=4*PI*(0._q,1._q)**QOYS(NT)%LPS(LL)*FAKT
      DO IPL=1,NPL
      CALL INTERPOLATE4(QOYS(NT)%PSMAXN,QO%NPSNL,QOYS(NT)%QW(:,LL),KGNORM(IPL),QW1)
      DO MM=1,2*QOYS(NT)%LPS(LL)+1
      LM=MM+QOYS(NT)%LPS(LL)**2
      DO NI=1,NITYP
        QO%WATOM(IPL,IQO+LMQ+MM+(NI-1)*LMQD)=CONJG(QW1*YLM(IPL,LM)*KGREXP(IPL,NIS+NI)*CFKT)  ! |LCAO>
      ENDDO ! I
      ENDDO; ENDDO ! MM,IPL
      LMQ=LMQ+2*QOYS(NT)%LPS(LL)+1
      ENDDO ! LL
      NIS=NIS+NITYP; IQO=IQO+NITYP*LMQD
      ENDDO ! NT

! Calc CPROJ
      CALL PHASE(WDES,NONL_S,NK)

! Set phase for NK. ONly exp(igr) since exp(ikr) would be cancelled considering |Pi>Qij<Pj|
      CALL SETWDES(WDES,WDES1,NK)

      DO IQO=1,QO%NQBND
      W1%CPTWFP=>QO%WATOM(:,IQO)
      W1%CPROJ =>QO%WAPRJ(:,IQO)
      CALL PROJ1(NONL_S,WDES1,W1)  ! Calc <P~|LCAO>
      ENDDO
      RETURN

      END SUBROUTINE GEN_LCAO

!++++++++++++++++++++++++++++++++++++++++++++++++++
! Given F(NMAX) and the maximum x, xmax, four-point 
! interpolation of f @ x. used in VASP
!++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE INTERPOLATE4(XMAX,NMAX,F,X,FX)
      USE prec
      IMPLICIT NONE
      INTEGER NMAX
      REAL(q) XMAX,F(0:NMAX),X,FX
! Local
      REAL(q) ARGSC,ARG,V1,V2,V3,V4,T0,T1,T2,T3,REM
      INTEGER NADDR

      ARGSC=NMAX/XMAX
      ARG=X*ARGSC+1
      NADDR=INT(ARG)
      IF(NADDR<NMAX-2)THEN
        REM=MOD(ARG,1.0_q)
        V1=F(NADDR-1)
        V2=F(NADDR  )
        V3=F(NADDR+1)
        V4=F(NADDR+2)
        T0=V2
        T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
        T2=(V1+V3-(2*V2))/2._q
        T3=(V4-V1+(3*(V2-V3)))/6._q
        FX= T0+REM*(T1+REM*(T2+REM*T3))
      ELSE
        FX=0
      ENDIF

      END SUBROUTINE INTERPOLATE4

!++++++++++++++++++++++++++++++++++++++++++++++++++
! Calc <WFP|1+Q|WF>
!++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE WPSW(WFP,WF,CPROJ,CPROW,NPL,NPRO,RES) 
      USE prec
      IMPLICIT NONE
      INTEGER    NPL,NPRO
      COMPLEX(q) WF(*),WFP(*),CPROJ(*),CPROW(*),RES
      
      RES=SUM( WF(1:NPL)*CONJG(WFP(1:NPL)) )   ! <WF'|WF> 
      RES=RES+DOT_PRODUCT(CPROJ(1:NPRO),CPROW(1:NPRO))   ! Q included in CPROW
      
      END SUBROUTINE WPSW

!=============================================================================
! FIND THE INVERSE OF A GENERAL MATRIX A. ON EXIT: A -> A^-1
!=============================================================================
      SUBROUTINE ZGETRI_(A,LDA)
      USE prec
      IMPLICIT NONE
      INTEGER LDA
      COMPLEX(q) A(LDA,LDA)
! Local
      INTEGER IPIV(LDA),INFO
      COMPLEX(q) WORK(LDA)

      CALL ZGETRF( LDA, LDA, A, LDA, IPIV, INFO )
      IF (INFO.NE.0)THEN
      WRITE(0,*)'INFO=',INFO; STOP 'DGETRI_:DGETRF ERROR!'
      ENDIF
      CALL ZGETRI( LDA, A, LDA, IPIV, WORK, LDA, INFO )
      IF (INFO.NE.0)THEN
      WRITE(0,*)'INFO=',INFO; STOP 'DGETRI_:DGETRI ERROR!'
      ENDIF

      END SUBROUTINE ZGETRI_

      SUBROUTINE DGETRI_(A,LDA)
      USE prec
      IMPLICIT NONE
      INTEGER LDA
      REAL(q) A(LDA,LDA)
! Local
      INTEGER IPIV(LDA),INFO
      REAL(q) WORK(LDA)

      CALL DGETRF( LDA, LDA, A, LDA, IPIV, INFO )
      IF (INFO.NE.0)THEN
      WRITE(0,*)'INFO=',INFO; STOP 'DGETRI_:DGETRF ERROR!'
      ENDIF
      CALL DGETRI( LDA, A, LDA, IPIV, WORK, LDA, INFO )
      IF (INFO.NE.0)THEN
      WRITE(0,*)'INFO=',INFO; STOP 'DGETRI_:DGETRI ERROR!'
      ENDIF

      END SUBROUTINE DGETRI_

!=============================================================================
! GENERATE A UNIT MATRIX
!=============================================================================
      SUBROUTINE MATUNIT(MAT,NDIM)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,I
      REAL(q) MAT(NDIM,NDIM)

      MAT=0
      DO I=1,NDIM
        MAT(I,I)=1._q
      ENDDO

      END SUBROUTINE MATUNIT

!=============================================================================
      SUBROUTINE CHK_SHELL(PAIRDIST,NIONS,CELLMAX,IO)
      USE prec
      IMPLICIT NONE
      INTEGER NIONS,CELLMAX,IO
      REAL(q) PAIRDIST(NIONS,NIONS,CELLMAX)
! LOCAL
      INTEGER,PARAMETER :: NNN=40
      INTEGER NI,I
      INTEGER,ALLOCATABLE::IDX(:)
      INTEGER CNN(NNN),INN
      REAL(q) RNN(NNN),DIST
      REAL(q),ALLOCATABLE::R(:)

      WRITE(*,'(" CHK_SHELL...")')
      ALLOCATE(IDX(NIONS*CELLMAX),R(NIONS*CELLMAX))
      DO NI=1,NIONS
        R=RESHAPE(PAIRDIST(NI,:,:),(/NIONS*CELLMAX/))
        CALL ORDER(CELLMAX*NIONS,R,IDX)
        DIST=0; INN=0; CNN=0; RNN=0
        DO I=2,CELLMAX*NIONS
          IF(ABS(R(IDX(I))-DIST).GT..05_q)THEN
            INN=INN+1
            IF(INN.GT.NNN)EXIT
            RNN(INN)=R(IDX(I))
            DIST=R(IDX(I))
          ENDIF
          CNN(INN)=CNN(INN)+1
        ENDDO
        WRITE(*,'(" NI=",I3," INN=",I3)')NI,INN
        WRITE(*,'(20F6.2)')RNN
        WRITE(*,'(20I6  )')CNN
        WRITE(IO,'(" NI=",I3," INN=",I3)')NI,INN
        WRITE(IO,'(20F6.2)')RNN
        WRITE(IO,'(20I6  )')CNN
      ENDDO
      DEALLOCATE(IDX,R)
      RETURN

      END SUBROUTINE CHK_SHELL

!**************** SUBROUTINE CALC_DLEN**********************************
! Calculate the distance between two points in direct/reciprocal coord.
! ) direct lattice      (BASIS must be equal to B reciprocal lattice)
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************
      SUBROUTINE CALC_DLEN(DLEN,A,B,BASIS,A0)
      USE prec
      REAL(q) DLEN,A(3),B(3),BASIS(3,3),A0
! Local
      REAL(q) C(3)

      C=MATMUL(BASIS,A-B)
      DLEN=SQRT(SUM(C**2))*A0 ! IN UNIT OF 2pi/scale
      RETURN

      END SUBROUTINE CALC_DLEN
