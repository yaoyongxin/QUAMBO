      MODULE PREC
      INTEGER, PARAMETER :: q = 8
      END MODULE PREC


      MODULE CONSTANT
      USE prec
      REAL(q),PARAMETER  :: PI =3.141592653589793238_q,TPI=2*PI
      COMPLEX(q),PARAMETER  :: CI=(0._q,1._q)
      COMPLEX(q),PARAMETER  :: CITPI = (0._q,1._q)*TPI
      END MODULE CONSTANT

      
      MODULE TB_MODULE
      USE prec
      IMPLICIT NONE

! TYPES FOR TB
      TYPE latt
         REAL(q) :: SCALE
         REAL(q) :: A(3,3),B(3,3)
         REAL(q) :: ANORM(3),BNORM(3)
         REAL(q) :: OMEGA
      END TYPE latt

      TYPE type_info
!only T_INFO
        CHARACTER*40 SZNAM2           ! name of poscar file
        INTEGER NTYP                  ! number of types
        INTEGER NIONS                 ! actual number of ions
        LOGICAL LDIRCO                ! positions in direct/recproc. lattice
        REAL(q), POINTER :: POSION(:,:)  ! positions usually same as DYN%POSION
        INTEGER, POINTER :: ITYP(:)   ! type for each ion
        INTEGER, POINTER :: NITYP(:)  ! number of ions for each type
      END TYPE type_info

      TYPE kpoints_struct
!only  KPOINTS,KPT_BAND,KPT_FS
        INTEGER NKPTS,NKPT(3)        ! actual number of k-points
        REAL(q),POINTER:: VKPT(:,:)  ! coordinate of k-point
        REAL(q),POINTER:: WTKPT(:)   ! symmetry weight-factor for each k-point
        CHARACTER*40  SZNAMK         ! name of k-points file
      END TYPE kpoints_struct

      TYPE dos_struct
!only  DOS
        REAL(q) EFERMI,NET,EFERMIN
        INTEGER ISMEAR               ! type of smearing
        REAL(q) SIGMA                ! type of smearing
        INTEGER NGRID
        REAL(q) EMIN                 ! minimal E for DOS
        REAL(q) EMAX                 ! maximal E for DOS
        REAL(q),POINTER :: T(:,:),I(:,:),PAR(:,:,:)
      END TYPE dos_struct

      TYPE wavespin
!only W
        INTEGER NBANDS,ISPIN,RSPIN
        COMPLEX(q),POINTER:: CPTWFP(:,:,:,:) ! wavefunction
        REAL(q),   POINTER:: FERWE(:,:,:)    ! fermi-weight for each band
        REAL(q),   POINTER:: CELEN(:,:,:)    ! eigenvalues
        REAL(q),   POINTER:: PAR  (:,:,:,:)  ! projected weight
      END TYPE wavespin

      TYPE QO_TYPE
        INTEGER ISPIN  ! ISPIN FOR TBH/S
        INTEGER DIM
        REAL(q) RCUT,MAXRCUT
        CHARACTER(3),POINTER::LMCHAR(:)
        INTEGER,POINTER :: NI(:)
        REAL(q),POINTER :: SR(:,:,:,:)
        OVERLAP,POINTER :: HR(:,:,:,:)
      END TYPE QO_TYPE

      TYPE CELL_TYPE
        INTEGER DIM,ORIG
        INTEGER,ALLOCATABLE :: R(:,:)
      END TYPE CELL_TYPE

      TYPE CTRL_TYPE
        LOGICAL LMULLIKEN,LROBY,LMAYER,LSPDM
      END TYPE CTRL_TYPE

      TYPE SHELL_TYPE
        INTEGER N  ! NUMBER OF ATOM IN THIS SHELL
        REAL(q) DIST
        INTEGER,POINTER :: NI(:) ! ATOM INDEX
        INTEGER,POINTER :: IC(:) ! PRIMITIVE CELL TRANSLATIONAL VECTOR INDEX
      END TYPE SHELL_TYPE

      TYPE ATOM_TYPE
        INTEGER NSL
        TYPE (SHELL_TYPE),ALLOCATABLE::SL(:)
        INTEGER IQMIN,IQMAX,NQ
      END TYPE ATOM_TYPE

! VARIABLES
      TYPE (latt)          ,SAVE :: LATT_CUR
      TYPE (type_info)     ,SAVE :: T_INFO
      TYPE (kpoints_struct),SAVE :: KPOINTS,KPT_BAND,KPT_FS
      TYPE (dos_struct)    ,SAVE :: DOS
      TYPE (wavespin)      ,SAVE :: W
      TYPE (CELL_TYPE)     ,SAVE :: CELL
      TYPE (QO_TYPE)       ,SAVE :: QB
      TYPE (CTRL_TYPE)     ,SAVE :: CL
      TYPE (ATOM_TYPE),ALLOCATABLE,SAVE :: YZ(:)

      CONTAINS

!***********************************************************************
      SUBROUTINE TB_ANALYSIS(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO

      CALL READ_POSCAR(IU,IO)
      CALL READ_KPOINTS(KPOINTS,'KPOINTS',IU,IO)
      CALL READ_TBHS(IU,IO)
      CALL READ_INCAR(IU,IO)
      CALL SET_ATOM_SHELL(IO)
      CALL SET_ATOM_QO(IO)
      CALL SET_BAND_DOS(IO)
      CALL OUT_TBHS(IO)
      CALL CALC_FULL_BAND(.FALSE.,IO)
      CALL SET_FERMI_WT(IO)
      CALL OUT_DOS(IU,IO)
      CALL CALC_BAND_STRU(.FALSE.,IU,IO)
      CALL SPDM_MULLIKEN_MAYER_ANALYSIS(IO)
      CALL ROBY_COV_ION_BO_ANALYSIS(IO)
      RETURN

      END SUBROUTINE TB_ANALYSIS

      SUBROUTINE READ_POSCAR(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      CHARACTER*255  INPLIN,INPWRK
      INTEGER, EXTERNAL :: NITEMS
      CHARACTER*1  CSEL
      INTEGER I,NT,NI,NSCALE
      REAL(q) SCALEX,SCALEY,SCALEZ
      LOGICAL LREF

      WRITE(IO,'(" READING POSCAR...")')
      WRITE(*,'(" READING POSCAR...")')
      INQUIRE(FILE='POSREF',EXIST=LREF)
      IF(LREF)THEN
       OPEN(IU,FILE='POSREF',STATUS='OLD')
      ELSE
       WRITE(IO,'( " POSREF FILE NOT EXIST, USE POSCAR!")')
       OPEN(IU,FILE='POSCAR',STATUS='OLD')
      ENDIF
      READ(IU,'(A40)') T_INFO%SZNAM2
      WRITE(IO,*)T_INFO%SZNAM2
! one scaling parameter or three
      READ(IU,'(A)') INPLIN
      NSCALE=NITEMS(INPLIN,INPWRK,.TRUE.,'F')
      IF (NSCALE==1) THEN
        READ(INPLIN,*) LATT_CUR%SCALE
        SCALEX=1; SCALEY=1; SCALEZ=1
      ELSE IF (NSCALE==3) THEN
        LATT_CUR%SCALE=1
        READ(INPLIN,*) SCALEX,SCALEY,SCALEZ
      ELSE
        WRITE(0,*)'ERROR: there must be 1 or 3 items on line 2 of POSCAR'
        STOP
      ENDIF
      DO I=1,3
        READ(IU,*) LATT_CUR%A(1,I),LATT_CUR%A(2,I),LATT_CUR%A(3,I)
      ENDDO
      IF (LATT_CUR%SCALE<0._q) THEN
!----alternatively give a volume (=abs(scale)) and adjust the lengths of
!----the three lattice vectors to get the correct desired volume ... :
         CALL LATTIC(LATT_CUR)
         LATT_CUR%SCALE=(ABS(LATT_CUR%SCALE)  &
     &                 / ABS(LATT_CUR%OMEGA))**(1._q/3._q)
      ENDIF

      LATT_CUR%A(1,:) =LATT_CUR%A(1,:)*SCALEX*LATT_CUR%SCALE
      LATT_CUR%A(2,:) =LATT_CUR%A(2,:)*SCALEY*LATT_CUR%SCALE
      LATT_CUR%A(3,:) =LATT_CUR%A(3,:)*SCALEZ*LATT_CUR%SCALE

      CALL LATTIC(LATT_CUR)

      IF (LATT_CUR%OMEGA<0) THEN
        WRITE(0,*)'ERROR: the triple product of the basis vectors ', &
     &     'is negative exchange two basis vectors'
        STOP
      ENDIF

! we are mainly interested in this (6th) line ...
      READ(IU,'(A)') INPLIN
! how many words/data items? --> number of ion types on file POSCAR!
      T_INFO%NTYP=NITEMS(INPLIN,INPWRK,.TRUE.,'I')
      ALLOCATE(T_INFO%NITYP(T_INFO%NTYP))
!-----number of atoms per type
      READ(INPLIN,*) (T_INFO%NITYP(NT),NT=1,T_INFO%NTYP)
! how many ions do we have on file POSCAR ... ?
      T_INFO%NIONS=0
      DO NI=1,T_INFO%NTYP
         T_INFO%NIONS=T_INFO%NIONS+T_INFO%NITYP(NI)
      END DO
      ALLOCATE(T_INFO%ITYP(T_INFO%NIONS))
!---- Set up the table from which we get type of each ion
      NI=1
      DO NT=1,T_INFO%NTYP
      DO NI=NI,T_INFO%NITYP(NT)+NI-1
        T_INFO%ITYP(NI)=NT
      ENDDO
      ENDDO
! posion
      READ(IU,'(A1)') CSEL
      IF (CSEL=='K'.OR.CSEL=='k'.OR. &
     &    CSEL=='C'.OR.CSEL=='c') THEN
        CSEL='K'
        WRITE(IO,*)'Positions in cartesian coordinates'
        T_INFO%LDIRCO=.FALSE.
      ELSE
        WRITE(IO,*)'Positions in direct lattice'
        T_INFO%LDIRCO=.TRUE.
      ENDIF
      ALLOCATE(T_INFO%POSION(3,T_INFO%NIONS)); T_INFO%POSION=0
      DO NI=1,T_INFO%NIONS
        READ(IU,*,ERR=400,END=400) T_INFO%POSION(:,NI)
      ENDDO

      IF (CSEL=='K') THEN
        T_INFO%POSION(1,:)=LATT_CUR%SCALE*T_INFO%POSION(1,:)*SCALEX
        T_INFO%POSION(2,:)=LATT_CUR%SCALE*T_INFO%POSION(2,:)*SCALEY
        T_INFO%POSION(3,:)=LATT_CUR%SCALE*T_INFO%POSION(3,:)*SCALEZ
        CALL KARDIR(T_INFO%NIONS,T_INFO%POSION,LATT_CUR%B)
      ENDIF
      CALL TOPRIM(T_INFO%NIONS,T_INFO%POSION)
      CLOSE(IU)
      RETURN

 400 CONTINUE
      WRITE(0,*)' No initial positions read in'
      STOP

      END SUBROUTINE READ_POSCAR

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE READ_KPOINTS(KPOINTS,FNAME,IU,IO)
      USE prec
      IMPLICIT NONE
      TYPE (kpoints_struct) KPOINTS
      INTEGER IO,IU
      CHARACTER(*) FNAME
      INTEGER NKP,INDX,NINTER,NKPX,NKPY,NKPZ,I,IX,IY,IZ,IRED
      INTEGER,PARAMETER :: NKDIM=100
      CHARACTER*1   CSEL,CLINE
      REAL(q) SHIFT(3),NKPIF(2,3),X(3)
! required for reallocation
      INTEGER IERR,N
      REAL(q),ALLOCATABLE   :: VKPT2(:,:)

      KPOINTS%NKPT=0
      WRITE(IO,'(" READING KPOINTS...")')
      WRITE(* ,'(" READING KPOINTS...")')
      OPEN(IU,FILE=FNAME,STATUS='OLD')
      READ(IU,'(A40)') KPOINTS%SZNAMK
      WRITE(IO,*) KPOINTS%SZNAMK
      READ(IU,*) NINTER
      READ(IU,'(A1)') CSEL
      IF (CSEL=='L'.OR.CSEL=='l') THEN
         CLINE='L'
         READ(IU,'(A1)') CSEL
         WRITE(IO,*)'Interpolating k-points between supplied coordinates'
      ELSE
         CLINE=" "
      ENDIF
      IF (CSEL=='K'.OR.CSEL=='k'.OR. &
     &    CSEL=='C'.OR.CSEL=='c') THEN
        CSEL='K'
        WRITE(IO,*)'k-points in cartesian coordinates'
      ELSE
        WRITE(IO,*)'k-points in reciprocal lattice'
      ENDIF
!=======================================================================
! read in a set of k-points and interpolate NKPTS between each
!=======================================================================
      IF (NINTER>0) THEN
      kr: IF (CLINE=='L') THEN
        ALLOCATE(VKPT2(3,NKDIM)); VKPT2=0
        NKP=0  ! counter for the number of k-points already read in
        DO
          NKP=NKP+1
          IF (NKP>NKDIM) THEN
            WRITE(0,*)'ERROR: MAIN: increase NKDIM'
            STOP
          ENDIF
          READ(IU,*,IOSTAT=IERR) VKPT2(:,NKP)
          IF (IERR/=0) EXIT
        ENDDO
        NKP=NKP-1
        KPOINTS%NKPTS=(NKP/2)*NINTER
        ALLOCATE(KPOINTS%VKPT(3,KPOINTS%NKPTS),KPOINTS%WTKPT(KPOINTS%NKPTS))
        IF (CSEL=='K') THEN
          VKPT2(:,1:NKP)=VKPT2(:,1:NKP)/LATT_CUR%SCALE
          CALL KARDIR(NKP,VKPT2,LATT_CUR%A)
        ENDIF

        INDX=0
        ! make NKPTS even
        NKP=(NKP/2)*2
        DO NKP=1,NKP-1,2
          SHIFT=(VKPT2(:,NKP+1)-VKPT2(:,NKP))/(NINTER-1)
          DO N=0,NINTER-1
          INDX=INDX+1
          KPOINTS%VKPT(:,INDX)=VKPT2(:,NKP)+SHIFT*N
          KPOINTS%WTKPT(INDX)=1._q/KPOINTS%NKPTS
          ENDDO
        ENDDO
      ELSE kr
        WRITE(0,'(" ERROR: INITIAL NKPTS>0 while not line mode! Not supported!")')
        STOP
      ENDIF kr
      ELSE
!=======================================================================
! Automatic generation of a mesh if NINTER<=0:
!=======================================================================
        WRITE(IO,*) 'Automatic generation of k-mesh.'
! k-lattice basis vectors in cartesian or reciprocal coordinates?
! Always use gamma center k-mesh.
        IF((CSEL.NE.'G').AND.(CSEL.NE.'g'))THEN
          WRITE(0 ,'(" CSEL \= G,g: Not supported! Switch to g-NO Inversion!")')
          WRITE(IO,'(" CSEL \= G,g: Not supported! Switch to g-NO Inversion!")')
          CSEL='g'
        ENDIF
! G: time reversal symmetry; g: no time reversal symmetry
! Here we give the Monkhorst-Pack conventions ... :
        READ(IU,*) NKPX,NKPY,NKPZ
        IF(FNAME.NE.'KPOINTS_FS')THEN
! For simplicity, adjust NKPX/Y/Z to odd
        NKPIF(:,1)=NKPX/2; NKPIF(:,2)=NKPY/2; NKPIF(:,3)=NKPZ/2
        NKPX=NKPX/2*2+1  ; NKPY=NKPY/2*2+1  ; NKPZ=NKPZ/2*2+1
        KPOINTS%NKPT(1)=NKPX; KPOINTS%NKPT(2)=NKPY; KPOINTS%NKPT(3)=NKPZ
        WRITE(IO,'(" Adjusted k-mesh:",3I3)')NKPX,NKPY,NKPZ
! Set up k-points
        IF(CSEL=='G')THEN
        IRED=1
        DO I=3,1,-1
          IF(NKPIF(1,I).EQ.0)CYCLE
          NKPIF(1,I)=0; IRED=I; EXIT
        ENDDO
        ENDIF
        ELSE ! KPOINTS_FS
        KPOINTS%NKPT(1)=NKPX; KPOINTS%NKPT(2)=NKPY; KPOINTS%NKPT(3)=NKPZ
        NKPIF(1,:)=0; NKPIF(2,:)=KPOINTS%NKPT(:)
        ENDIF

        KPOINTS%NKPTS=1
        DO I=1,3
          KPOINTS%NKPTS=KPOINTS%NKPTS*(NKPIF(2,I)+NKPIF(1,I)+1)
        ENDDO
        WRITE(IO,'(" Total number of k-points:",I4)')KPOINTS%NKPTS
        ALLOCATE(KPOINTS%VKPT(3,KPOINTS%NKPTS),KPOINTS%WTKPT(KPOINTS%NKPTS))
        NKP=0
! FOR XCRYSDEN
        DO IX=-NKPIF(1,1),NKPIF(2,1); X(1)=FLOAT(IX)
        DO IY=-NKPIF(1,2),NKPIF(2,2); X(2)=FLOAT(IY)
        DO IZ=-NKPIF(1,3),NKPIF(2,3); X(3)=FLOAT(IZ)
          NKP=NKP+1
          KPOINTS%VKPT(1,NKP)=X(1)/NKPX
          KPOINTS%VKPT(2,NKP)=X(2)/NKPY
          KPOINTS%VKPT(3,NKP)=X(3)/NKPZ
          IF(X(IRED).EQ.0._q.OR.CSEL=='g')THEN
            KPOINTS%WTKPT(NKP)=1._q/(NKPX*NKPY*NKPZ)
          ELSE
            KPOINTS%WTKPT(NKP)=2._q/(NKPX*NKPY*NKPZ)
          ENDIF
        ENDDO; ENDDO; ENDDO
      ENDIF
      CLOSE(IU)
      WRITE(IO,'(" KPOINTS (5,Direct):")')
      WRITE(IO,'(" NK   K    WT")')
      DO NKP=1,MIN(KPOINTS%NKPTS,5)
        WRITE(IO,'(I4,3F10.5,F10.5)')NKP,KPOINTS%VKPT(:,NKP),KPOINTS%WTKPT(NKP)
      ENDDO
      WRITE(IO,'(" KPOINTS (5,Cartesian):")')
      WRITE(IO,'(" NK   K    WT")')
      DO NKP=1,MIN(KPOINTS%NKPTS,5)
        WRITE(IO,'(I4,3F10.5,F10.5)')NKP,MATMUL(LATT_CUR%B,KPOINTS%VKPT(:,NKP)),KPOINTS%WTKPT(NKP)
      ENDDO
      RETURN

      END SUBROUTINE READ_KPOINTS

!***********************************************************************
      SUBROUTINE READ_INCAR(IU5,IU0)
      USE prec
      IMPLICIT NONE
      INTEGER IU5,IU0
! LOCAL
      INTEGER IDUM, N, IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LOPEN,LDUM
      CHARACTER*1 :: CHARAC

      LOPEN=.FALSE.
      OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')

      CL%LMULLIKEN = .TRUE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOLMULLIKEN','=','#',';','L', &
     &            IDUM,RDUM,CDUM,CL%LMULLIKEN,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOLMULLIKEN'' from file INCAR.'
         CL%LMULLIKEN = .FALSE.
      ENDIF

      CL%LROBY = .FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOLROBY','=','#',';','L', &
     &            IDUM,RDUM,CDUM,CL%LROBY,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOLROBY'' from file INCAR.'
         CL%LROBY = .FALSE.
      ENDIF

      CL%LSPDM = .TRUE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOLSPDM','=','#',';','L', &
     &            IDUM,RDUM,CDUM,CL%LSPDM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOLSPDM'' from file INCAR.'
         CL%LSPDM = .FALSE.
      ENDIF

      CL%LMAYER = .TRUE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'QOLMAYER','=','#',';','L', &
     &            IDUM,RDUM,CDUM,CL%LMAYER,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''QOLMAYER'' from file INCAR.'
         CL%LMAYER = .FALSE.
      ENDIF

      DOS%NGRID=301
      CALL RDATAB(LOPEN,'INCAR',IU5,'NEDOS','=','#',';','I', &
     &            DOS%NGRID,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NEDOS'' from file INCAR.'
         DOS%NGRID=301
      ENDIF

      DOS%ISMEAR=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'ISMEAR','=','#',';','I', &
     &            DOS%ISMEAR,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ISMEAR'' from file INCAR.'
         DOS%ISMEAR=1
      ENDIF

      DOS%SIGMA=0.2_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'SIGMA','=','#',';','F', &
     &            IDUM,DOS%SIGMA,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''SIGMA'' from file INCAR.'
         DOS%SIGMA=0.2_q
      ENDIF

      DOS%EMIN=10._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'EMIN','=','#',';','F', &
     &            IDUM,DOS%EMIN,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''EMIN'' from file INCAR.'
         DOS%EMIN=10._q
      ENDIF

      DOS%EMAX=-10._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'EMAX','=','#',';','F', &
     &            IDUM,DOS%EMAX,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''EMAX'' from file INCAR.'
         DOS%EMAX=-10._q
      ENDIF
      CLOSE(IU5)
      RETURN
   
      END SUBROUTINE READ_INCAR

!**************** SUBROUTINE LATTIC  ***********************************
!  subroutine for calculating the reciprocal lattice from the direct
!  lattice, in addition the norm of the lattice-vectors and the volume of
!  the basis-cell is calculated
!***********************************************************************
      SUBROUTINE LATTIC(Mylatt)
      USE prec
      IMPLICIT NONE

      TYPE(LATT) Mylatt
      REAL(q) Omega
      INTEGER I,J
      INTRINSIC SUM

      CALL EXPRO(Mylatt%B(1:3,1),Mylatt%A(1:3,2),Mylatt%A(1:3,3))
      CALL EXPRO(Mylatt%B(1:3,2),Mylatt%A(1:3,3),Mylatt%A(1:3,1))
      CALL EXPRO(Mylatt%B(1:3,3),Mylatt%A(1:3,1),Mylatt%A(1:3,2))

      Omega =Mylatt%B(1,1)*Mylatt%A(1,1)+Mylatt%B(2,1)*Mylatt%A(2,1) &
     &      +Mylatt%B(3,1)*Mylatt%A(3,1)

      DO I=1,3; DO J=1,3
        Mylatt%B(I,J)=Mylatt%B(I,J)/Omega
      ENDDO; ENDDO

      DO I=1,3
        Mylatt%ANORM(I)=SQRT(SUM(Mylatt%A(:,I)*Mylatt%A(:,I)))
        Mylatt%BNORM(I)=SQRT(SUM(Mylatt%B(:,I)*Mylatt%B(:,I)))
      ENDDO
      Mylatt%Omega=Omega
      RETURN
      END SUBROUTINE LATTIC

!***********************************************************************
      SUBROUTINE READ_TBHS(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      INTEGER IC,IQ1,IQ2,ISP
      INTEGER NQ,NC,ISPIN
      REAL(q) RES,RET

      WRITE(IO,'(" READ_TBHS...")')
      WRITE(* ,'(" READ_TBHS...")')
      OPEN(IU,FILE='TBHC.INP',STATUS='OLD')
! DOS INFO
      READ(IU,*)
      READ(IU,*)DOS%NET,DOS%EFERMI
      DOS%EFERMIN=DOS%EFERMI
! BV space CELL info
      READ(IU,*)
      READ(IU,*)NC
      ALLOCATE(CELL%R(3,NC))
      DO IC=1,NC; READ(IU,*)CELL%R(:,IC); ENDDO
! QO info
      READ(IU,*)
      READ(IU,*)NQ,ISPIN,W%RSPIN,QB%RCUT
      ALLOCATE(QB%NI(NQ),QB%LMCHAR(NQ))
      READ(IU,*)
      READ(IU,*)QB%NI
      READ(IU,*)QB%LMCHAR
! TBH
      READ(IU,*)
      ALLOCATE(QB%HR(NQ,NQ,NC,ISPIN)); QB%HR=0
      DO
#ifdef realmode
        READ(IU,*,ERR=102,END=102)ISP,IQ1,IQ2,IC,RES
#else
        READ(IU,*,ERR=102,END=102)ISP,IQ1,IQ2,IC,RES,RET
#endif
        IF(ABS(QB%HR(IQ1,IQ2,IC,ISP)).NE.0)THEN
        WRITE(0,'(" ISP,IQ1,IQ2,IC",3I3,I6," ENTRY REPEATED IN TBHC.INP!")')ISP,IQ1,IQ2,IC
        STOP
        ENDIF
#ifdef realmode
        QB%HR(IQ1,IQ2,IC,ISP)=RES
#else
        QB%HR(IQ1,IQ2,IC,ISP)=CMPLX(RES,RET)
#endif
      ENDDO
102   CLOSE(IU)
      OPEN(IU,FILE='TBSC.INP',STATUS='OLD')
      ALLOCATE(QB%SR(NQ,NQ,NC,ISPIN)); QB%SR=0
      READ(IU,*)
      DO 
        READ(IU,*,ERR=101,END=101)ISP,IQ1,IQ2,IC,RES
        IF(QB%SR(IQ1,IQ2,IC,ISP).NE.0)THEN
        WRITE(0,'(" ISP,IQ1,IQ2,IC",3I3,I6," ENTRY REPEATED IN TBHS.INP!")')ISP,IQ1,IQ2,IC
        STOP
        ENDIF
        QB%SR(IQ1,IQ2,IC,ISP)=RES
      ENDDO
101   CONTINUE
      CLOSE(IU)
! Locate origin
      DO IC=1,NC
        IF(CELL%R(1,IC).EQ.0.AND.CELL%R(2,IC).EQ.0.AND.CELL%R(3,IC).EQ.0)THEN
          CELL%ORIG=IC
          GOTO 100
        ENDIF
      ENDDO
      STOP 'NO ORIGIN IN CELL!'
100   CONTINUE
#ifndef realmode
      WRITE(IO,'(" MAXIMAL IMAGINARY H:",F12.6)')MAXVAL(ABS(AIMAG(QB%HR)))
#endif
      QB%DIM=NQ; CELL%DIM=NC; QB%ISPIN=ISPIN
      W%ISPIN=ISPIN
      RETURN

      END SUBROUTINE READ_TBHS

!***********************************************************************
      SUBROUTINE SET_BAND_DOS(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO

      W%NBANDS=QB%DIM
! WAVE
      ALLOCATE(W%CPTWFP(QB%DIM,W%NBANDS,KPOINTS%NKPTS,W%ISPIN), &
              &W%CELEN (       W%NBANDS,KPOINTS%NKPTS,W%ISPIN), &
              &W%FERWE (       W%NBANDS,KPOINTS%NKPTS,W%ISPIN), &
              &W%PAR   (QB%DIM,W%NBANDS,KPOINTS%NKPTS,W%ISPIN))
      W%CPTWFP=0; W%CELEN=0; W%FERWE=0; W%PAR=0
! DOS
      IF(DOS%ISMEAR.LT.0)THEN
        WRITE(IO,'(" NOT SUPPORT ISMEAR<0! SWITCH TO GAUSSIAN SMEARING!")')
        DOS%ISMEAR=0; DOS%SIGMA=0.01_q
      ENDIF
      ALLOCATE(DOS%T(DOS%NGRID,W%ISPIN),DOS%I(DOS%NGRID,W%ISPIN),DOS%PAR(DOS%NGRID,QB%DIM,W%ISPIN))
      DOS%T=0; DOS%I=0; DOS%PAR=0
      RETURN

      END SUBROUTINE SET_BAND_DOS

!***********************************************************************
      SUBROUTINE SET_ATOM_SHELL(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER,PARAMETER::NSL=20,NIMAX=100
      INTEGER NI1,NI2,IC,NDIM,INN,I
      REAL(q) DPOS(3),DIST
      REAL(q),ALLOCATABLE::PAIRDIST(:,:,:),R(:),RNN(:)
      INTEGER,ALLOCATABLE::IDX(:),NNN(:),ICNN(:,:),NINN(:,:)
   
      WRITE(IO,'(" SET_ATOM_SHELL...")')
      WRITE(* ,'(" SET_ATOM_SHELL...")')
      ALLOCATE(PAIRDIST(T_INFO%NIONS,T_INFO%NIONS,CELL%DIM)); PAIRDIST=0
      DO NI1=1,T_INFO%NIONS; DO NI2=1,T_INFO%NIONS; DO IC=1,CELL%DIM
      DPOS=CELL%R(:,IC)+T_INFO%POSION(:,NI2)-T_INFO%POSION(:,NI1)
      DPOS=MATMUL(LATT_CUR%A,DPOS)
      PAIRDIST(NI1,NI2,IC)=SQRT(SUM(DPOS**2))
      ENDDO; ENDDO; ENDDO ! IC,NI2,NI1

      ALLOCATE(YZ(T_INFO%NIONS))
      ALLOCATE(RNN(NSL),NNN(NSL),ICNN(NIMAX,NSL),NINN(NIMAX,NSL))
      NDIM=T_INFO%NIONS*CELL%DIM
      ALLOCATE(IDX(NDIM),R(NDIM))
      DO NI1=1,T_INFO%NIONS
        R=RESHAPE(PAIRDIST(NI1,:,:),(/NDIM/))
        CALL ORDER(NDIM,R,IDX)
        DIST=0; INN=0; NNN=0; RNN=0; ICNN=0
        DO I=2,NDIM
          IF(R(IDX(I)).GT.QB%RCUT)EXIT
          IF(ABS(R(IDX(I))-DIST).GT..05_q)THEN
            IF(INN.GE.NSL)THEN
              WRITE(IO,'(" WARNING: REACH MAX N.SHELL!")'); EXIT
            ENDIF
            INN=INN+1
            RNN(INN)=R(IDX(I))
            DIST=R(IDX(I))
          ENDIF
          IF(NNN(INN).GE.NIMAX)THEN
            WRITE(0,'(" NI=",I3," SHELL",I3," HAS ATOMS >",I3)')NI1,INN,NIMAX; STOP
          ENDIF
          NNN(INN)=NNN(INN)+1
          NI2=MOD(IDX(I)-1,T_INFO%NIONS)+1; IC=(IDX(I)-1)/T_INFO%NIONS+1
          NINN(NNN(INN),INN)=NI2; ICNN(NNN(INN),INN)=IC
        ENDDO
        YZ(NI1)%NSL=INN; ALLOCATE(YZ(NI1)%SL(0:INN))
        YZ(NI1)%SL(0)%N=1
        YZ(NI1)%SL(0)%DIST=0._q
        ALLOCATE(YZ(NI1)%SL(0)%NI(1),YZ(NI1)%SL(0)%IC(1))
        YZ(NI1)%SL(0)%NI=NI1
        YZ(NI1)%SL(0)%IC=CELL%ORIG
        DO I=1,INN
          YZ(NI1)%SL(I)%N=NNN(I)
          YZ(NI1)%SL(I)%DIST=RNN(I)
          ALLOCATE(YZ(NI1)%SL(I)%NI(NNN(I)),YZ(NI1)%SL(I)%IC(NNN(I)))
          YZ(NI1)%SL(I)%NI=NINN(1:NNN(I),I)
          YZ(NI1)%SL(I)%IC=ICNN(1:NNN(I),I)
        ENDDO
        WRITE(IO,'(" NI=",I3," HAS",I3," SHELLS:")')NI1,INN
        DO I=1,INN
          WRITE(IO,'("  SHELL ",I3," DIST=",F6.2," HAS",I3," ATOMS.")')I,RNN(I),NNN(I)
          DO NI2 =1,NNN(I)
          WRITE(IO,'("   ATOM",I3," CELL%R=",3I3)')NINN(NI2,I),CELL%R(:,ICNN(NI2,I))
          ENDDO
        ENDDO
      ENDDO
      WRITE(IO,'(" MAX. # SHELLS:",I3)')MAXVAL(YZ(:)%NSL)
      CALL CALC_MAX_RCUT(CELL%DIM,CELL%R,LATT_CUR%A,QB%MAXRCUT)
      WRITE(IO,'(" QB%RCUT=",F10.5," QB%MAXRCUT=",F10.5)')QB%RCUT,QB%MAXRCUT
      RETURN

      END SUBROUTINE SET_ATOM_SHELL

!***********************************************************************
      SUBROUTINE OUT_TBHS(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER IMODE,ISP

      WRITE(*,'(" OUT_TBHS...")')
      DO ISP=1,W%ISPIN; DO IMODE=1,2
        CALL OUT_TBHSS(IO,IMODE,ISP)
      ENDDO; ENDDO
      RETURN

      END SUBROUTINE OUT_TBHS

!***********************************************************************
      SUBROUTINE OUT_TBHSS(IO,IMODE,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER IO,IMODE,ISP
! LOCAL
      INTEGER NI1,NI2,I,J,IA,IC
      
      IF(IMODE.EQ.1)THEN
        WRITE(IO,'(/," PAIR HOPPING DATA WITH ISP=",I1)')ISP
      ELSE
        WRITE(IO,'(/," PAIR OVERLAPPING DATA WITH ISP=",I1)')ISP
      ENDIF
      DO NI1=1,T_INFO%NIONS
      DO I=1,MIN(YZ(NI1)%NSL,9) ! Upto 10th NN shell
      WRITE(IO,'(" ATOM",I3," NN-",I1," DIST=",F6.2)')NI1,I,YZ(NI1)%SL(I)%DIST
      DO J=1,YZ(NI1)%SL(I)%N; NI2=YZ(NI1)%SL(I)%NI(J); IC=YZ(NI1)%SL(I)%IC(J)
      WRITE(IO,'(" ATOM",I3,"|--ATOM",I3," UNIT CELL=",3I3," ISP=",I1)')NI1,NI2,CELL%R(:,IC),ISP
      WRITE(IO,101)QB%LMCHAR(YZ(NI2)%IQMIN:YZ(NI2)%IQMAX)
      DO IA= YZ(NI1)%IQMIN,YZ(NI1)%IQMAX
      IF(IMODE.EQ.1)THEN
        WRITE(IO,100)QB%LMCHAR(IA),QB%HR(IA,YZ(NI2)%IQMIN:YZ(NI2)%IQMAX,IC,ISP)
      ELSE
        WRITE(IO,100)QB%LMCHAR(IA),QB%SR(IA,YZ(NI2)%IQMIN:YZ(NI2)%IQMAX,IC,ISP)
      ENDIF
      ENDDO ! IA
      ENDDO; ENDDO ! I,J
      ENDDO ! NI1
100   FORMAT(A3,1X,<YZ(NI2)%NQ*2>F15.10)
101   FORMAT(4X,<YZ(NI2)%NQ>(A7))
      RETURN

      END SUBROUTINE OUT_TBHSS

!***********************************************************************
      SUBROUTINE SET_ATOM_QO(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER IQ,NI
      
      YZ(:)%IQMIN=0; YZ(:)%IQMAX=0
      DO IQ=1,QB%DIM; NI=QB%NI(IQ)
        IF(YZ(NI)%IQMIN.EQ.0)THEN
          YZ(NI)%IQMIN=IQ; YZ(NI)%IQMAX=IQ
        ELSE
          YZ(NI)%IQMAX=IQ
        ENDIF
      ENDDO
      YZ(:)%NQ=YZ(:)%IQMAX-YZ(:)%IQMIN+1
      RETURN

      END SUBROUTINE SET_ATOM_QO

!***********************************************************************
      SUBROUTINE SET_HS1K(VKPT,AK,AR)
      USE prec; USE constant
      IMPLICIT NONE
      REAL(q) VKPT(3)
      OVERLAP AR(QB%DIM,QB%DIM,CELL%DIM)
      COMPLEX(q) AK(QB%DIM,QB%DIM)
! LOCAL
      INTEGER IC,IQ1,IQ2
      COMPLEX(q) PHASEK

      AK=0
      DO IC=1,CELL%DIM
      PHASEK=EXP(-CITPI*SUM(CELL%R(:,IC)*VKPT) ) ! EXP(-iK.R)
      DO IQ1=1,QB%DIM; DO IQ2=IQ1,QB%DIM
        AK(IQ1,IQ2)=AK(IQ1,IQ2)+AR(IQ1,IQ2,IC)*PHASEK
      ENDDO; ENDDO; ENDDO
      DO IQ1=1,QB%DIM; DO IQ2=1,IQ1-1
        AK(IQ1,IQ2)=CONJG(AK(IQ2,IQ1))
      ENDDO; ENDDO
      RETURN

      END SUBROUTINE SET_HS1K

!***********************************************************************
      SUBROUTINE CALC_FULL_BAND(LORTH,IO)
      USE prec
      IMPLICIT NONE
      LOGICAL LORTH
      INTEGER IO
! LOCAL
      INTEGER ISP,K
      REAL(q) KPT(3)

      WRITE(IO,'(" CALC_FULL_BAND...")')
      WRITE(* ,'(" CALC_FULL_BAND...")')
      DO ISP=1,W%ISPIN
      DO K=1,KPOINTS%NKPTS
        KPT=KPOINTS%VKPT(:,K)
        CALL EIGEN_SOL(KPT,K,W%CELEN(:,K,ISP),W%CPTWFP(:,:,K,ISP),QB%DIM,'V',ISP,LORTH,W%PAR(:,:,K,ISP))
      ENDDO; ENDDO ! ISP,K
      RETURN

      END SUBROUTINE CALC_FULL_BAND

!***********************************************************************
      SUBROUTINE EIGEN_SOL(KPT,K,W,A,N,JOBZ,ISP,LORTH,PJ)
      USE prec
      IMPLICIT NONE
      INTEGER K,N,ISP
      REAL(q) W(N),KPT(3)
      COMPLEX(q) A(N,N)
      CHARACTER*1 JOBZ
      LOGICAL LORTH
      REAL(q),OPTIONAL,INTENT(OUT)::PJ(N,N)
! Local
      INTEGER LWORK,INFO,IV,IVP
      REAL(q) RWORK(3*N-2)
      COMPLEX(q),ALLOCATABLE :: B(:,:),WORK(:),BP(:,:),SV(:)

      CALL  GET_EFF_HAMILT(A,N,KPT,K,ISP,LORTH)
      LWORK=32*N
      ALLOCATE(WORK(LWORK)); WORK=0
      IF(.NOT.LORTH)THEN
      ALLOCATE(B(N,N)); B=0
      CALL  GET_EFF_OVERLAP(B,N,KPT,K,ISP)
      IF(PRESENT(PJ))THEN
      ALLOCATE(BP(N,N),SV(N)); BP=B
      ENDIF
      CALL ZHEGV(1,JOBZ,'L',N,A,N,B,N,W,WORK,LWORK,RWORK,INFO)
      IF(PRESENT(PJ))THEN
      DO IV=1,N
        SV=MATMUL(BP,A(:,IV))
        DO IVP=1,N
        PJ(IVP,IV)=CONJG(A(IVP,IV))*SV(IVP)
      ENDDO; ENDDO
      ENDIF
      ELSE
      CALL ZHEEV(    JOBZ,'L',N,A,N,    W,WORK,LWORK,RWORK,INFO)
      IF(PRESENT(PJ))THEN
        DO IV=1,N; DO IVP=1,N
        PJ(IVP,IV)=CONJG(A(IVP,IV))*A(IVP,IV)
        ENDDO; ENDDO
      ENDIF
      ENDIF
      DEALLOCATE(WORK)
      IF(.NOT.LORTH)THEN
      DEALLOCATE(B)
      IF(PRESENT(PJ))DEALLOCATE(BP)
      ENDIF
      IF(INFO.NE.0)THEN
        WRITE(0,'(" EIGEN_SOL: INFO=",I4)')INFO
        STOP
      ENDIF
      RETURN

      END SUBROUTINE EIGEN_SOL

!***********************************************************************
      SUBROUTINE GET_EFF_HAMILT(HMAT,NDIM,KPT,K,ISP,LORTH)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(q) KPT(3)
      COMPLEX(q) HMAT(NDIM,NDIM)
      LOGICAL LORTH

      CALL GET_UCTBH(HMAT,NDIM,KPT,K,ISP,LORTH)
      RETURN

      END SUBROUTINE GET_EFF_HAMILT

!***********************************************************************
      SUBROUTINE GET_UCTBH(HMAT,NDIM,KPT,K,ISP,LORTH)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(q) KPT(3)
      COMPLEX(q) HMAT(NDIM,NDIM)
      LOGICAL LORTH
! LOCAL
      INTEGER IQ1,IQ2
      COMPLEX(q),ALLOCATABLE :: SK(:,:)
      OVERLAP,ALLOCATABLE::CR(:,:,:)

      CALL SET_HS1K(KPT,HMAT,QB%HR(:,:,:,ISP))
      IF(LORTH)THEN
      ALLOCATE(SK(QB%DIM,QB%DIM)); SK=0
      ALLOCATE(CR(QB%DIM,QB%DIM,CELL%DIM)); CR=QB%SR(:,:,:,ISP)
      CALL SET_HS1K(KPT,SK,CR)
      CALL ORTH_HS1K(HMAT,SK)
      DEALLOCATE(SK,CR)
      ENDIF
      RETURN

      END SUBROUTINE GET_UCTBH

!***********************************************************************
      SUBROUTINE GET_EFF_OVERLAP(SMAT,NDIM,KPT,K,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(q) KPT(3)
      COMPLEX(q) SMAT(NDIM,NDIM)

      CALL GET_UCTBS(SMAT,NDIM,KPT,K,ISP)
      RETURN

      END SUBROUTINE GET_EFF_OVERLAP

!***********************************************************************
      SUBROUTINE GET_UCTBS(SMAT,NDIM,KPT,K,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(q) KPT(3)
      COMPLEX(q) SMAT(NDIM,NDIM)
      OVERLAP,ALLOCATABLE::CR(:,:,:)

      ALLOCATE(CR(QB%DIM,QB%DIM,CELL%DIM)); CR=QB%SR(:,:,:,ISP)
      CALL SET_HS1K(KPT,SMAT,CR)
      DEALLOCATE(CR)
      RETURN

      END SUBROUTINE GET_UCTBS

!***********************************************************************
      SUBROUTINE SET_FERMI_WT(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER ISP,K,N
      REAL(q) EADD,E

      IF(DOS%EMIN.GT.DOS%EMAX)THEN
        DO ISP=1,W%ISPIN
        DO K=1,KPOINTS%NKPTS
        DO N=1,W%NBANDS
          E=W%CELEN(N,K,ISP)
          DOS%EMAX=MAX(DOS%EMAX,E)
          DOS%EMIN=MIN(DOS%EMIN,E)
        ENDDO; ENDDO; ENDDO
        EADD=(DOS%EMAX-DOS%EMIN)*0.05_q
        EADD=MAX(EADD,10*ABS(DOS%SIGMA))
        DOS%EMIN=DOS%EMIN-EADD
        DOS%EMAX=DOS%EMAX+EADD
      ENDIF
! Determine fermi level
      CALL DENMP(W%ISPIN,W%RSPIN,DOS%EMIN,DOS%EMAX,DOS%NET,DOS%EFERMI, &
             & DOS%ISMEAR,DOS%SIGMA,W%NBANDS,KPOINTS%NKPTS,W%FERWE,W%CELEN, &
             & DOS%NGRID,DOS%T,DOS%I,W%PAR,DOS%PAR,KPOINTS%WTKPT)
      WRITE(IO,'(" EFERMI=",F10.5," EFERMI_INPUT=",F10.5)')DOS%EFERMI,DOS%EFERMIN
      RETURN

      END SUBROUTINE SET_FERMI_WT

!********************************************************************
! PRINT OUT TDOS
!********************************************************************
      SUBROUTINE OUT_DOS(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      INTEGER I,NI,ISP
      REAL(q) DELTAE

      WRITE(IO,'(" OUT_DOS...")')
      WRITE(* ,'(" OUT_DOS...")')
      DELTAE=(DOS%EMAX-DOS%EMIN)/(DOS%NGRID-1)
      OPEN(IU,FILE='TB_DOS',STATUS='REPLACE')
      WRITE(IU,'(" TB_TDOS")')
      WRITE(IU,'(F10.4,F9.3)')DOS%EFERMI,0_q
      WRITE(IU,'(F10.4,F9.3)')DOS%EFERMI,MAXVAL(DOS%T)
      DO I=1,DOS%NGRID
        WRITE(IU,'(F10.4,4F9.3)')DOS%EMIN+(I-1)*DELTAE,DOS%T(I,:),DOS%I(I,:)
      ENDDO
      DO NI=1,T_INFO%NIONS
      WRITE(IU,'(" TB_PAR, NI=",I4)')NI
      WRITE(IU,'(F10.4,F9.3)')DOS%EFERMI,0_q
      WRITE(IU,'(F10.4,F9.3)')DOS%EFERMI,MAXVAL(DOS%PAR(:,YZ(NI)%IQMIN:YZ(NI)%IQMAX,:))
      WRITE(IU,'(10X,32(4X,A3))')(QB%LMCHAR(YZ(NI)%IQMIN:YZ(NI)%IQMAX),ISP=1,W%ISPIN)
      DO I=1,DOS%NGRID
        WRITE(IU,'(F10.4,32F7.3)')DOS%EMIN+(I-1)*DELTAE,DOS%PAR(I,YZ(NI)%IQMIN:YZ(NI)%IQMAX,:)
      ENDDO
      ENDDO
      CLOSE(IU)
      RETURN

      END SUBROUTINE OUT_DOS

!********************************************************************
      SUBROUTINE ORTH_HS1K(HK,SK)
      USE prec
      IMPLICIT NONE
      COMPLEX(q) HK(QB%DIM,QB%DIM),SK(QB%DIM,QB%DIM)
! LOCAL
      INTEGER LWORK,IERR,IQ,IQ1,IQ2
      COMPLEX(q) SNHK(QB%DIM,QB%DIM),SK_(QB%DIM,QB%DIM)
      REAL(q) RWORK(3*QB%DIM),W(QB%DIM)
      COMPLEX(q) WORK(32*QB%DIM)

      LWORK=32*QB%DIM; SK_=SK
      CALL ZHEEV('V','U',QB%DIM,SK,QB%DIM,W,WORK,LWORK,RWORK,IERR)
      IF(IERR.NE.0)STOP 'ZHEEV ERROR-I IN RTH_HS1K!'
      IF(MINVAL(W)/MAXVAL(W).LT.1.E-5_q)THEN ! Condition number
      WRITE(0,*)
      WRITE(0,*)'EIGENVALUES IN ORTH_HS1K:'
      WRITE(0,'(10F10.5)')W
      STOP 'ZHEEV ERROR-II IN RTH_HS1K'
      ENDIF
! Set up S^(-1/2)
      SNHK=0
      DO IQ1=1,QB%DIM; DO IQ2=1,QB%DIM; DO IQ=1,QB%DIM
        SNHK(IQ2,IQ1)=SNHK(IQ2,IQ1)+SK(IQ2,IQ)*CONJG(SK(IQ1,IQ))/SQRT(W(IQ))
      ENDDO; ENDDO; ENDDO
      HK=MATMUL(CONJG(TRANSPOSE(SNHK)),MATMUL(HK ,SNHK))
      SK=MATMUL(CONJG(TRANSPOSE(SNHK)),MATMUL(SK_,SNHK))
      RETURN

      END SUBROUTINE ORTH_HS1K

!***********************************************************************
      SUBROUTINE CALC_BAND_STRU(LORTH,IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
      LOGICAL LORTH,LEXIST

      WRITE(IO,'(" CALC_BAND_STRU...")')
      WRITE(* ,'(" CALC_BAND_STRU...")')
      INQUIRE(FILE='KPOINTS_BAND',EXIST=LEXIST)
      IF(.NOT.LEXIST)THEN
        WRITE(IO,'(" FILE KPOINTS_BAND NOT EXIST, SKIPPED.")')
        WRITE(* ,'(" FILE KPOINTS_BAND NOT EXIST, SKIPPED.")')
        RETURN
      ENDIF
      CALL READ_KPOINTS(KPT_BAND,'KPOINTS_BAND',IU,IO)
      CALL CALC_BAND_STRU_(LORTH,IU)
      RETURN

      END SUBROUTINE CALC_BAND_STRU

!***********************************************************************
      SUBROUTINE CALC_BAND_STRU_(LORTH,IU)
      USE prec
      IMPLICIT NONE
      INTEGER IU
      LOGICAL LORTH
! Local
      INTEGER ISP,NK
      REAL(q) EVAL(QB%DIM),DKLEN,KLEN
      COMPLEX(q) EVEC(QB%DIM,QB%DIM)

      DO ISP=1,W%ISPIN
      IF(ISP.EQ.1)THEN
        OPEN(IU,FILE='TBBAND.1',STATUS='REPLACE')
      ELSE
        OPEN(IU,FILE='TBBAND.2',STATUS='REPLACE')
      ENDIF
      KLEN=0
      DO NK=1,KPT_BAND%NKPTS
        CALL EIGEN_SOL(KPT_BAND%VKPT(:,NK),-NK,EVAL,EVEC,QB%DIM,'N',ISP,LORTH)
        IF(NK.EQ.1)THEN
          DKLEN=0
        ELSE
          CALL CALC_DLEN(DKLEN,KPT_BAND%VKPT(:,NK-1),KPT_BAND%VKPT(:,NK),LATT_CUR%B,LATT_CUR%SCALE)
        ENDIF
        KLEN=KLEN+DKLEN
        WRITE(IU,'(100F10.4)')KLEN,EVAL-DOS%EFERMI
      ENDDO
      CLOSE(IU)
      ENDDO ! ISPIN
      RETURN

      END SUBROUTINE CALC_BAND_STRU_

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

!********************************************************************
      SUBROUTINE SPDM_MULLIKEN_MAYER_ANALYSIS(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      REAL(q),POINTER :: PIAJB(:,:,:,:) ! Single particle density matrix 

      IF(.NOT.CL%LMULLIKEN.AND..NOT.CL%LMAYER.AND..NOT.CL%LSPDM)RETURN
      WRITE(IO,'(" SPDM_MULLIKEN_MAYER_ANALYSIS...")')
      WRITE(* ,'(" SPDM_MULLIKEN_MAYER_ANALYSIS...")')
      ALLOCATE(PIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)); PIAJB=0
      CALL CALC_FULL_PIAJB(PIAJB)
      CALL OUT_SPDM(PIAJB,IO)
      CALL CALC_MULLIKEN(PIAJB,IO)
      CALL CALC_MAYER_BO(PIAJB,IO)
      DEALLOCATE(PIAJB)
      RETURN

      END SUBROUTINE SPDM_MULLIKEN_MAYER_ANALYSIS

!********************************************************************
      SUBROUTINE ROBY_COV_ION_BO_ANALYSIS(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER NI,NIP,NIP_,ISL,NBASE,IQ,NIP__
      INTEGER NG(2),NQG(2),NG12,NQG12
      INTEGER,ALLOCATABLE::NIG(:),CRG(:,:),IQG(:),IQNIG(:)

      IF(.NOT.CL%LROBY)RETURN
      WRITE(IO,*)
      WRITE(IO,'(" ROBY_COV_ION_BO_ANALYSIS...")')
      WRITE(* ,'(" ROBY_COV_ION_BO_ANALYSIS...")')
      DO NI=1,T_INFO%NIONS
      NG(1)=1; NG(2)=SUM(YZ(NI)%SL(:)%N); NG12=NG(1)+NG(2)
      ALLOCATE(NIG(NG12),CRG(3,NG12)); CRG=0
      NIG(1)=NI; NQG(1)=YZ(NI)%NQ; NBASE=1; NQG(2)=0
      DO ISL=1,YZ(NI)%NSL; DO NIP_=1,YZ(NI)%SL(ISL)%N
      NIP=YZ(NI)%SL(ISL)%NI(NIP_)
      NBASE=NBASE+1
      NIG(NBASE)=NIP
      CRG(:,NBASE)=CELL%R(:,YZ(NI)%SL(ISL)%IC(NIP_))
      NQG(2)=NQG(2)+YZ(NIP)%NQ
      ENDDO; ENDDO
      NQG12=NQG(1)+NQG(2)
      ALLOCATE(IQG(NQG12),IQNIG(NQG12))
      DO IQ=1,YZ(NI)%NQ; IQG(IQ)=YZ(NI)%IQMIN-1+IQ; IQNIG(IQ)=1; ENDDO
      NBASE=YZ(NI)%NQ; NIP__=1
      DO ISL=1,YZ(NI)%NSL; DO NIP_=1,YZ(NI)%SL(ISL)%N
      NIP=YZ(NI)%SL(ISL)%NI(NIP_); NIP__=NIP__+1
      DO IQ=1,YZ(NIP)%NQ
      NBASE=NBASE+1
      IQG(NBASE)=YZ(NIP)%IQMIN-1+IQ; IQNIG(NBASE)=NIP__
      ENDDO; ENDDO; ENDDO
      WRITE(IO,'(" FOR NI=",I3," WITH NNs:")')NI
      CALL CALC_ROBY_INDEX(NG,NQG,NG12,NQG12,CRG,IQG,IQNIG,IO)
      DEALLOCATE(NIG,CRG,IQG,IQNIG)
      ENDDO ! NI
      RETURN

      END SUBROUTINE ROBY_COV_ION_BO_ANALYSIS

!********************************************************************
      SUBROUTINE CALC_ROBY_INDEX(NG,NQG,NG12,NQG12,CRG,IQG,IQNIG,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
      INTEGER NG(2),NQG(2),NG12,NQG12
      INTEGER CRG(3,NG12),IQG(NQG12),IQNIG(NQG12)
! LOCAL
      INTEGER ISP,IQ1,IQ2,I,J,IP,JP,J_,JP_
      REAL(q) CAB(W%ISPIN),IAB(W%ISPIN),TAB(W%ISPIN)
      REAL(q),ALLOCATABLE::S(:,:),SAINV(:,:),SBINV(:,:),A(:,:),SINV(:,:),PA(:,:),PB(:,:),PAB(:,:),B(:,:)
      REAL(q),ALLOCATABLE::PIAJB(:,:)

      ALLOCATE(S(NQG12,NQG12),A(NQG12,NQG12),B(NQG12,NQG12),SINV(NQG12,NQG12),SAINV(NQG(1),NQG(1)), &
              &SBINV(NQG(2),NQG(2)),PA(NQG12,NQG12),PB(NQG12,NQG12),PAB(NQG12,NQG12),PIAJB(NQG12,NQG12))
      DO ISP=1,W%ISPIN
      WRITE(IO,'(" ISP=",I2)')ISP
      CALL MATCH_TBS_SIJ(S,NG12,NQG,IQG,NQG12,CRG,IQNIG,ISP)
      CALL CALC_SELECTIVE_PIAJB(PIAJB,NQG12,IQG,IQNIG,CRG,NQG,NG12,ISP)
      SINV =S                               ; CALL DGETRI_(SINV ,NQG12 )
      SAINV=S(1:NQG(1)      ,1:NQG(1)      ); CALL DGETRI_(SAINV,NQG(1))
      SBINV=S(1+NQG(1):NQG12,1+NQG(1):NQG12); CALL DGETRI_(SBINV,NQG(2))
! SET-UP MATRICES: PA,PB,PAB
!      PA=0; PB=0; PAB=0
!      DO IQ1=1,NQG12; DO IQ2=1,NQG12
!      J_ =0
!      DO I =1,2; DO J =1,NQG(I ); J_ =J_ +1
!      JP_=0
!      DO IP=1,2; DO JP=1,NQG(IP); JP_=JP_+1
!      IF(I.EQ.1.AND.IP.EQ.1)PA(IQ1,IQ2)=PA(IQ1,IQ2)+S(IQ1,J_)*SAINV(J,JP)*S(JP_,IQ2)
!      IF(I.EQ.2.AND.IP.EQ.2)PB(IQ1,IQ2)=PB(IQ1,IQ2)+S(IQ1,J_)*SBINV(J,JP)*S(JP_,IQ2)
!      PAB(IQ1,IQ2)=PAB(IQ1,IQ2)+S(IQ1,J_)*SINV(J_,JP_)*S(JP_,IQ2)
!      ENDDO; ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
      PA=0; A=0
      CALL DGEMM('N','N',NQG(1),NQG12,NQG(1),1._q,SAINV              ,NQG(1),S(1:NQG(1),1:NQG12),NQG(1), &
                &0._q,A(1:NQG(1),1:NQG12),NQG(1))
      CALL DGEMM('N','N',NQG12 ,NQG12,NQG(1),1._q,S(1:NQG12,1:NQG(1)),NQG12,A(1:NQG(1),1:NQG12),NQG(1), &
                &0._q,PA                 ,NQG12)
      PB=0; A=0
      CALL DGEMM('N','N',NQG(2),NQG12,NQG(2),1._q,SBINV                    ,NQG(2),S(NQG(1)+1:NQG12,1:NQG12),NQG(2), &
                &0._q,A(1:NQG(2),1:NQG12),NQG(2))
      CALL DGEMM('N','N',NQG12 ,NQG12,NQG(2),1._q,S(1:NQG12,1+NQG(1):NQG12),NQG12,A(1:NQG(2),1:NQG12)      ,NQG(2), &
                &0._q,PB                 ,NQG12)
      PAB=0; A=0
      CALL DGEMM('N','N',NQG12,NQG12,NQG12,1._q,SINV,NQG12,S,NQG12,0._q,A  ,NQG12)
      CALL DGEMM('N','N',NQG12,NQG12,NQG12,1._q,S   ,NQG12,A,NQG12,0._q,PAB,NQG12)

! Covalent operator: PA+PB-PAB
      A=PA+PB-PAB; B=S
      CALL CALC_BOND_ORDER_ROBY_(A,B,CAB(ISP),PIAJB,NQG12,NQG,IO)
      WRITE(IO,'(" CAB=",F10.5)')CAB(ISP)
! Ionic operator
      A=PA-PB; B=S
      CALL CALC_BOND_ORDER_ROBY_(A,B,IAB(ISP),PIAJB,NQG12,NQG,IO)
      WRITE(IO,'(" IAB=",F10.5)')IAB(ISP)
      TAB(ISP)=SQRT(CAB(ISP)**2+IAB(ISP)**2)
      WRITE(IO,'(" BOND INDEX TAB=",F10.5)')TAB(ISP)
      WRITE(IO,'(" COVALENCY=",F6.1,"%, IONICITY=",F6.1,"%.")')(CAB(ISP)/TAB(ISP))**2*100,(IAB(ISP)/TAB(ISP))**2*100
      ENDDO ! ISP
      CAB(1)=SUM(CAB)/W%ISPIN; IAB(1)=SUM(IAB)/W%ISPIN
      TAB(1)=SQRT(CAB(1)**2+IAB(1)**2)
      WRITE(IO,'(" SPIN-AVERAGED:")')
      WRITE(IO,'(" BOND INDEX TAB=",F10.5)')TAB(1)
      WRITE(IO,'(" COVALENCY=",F6.1,"%, IONICITY=",F6.1,"%.")')(CAB(1)/TAB(1))**2*100,(IAB(1)/TAB(1))**2*100
      DEALLOCATE(S,A,B,SINV,SAINV,SBINV,PA,PB,PAB)

      RETURN

      END SUBROUTINE CALC_ROBY_INDEX

!=======================================================================
      SUBROUTINE CALC_BOND_ORDER_ROBY_(A,B,CIAB,PIAJB,N,NG,IO)
      IMPLICIT NONE
      INTEGER N,NG(2),IO
      REAL(q) CIAB,PIAJB(N,N),A(N,N),B(N,N)
! LOCAL
      INTEGER INFO,LDWORK,IPAIR(2,N/2),NPAIR,I,NN,IP,J
      REAL(q) FPAIR(2,2,N/2),FRAC1,FRAC2,OCC(2),EPAIR(N/2)
      REAL(q),ALLOCATABLE::WORK(:),W(:)

      LDWORK=8*N; CIAB=0
      ALLOCATE(W(N),WORK(LDWORK)); W=0; WORK=0
      CALL DSYGV(1,'V','U',N,A,N,B,N,W,WORK,LDWORK,INFO)
      IF(INFO.NE.0)THEN
        WRITE(0,'(" INFO=",I5)')INFO
        STOP'FETAL ERROR: DSPGV IN CALC_BOND_ORDER_ROBY_!'
      ENDIF
      NPAIR=0; FPAIR=0
      DO I=1,N; IF(W(I).GT.0)EXIT; ENDDO; NN=I-1
      DO I=1,NN
        FRAC1=SUM(ABS(A(1:NG(1)  ,I))); IF(FRAC1.LT..001_q)CYCLE
        FRAC2=SUM(ABS(A(1+NG(1):N,I))); IF(FRAC2.LT..001_q)CYCLE
        IF(ABS(W(I)).GT..999_q)CYCLE
        NPAIR=NPAIR+1
        FPAIR(1,1,NPAIR)=FRAC1; FPAIR(2,1,NPAIR)=FRAC2
        EPAIR(NPAIR)=-W(I)
        IPAIR(1,NPAIR)=I
        DO IP=NN+1,N; IF(ABS(W(IP)+W(I)).LT.1.E-8_q)GOTO 100; ENDDO
        STOP'FAILED TO FIND PAIR!'
100     IPAIR(2,NPAIR)=IP
        FPAIR(1,2,NPAIR)=SUM(ABS(A(1:NG(1)  ,IP)))
        FPAIR(2,2,NPAIR)=SUM(ABS(A(1+NG(1):N,IP)))
      ENDDO ! I
      WRITE(IO,'(" ROBY: FOUND",I3," PAIRS WITH",I3," POSSIBLE LONE ORBS")')NPAIR,N-2*NPAIR
      DO I=1,NPAIR
      DO J=1,2
      OCC(J:J)=RESHAPE(MATMUL(TRANSPOSE(A(:,IPAIR(J,I):IPAIR(J,I))), &
                      &MATMUL(PIAJB,A(:,IPAIR(J,I):IPAIR(J,I)))),(/1/))
      ENDDO
      WRITE(IO,'(" PAIR OCC:",2F10.5," A/B COMP. FOR BONDING:",2F6.2, &
                  &" FOR ANTI-B:",2F6.2," E-VALUE=",F10.6)')OCC,FPAIR(:,:,I),EPAIR(I)
      CIAB=OCC(2)-OCC(1)+CIAB
      ENDDO
      RETURN

      END SUBROUTINE CALC_BOND_ORDER_ROBY_

!=======================================================================
      SUBROUTINE MATCH_TBS_SIJ(S,NG12,NQG,IQG,NQG12,CRG,IQNIG,ISP)
      IMPLICIT NONE
      INTEGER NQG12,ISP,NG12,NQG(2),CRG(3,NG12),IQG(NQG12),IQNIG(NQG12)
      REAL(q) S(NQG12,NQG12)
! LOCAL
      INTEGER I1,I2,J1,J2,NI1,NI2,IQ1,IQ2,IC,J1_,J2_
      INTEGER CR(3)

      J1_=0
      DO I1=1,2; DO J1=1,NQG(I1); J1_=J1_+1
      IQ1=IQG(J1_); NI1=IQNIG(J1_)
      J2_=0
      DO I2=1,2; DO J2=1,NQG(I2); J2_=J2_+1
      IQ2=IQG(J2_); NI2=IQNIG(J2_)
      CR=CRG(:,NI2)-CRG(:,NI1)
      CALL SEARCH_INDEX_IC(CR,IC)
      IF(IC.EQ.0)CYCLE
      S(J1_,J2_)=QB%SR(IQ1,IQ2,IC,ISP)
      ENDDO; ENDDO; ENDDO; ENDDO
      RETURN

      END SUBROUTINE MATCH_TBS_SIJ



!********************************************************************
      SUBROUTINE OUT_SPDM(PIAJB,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
      REAL(q) PIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)

      IF(.NOT.CL%LSPDM)RETURN
      WRITE(IO,'(" PRINT SINGLE PARTICLE DENSITY MATRIX ELEMENTS (BONG ORDER) FOR NN PAIRS WITH NIP-->:")')
      CALL PRINT_PAIR_ORB_INFO(PIAJB,IO,1)
      RETURN

      END SUBROUTINE OUT_SPDM

!********************************************************************
      SUBROUTINE CALC_MULLIKEN(PIAJB,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
      REAL(q) PIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)
! LOCAL
      INTEGER ISP,NI,NIP_,IQ1,IC
      REAL(q),ALLOCATABLE :: PSIAJB(:,:,:,:)

      IF(.NOT.CL%LMULLIKEN)RETURN
      ALLOCATE(PSIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN))
      PSIAJB=PIAJB*QB%SR
      WRITE(IO,*)
      WRITE(IO,'(" PRINT MULLIKEN FOR NN ORB PAIRS WITH JB-->:")')
      CALL PRINT_PAIR_ORB_INFO(PSIAJB,IO,2)
      WRITE(IO,*)
      WRITE(IO,'(" PRINT MULLIKEN FOR NN ATM PAIRS WITH J -->:")')
      CALL PRINT_PAIR_ATM_INFO(PSIAJB,IO)
      RETURN

      END SUBROUTINE CALC_MULLIKEN

!********************************************************************
      SUBROUTINE CALC_MAYER_BO(PIAJB,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
      REAL(q) PIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)
! LOCAL
      REAL(q),ALLOCATABLE :: PSIAJB(:,:,:,:),SPIAJB(:,:,:,:)

      IF(.NOT.CL%LMAYER)RETURN
      WRITE(*,'(" CALC_MAYER_BO...")')
      ALLOCATE(PSIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)); PSIAJB=0
      CALL CALC_XYIAJB(PIAJB,QB%SR,PSIAJB)
      ALLOCATE(SPIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)); SPIAJB=0
      CALL CALC_XYIAJB(QB%SR,PIAJB,SPIAJB)
      PSIAJB=PSIAJB*SPIAJB*2 ! 2--according to calibration of Qian
      DEALLOCATE(SPIAJB)
      WRITE(IO,*)
      WRITE(IO,'(" PRINT MAYER BO FOR NN ORB PAIRS WITH JB-->:")')
      CALL PRINT_PAIR_ORB_INFO(PSIAJB,IO,2)
      WRITE(IO,*)
      WRITE(IO,'(" PRINT MAYER BO FOR NN ATM PAIRS WITH J -->:")')
      CALL PRINT_PAIR_ATM_INFO(PSIAJB,IO)
      DEALLOCATE(PSIAJB)
      RETURN

      END SUBROUTINE CALC_MAYER_BO

!********************************************************************
      SUBROUTINE PRINT_PAIR_ORB_INFO(XIAJB,IO,MODE)
      USE prec
      IMPLICIT NONE
      REAL(q) XIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)
      INTEGER IO,MODE
! LOCAL
      INTEGER NI,NIP,NIP_,IQ1,IC,ISP
      REAL(q) XIA(QB%DIM,W%ISPIN),RES1,RES2
      
      XIA=0
      DO NI=1,T_INFO%NIONS
! 1NNs
      DO NIP_=1,YZ(NI)%SL(1)%N
      NIP=YZ(NI)%SL(1)%NI(NIP_); IC=YZ(NI)%SL(1)%IC(NIP_)
      WRITE(IO,'(" NI:",I3," NIP:",I3," DIST=",F6.2," CELL%R",3I3)')NI,NIP,YZ(NI)%SL(1)%DIST,CELL%R(:,IC)
      DO ISP=1,W%ISPIN
      WRITE(IO,'(" ISP=",I2)')ISP
      WRITE(IO,'(5X,16A7)')QB%LMCHAR(YZ(NIP)%IQMIN:YZ(NIP)%IQMAX)
      DO IQ1=YZ(NI)%IQMIN,YZ(NI)%IQMAX
        WRITE(IO,'(1X,A3,1X,16F10.5)')QB%LMCHAR(IQ1),XIAJB(IQ1,YZ(NIP)%IQMIN:YZ(NIP)%IQMAX,IC,ISP)
        XIA(IQ1,ISP)=XIA(IQ1,ISP)+SUM(XIAJB(IQ1,YZ(NIP)%IQMIN:YZ(NIP)%IQMAX,IC,ISP))
      ENDDO; ENDDO; ENDDO ! IQ1,ISP,NIP_
! 2NNs
      DO NIP_=1,YZ(NI)%SL(2)%N
      NIP=YZ(NI)%SL(2)%NI(NIP_); IC=YZ(NI)%SL(2)%IC(NIP_)
      WRITE(IO,'(" NI:",I3," NIP:",I3," DIST=",F6.2," CELL%R",3I3)')NI,NIP,YZ(NI)%SL(2)%DIST,CELL%R(:,IC)
      DO ISP=1,W%ISPIN
      WRITE(IO,'(" ISP=",I2)')ISP
      WRITE(IO,'(5X,16A7)')QB%LMCHAR(YZ(NIP)%IQMIN:YZ(NIP)%IQMAX)
      DO IQ1=YZ(NI)%IQMIN,YZ(NI)%IQMAX
        WRITE(IO,'(1X,A3,1X,16F10.5)')QB%LMCHAR(IQ1),XIAJB(IQ1,YZ(NIP)%IQMIN:YZ(NIP)%IQMAX,IC,ISP)
      ENDDO; ENDDO; ENDDO ! IQ1,ISP,NIP_

      IF(MODE.GT.1)THEN
      DO IQ1=YZ(NI)%IQMIN,YZ(NI)%IQMAX; DO ISP=1,W%ISPIN
      RES1=SUM(XIAJB(IQ1,:,:,ISP)); RES2=RES1-XIAJB(IQ1,IQ1,CELL%ORIG,ISP)
      WRITE(IO,'(" ISP=",I1," TOT:",F10.5," SHARING:",F10.5," SHARING WITH NNs:",F10.5)')ISP,RES1,RES2,XIA(IQ1,ISP)
      ENDDO
      RES1=SUM(XIAJB(IQ1,:,:,:))*W%RSPIN; RES2=RES1-SUM(XIAJB(IQ1,IQ1,CELL%ORIG,:))*W%RSPIN
      WRITE(IO,'(     " UP+DN TOT:",F10.5," SHARING:",F10.5," SHARING WITH NNs:",F10.5)')    RES1,RES2,SUM(XIA(IQ1,:))*W%RSPIN
      ENDDO ! IQ1
      ENDIF ! MODE.GT.1
      ENDDO ! NI
      RETURN

      END SUBROUTINE PRINT_PAIR_ORB_INFO

!********************************************************************
      SUBROUTINE PRINT_PAIR_ATM_INFO(XIAJB,IO)
      USE prec
      IMPLICIT NONE
      REAL(q) XIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)
      INTEGER IO
! LOCAL
      INTEGER NI,NIP,NIP_,IQ1,IQ2,IC,ISP
      REAL(q) XI(T_INFO%NIONS,W%ISPIN),RES1,RES2
      REAL(q),ALLOCATABLE::XIJ(:,:,:,:)

      ALLOCATE(XIJ(T_INFO%NIONS,T_INFO%NIONS,CELL%DIM,W%ISPIN)); XIJ=0
      DO NI =1,T_INFO%NIONS; DO IQ1=YZ(NI )%IQMIN,YZ(NI )%IQMAX
      DO NIP=1,T_INFO%NIONS; DO IQ2=YZ(NIP)%IQMIN,YZ(NIP)%IQMAX
        XIJ(NI,NIP,:,:)=XIJ(NI,NIP,:,:)+XIAJB(IQ1,IQ2,:,:)
      ENDDO; ENDDO; ENDDO; ENDDO
      XI=0
      DO NI=1,T_INFO%NIONS
      WRITE(IO,'(" NI=",I3)')NI
! 1NN
      WRITE(IO,'(" NN IDX(NJ):",20(I3,"(",I2")"))')(NIP,YZ(NI)%SL(1)%NI(NIP),NIP=1,YZ(NI)%SL(1)%N)
      DO ISP=1,W%ISPIN
      WRITE(IO,'(" ISP=",I2,20F10.5)')ISP,(XIJ(NI,YZ(NI)%SL(1)%NI(NIP),YZ(NI)%SL(1)%IC(NIP),ISP),NIP=1,YZ(NI)%SL(1)%N)
      DO NIP=1,YZ(NI)%SL(1)%N; XI(NI,ISP)=XI(NI,ISP)+XIJ(NI,YZ(NI)%SL(1)%NI(NIP),YZ(NI)%SL(1)%IC(NIP),ISP); ENDDO
      ENDDO ! ISP
      WRITE(IO,'("  UP+DN",20F10.5)')(SUM(XIJ(NI,YZ(NI)%SL(1)%NI(NIP),YZ(NI)%SL(1)%IC(NIP),:))*W%RSPIN,NIP=1,YZ(NI)%SL(1)%N)
! 2NN
      WRITE(IO,'(" 2NN IDX(NJ):",20(I3,"(",I2")"))')(NIP,YZ(NI)%SL(2)%NI(NIP),NIP=1,YZ(NI)%SL(2)%N)
      DO ISP=1,W%ISPIN
      WRITE(IO,'(" ISP=",I2,20F10.5)')ISP,(XIJ(NI,YZ(NI)%SL(2)%NI(NIP),YZ(NI)%SL(2)%IC(NIP),ISP),NIP=1,YZ(NI)%SL(2)%N)
!      DO NIP=1,YZ(NI)%SL(2)%N; XI(NI,ISP)=XI(NI,ISP)+XIJ(NI,YZ(NI)%SL(2)%NI(NIP),YZ(NI)%SL(2)%IC(NIP),ISP); ENDDO
      ENDDO ! ISP
      WRITE(IO,'("  UP+DN",20F10.5)')(SUM(XIJ(NI,YZ(NI)%SL(2)%NI(NIP),YZ(NI)%SL(2)%IC(NIP),:))*W%RSPIN,NIP=1,YZ(NI)%SL(2)%N)

      DO ISP=1,W%ISPIN
      RES1=SUM(XIJ(NI,:,:,ISP)); RES2=RES1-XIJ(NI,NI,CELL%ORIG,ISP)
      WRITE(IO,'(" ISP=",I1," TOT:",F10.5," SHARING:",F10.5," SHARING WITH NNs:",F10.5)')ISP,RES1,RES2,XI(NI,ISP)
      ENDDO ! ISP
      RES1=SUM(XIJ(NI,:,:,:))*W%RSPIN; RES2=RES1-SUM(XIJ(NI,NI,CELL%ORIG,:))*W%RSPIN
      WRITE(IO,'(     " UP+DN TOT:",F10.5," SHARING:",F10.5," SHARING WITH NNs:",F10.5)')    RES1,RES2,SUM(XI(NI,:))*W%RSPIN
      ENDDO ! NI
      WRITE(IO,'(" SUM-RULE:",F10.5," NET=",F10.5)')SUM(XIAJB)*W%RSPIN,DOS%NET
      RETURN

      END SUBROUTINE PRINT_PAIR_ATM_INFO

!=======================================================================
      SUBROUTINE CALC_FULL_PIAJB(PIAJB)
      USE prec
      IMPLICIT NONE
      REAL(q) PIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)
! LOCAL
      INTEGER IC,IQ1,IQ2,ISP

      DO IQ1=1,QB%DIM; DO IQ2=1,QB%DIM; DO IC=1,CELL%DIM; DO ISP=1,W%ISPIN
      CALL CALC_PIAJB_(PIAJB(IQ1,IQ2,IC,ISP),W%CPTWFP(:,:,:,ISP),CELL%R(:,IC),IQ1,IQ2,ISP)
      ENDDO; ENDDO; ENDDO; ENDDO
      RETURN

      END SUBROUTINE CALC_FULL_PIAJB

!=======================================================================
      SUBROUTINE CALC_SELECTIVE_PIAJB(PIAJB,NQG12,IQG,IQNIG,CRG,NQG,NG12,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NQG12,ISP,NG12,NQG(2),IQG(NQG12),IQNIG(NQG12),CRG(3,NG12)
      REAL(q) PIAJB(NQG12,NQG12)
! LOCAL
      INTEGER I1,I2,J1,J2,NI1,NI2,J1_,J2_,IQ1,IQ2,CR(3),NK
      COMPLEX(q),ALLOCATABLE::WQO(:,:,:)

      ALLOCATE(WQO(QB%DIM,W%NBANDS,KPOINTS%NKPTS))
      WQO=W%CPTWFP(:,:,:,ISP)
      DO NK=1,KPOINTS%NKPTS
        CALL ZGETRI_(WQO(:,:,NK),QB%DIM); WQO(:,:,NK)=CONJG(TRANSPOSE(WQO(:,:,NK)))
      ENDDO ! NK
      J1_=0
      DO I1=1,2; DO J1=1,NQG(I1); J1_=J1_+1; NI1=IQNIG(J1_); IQ1=IQG(J1_)
      J2_=0
      DO I2=1,2; DO J2=1,NQG(I2); J2_=J2_+1; NI2=IQNIG(J2_); IQ2=IQG(J2_)
      CR=CRG(:,NI2)-CRG(:,NI1)
      CALL CALC_PIAJB_(PIAJB(J1_,J2_),WQO,CR,IQ1,IQ2,ISP)
      ENDDO; ENDDO; ENDDO; ENDDO
      RETURN

      END SUBROUTINE CALC_SELECTIVE_PIAJB

!=======================================================================
      SUBROUTINE CALC_PIAJB_(PIAJB_,WQO,CR,IQ1,IQ2,ISP)
      USE constant; USE prec
      IMPLICIT NONE
      REAL(q) PIAJB_
      INTEGER CR(3),IQ1,IQ2,ISP
      COMPLEX(q) WQO(QB%DIM,W%NBANDS,KPOINTS%NKPTS)
! LOCAL
      INTEGER NK,NB
      COMPLEX(q) CC,PHASEK(KPOINTS%NKPTS)

      PIAJB_=0
      DO NK=1,KPOINTS%NKPTS; PHASEK(NK)=EXP(-CITPI*SUM(CR*KPOINTS%VKPT(:,NK))); ENDDO
      CC=0
      DO NB=1,W%NBANDS
      CC=CC+SUM(W%FERWE(NB,:,ISP)*KPOINTS%WTKPT*CONJG(WQO(IQ1,NB,:)) &
          &*WQO(IQ2,NB,:)*PHASEK)
      ENDDO ! NB
      PIAJB_=REAL(CC,q)
      RETURN

      END SUBROUTINE CALC_PIAJB_

!=======================================================================
      SUBROUTINE CALC_XYIAJB(XIAJB,YIAJB,XYIAJB)
      USE prec
      IMPLICIT NONE
      REAL(q) XIAJB (QB%DIM,QB%DIM,CELL%DIM,W%ISPIN),YIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)
      REAL(q) XYIAJB(QB%DIM,QB%DIM,CELL%DIM,W%ISPIN)
! LOCAL
      INTEGER NI1,NI2,IS1,NI3,NI2_,NI3_,IC3,IS3,IC2,IQ1,IQ2,ISP,IC2T
      INTEGER CR(3)

      XYIAJB=0
      DO NI1=1,T_INFO%NIONS
      DO IS1=0,YZ(NI1)%NSL; DO NI3_=1,YZ(NI1)%SL(IS1)%N
      NI3=YZ(NI1)%SL(IS1)%NI(NI3_); IC3=YZ(NI1)%SL(IS1)%IC(NI3_)
      DO IS3=0,YZ(NI3)%NSL; DO NI2_=1,YZ(NI3)%SL(IS3)%N
      NI2=YZ(NI3)%SL(IS3)%NI(NI2_); IC2=YZ(NI3)%SL(IS3)%IC(NI2_)
      CR=CELL%R(:,IC3)+CELL%R(:,IC2)
      CALL SEARCH_INDEX_IC(CR,IC2T)
      IF(IC2T.EQ.0)CYCLE
      DO IQ1=YZ(NI1)%IQMIN,YZ(NI1)%IQMAX; DO IQ2=YZ(NI2)%IQMIN,YZ(NI2)%IQMAX; DO ISP=1,W%ISPIN
        XYIAJB(IQ1,IQ2,IC2T,ISP)=XYIAJB(IQ1,IQ2,IC2T,ISP)+ &
                                &SUM(XIAJB(IQ1,YZ(NI3)%IQMIN:YZ(NI3)%IQMAX,IC3,ISP)* &
                                &    YIAJB(YZ(NI3)%IQMIN:YZ(NI3)%IQMAX,IQ2,IC2,ISP))
      ENDDO; ENDDO; ENDDO
      ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
      RETURN

      END SUBROUTINE CALC_XYIAJB

!=======================================================================
      SUBROUTINE SEARCH_INDEX_IC(CR,IC)
      IMPLICIT NONE
      INTEGER CR(3),IC
! LOCAL
      INTEGER I

      IC=0
      DO I=1,CELL%DIM
      IF(CR(1).NE.CELL%R(1,I))CYCLE
      IF(CR(2).NE.CELL%R(2,I))CYCLE
      IF(CR(3).NE.CELL%R(3,I))CYCLE
      IC=I; EXIT
      ENDDO
      RETURN

      END SUBROUTINE SEARCH_INDEX_IC


      END MODULE TB_MODULE

! OPEN SUBROUTINES
!**************** SUBROUTINE EXPRO   ***********************************
! EXPRO: caclulates the x-product of two vectors
!***********************************************************************
      SUBROUTINE EXPRO(H,U1,U2)
      USE prec
      IMPLICIT NONE
      REAL(q) H(3),U1(3),U2(3)

      H(1)=U1(2)*U2(3)-U1(3)*U2(2)
      H(2)=U1(3)*U2(1)-U1(1)*U2(3)
      H(3)=U1(1)*U2(2)-U1(2)*U2(1)

      RETURN
      END SUBROUTINE

!***********************SUBROUTINE DENNP *******************************
! if ISMEAR=0
! subroutine DENSTA calculates a continuous density of states in the
! interval (EMIN,EMAX) by applying a gaussian broadening to the discrete
! eigenvalue spectrum contained in CELEN(NBANDS,NKPTS). The width of the
! gaussians is SIGMA. The fermi energy EFERMI is calculated from the
! integrated dos
! correction to the variational energy is calculated (-SIGMA S)
! according to A.de Vita
! if ISMEAR>0 the generalized form of Methfessel and Paxton of order
!        N=ISMEAR will be used instead of Gaussians to get the dos ...
! routine is parallelized to get full speed ...
! initially it took 2 seconds to calculate occupancies (gK)
!***********************************************************************
      SUBROUTINE DENMP(ISPIN,RSPIN,EMIN,EMAX,NELECT,EFERMI, &
             & ISMEAR,SIGMA,NBANDS,NKPTS,FERWE,CELEN, &
             & NEDOS,DOS,DOSI,PAR,DOSPAR,WTKPT)
      USE prec; USE constant
      IMPLICIT NONE
      INTEGER ISPIN,RSPIN,ISMEAR,NBANDS,NKPTS,NEDOS
      REAL(q) DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN)
      REAL(q) NELECT,EMIN,EMAX,EFERMI,SIGMA
      REAL(q) FERWE(NBANDS,NKPTS,ISPIN),CELEN(NBANDS,NKPTS,ISPIN)
      REAL(q) PAR(NBANDS,NBANDS,NKPTS,ISPIN),DOSPAR(NEDOS,NBANDS,ISPIN)
      REAL(q) WTKPT(NKPTS)
! Local
      INTEGER ISP,K,N,NELOW,NEHIG,I,NITER
      REAL(q) DELTAE,EPS,WEIGHT,SFUN_DONE,E,EPSDOS,DFUN,SFUN,DOSTOT
      REAL(q) EF1,EF2,ELECT,X1
      LOGICAL LOWB,HIGHB

      DELTAE=(EMAX-EMIN)/(NEDOS-1)
      DOS =0; DOSI=0; DOSPAR=0
!=======================================================================
! accumulate dos and integrated dos
!=======================================================================
      DO ISP=1,ISPIN
      DO K=1,NKPTS
      DO N=1,NBANDS
        EPS=CELEN(N,K,ISP)
        WEIGHT= RSPIN*WTKPT(K)
        NELOW=(EPS-8._q*SIGMA-EMIN)/DELTAE+1
        NEHIG=(EPS+8._q*SIGMA-EMIN)/DELTAE+1
        IF (NELOW<1)     NELOW=1
        IF (NELOW>NEDOS) NELOW=NEDOS
        IF (NEHIG<1)     NEHIG=1
        IF (NEHIG>NEDOS) NEHIG=NEDOS

        SFUN_DONE=0
        DO I=NELOW,NEHIG
          E=EMIN+DELTAE*(I-1)-EPS
          CALL DELSTP(ISMEAR,(E/SIGMA),DFUN,SFUN)
          EPSDOS=DFUN/SIGMA
!gK fix the DOS so that the integrated DOS yields accurate results
          EPSDOS=(SFUN-SFUN_DONE)/DELTAE
          SFUN_DONE=SFUN

          DOS(I,ISP) =DOS(I,ISP) +(WEIGHT*EPSDOS)
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT*SFUN
          DOSPAR(I,:,ISP)=DOSPAR(I,:,ISP)+(WEIGHT*EPSDOS)*PAR(:,N,K,ISP)
        ENDDO
        DO I=NEHIG+1,NEDOS
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT
        ENDDO
      ENDDO; ENDDO; ENDDO
!=======================================================================
! calculate approximated fermi energy
!=======================================================================
      DO I=1,NEDOS
        DOSTOT=SUM(DOSI(I,:))
        IF (ABS(DOSTOT-NELECT)<0.01_q .OR.DOSTOT>NELECT) EXIT
      ENDDO
      EFERMI= EMIN+(I-1)*DELTAE
      IF (SIGMA<1E-5_q) RETURN
!=======================================================================
! search now for exact Fermi-level using bisectioning
!=======================================================================
      EF1= EMIN+(I-2)*DELTAE
      LOWB =.FALSE.
      EF2= EMIN+(I-1)*DELTAE
      HIGHB=.FALSE.
      NITER=0

      setfermi: DO

      EFERMI=(EF1+EF2)/2
      NITER=NITER+1

      ELECT=0
      DO ISP=1,ISPIN
      DO K=1,NKPTS
      DO N=1,NBANDS
        EPS=CELEN(N,K,ISP)
        X1=(EFERMI-EPS)/SIGMA
        CALL DELSTP(ISMEAR,X1,DFUN,SFUN)
        FERWE(N,K,ISP)=SFUN
        ELECT=ELECT+FERWE(N,K,ISP)*WTKPT(K)
      ENDDO; ENDDO; ENDDO
      ELECT=ELECT*RSPIN

      ! compare with number of electrons

      IF ( ABS(ELECT-NELECT)<1E-10_q) RETURN
      IF ( (ABS(EF1-EF2)/(ABS(EFERMI)+1.E-10_q))<1E-14_q) EXIT
      IF ( ELECT>NELECT) THEN
        IF (.NOT.LOWB)  EF1=EF1-DELTAE
        HIGHB=.TRUE.
        EF2  =EFERMI
      ELSE
        IF (.NOT.HIGHB) EF2=EF2+DELTAE
        LOWB =.TRUE.
        EF1  =EFERMI
      ENDIF
      ENDDO setfermi
      WRITE(0,*)'WARNING: DENMP: can''t reach specified precision'
      WRITE(0,*)'Number of Electrons is NELECT =',ELECT

      END SUBROUTINE DENMP

!******************** DELSTP    ****************************************
! Returns generalised delta and step functions (Methfessel & Paxton)
!
!  Input:
!      n > -1 : order of approximant; x : argument
!  Output:
!      D_n (x) ,  S_n (x)
!  Remarks:
!      D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
!      S_n (x) = (1 - erf x)/2 + exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
!      where H is a Hermite polynomial and
!      A_i = (-1)^i / ( i! 4^i sqrt(pi) )
!***********************************************************************
      SUBROUTINE DELSTP(N,X,D,S)
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)

      IF (X<-1.E5_q) THEN
         D=0._q
         S=0._q
         RETURN
      END IF
      IF (X>1.E5_q) THEN
         D=0._q
         S=1._q
         RETURN
      END IF
!=======================================================================
!  If n < 0 : assume Gaussian type smearing
!  (must return  same as N=0 or ... )
!=======================================================================
      IF (N<0) THEN
         D=EXP(-(X*X))/SQRT(PI)
         S=0.5_q+0.5_q*ERRF(X)
         RETURN
      END IF
!=======================================================================
! Methfessel & Paxton
!=======================================================================
      EX2=EXP(-(X*X))
      S0=0.5_q*ERRF(X)
      A=1._q/SQRT(PI)
      K=0
      H1=1._q
      H2=2._q*X
      S=0._q
      D=A
      DO I=1,N
         A=A/((-4._q)*I)
         K=K+1
         H3=H1
         H1=H2
         H2=2._q*X*H2-2*K*H3
         S=S+A*H1
         K=K+1
         H3=H1
         H1=H2
         H2=2._q*X*H2-2*K*H3
         D=D+A*H1
      ENDDO
      D=D*EX2
      S=0.5_q+S0-S*EX2
      RETURN
      END SUBROUTINE DELSTP

!**************** SUBROUTINE KARDIR ************************************
! transform a set of vectors from cartesian coordinates to
! ) direct lattice      (BASIS must be equal to B reciprocal lattice)
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************
      SUBROUTINE KARDIR(NMAX,V,BASIS)
      USE prec
      IMPLICIT NONE
      INTEGER N,NMAX
      REAL(q) V(3,NMAX),BASIS(3,3),V1,V2,V3

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(2,1)+V(3,N)*BASIS(3,1)
        V2=V(1,N)*BASIS(1,2)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(3,2)
        V3=V(1,N)*BASIS(1,3)+V(2,N)*BASIS(2,3)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE

!**************** SUBROUTINE TOPRIM ************************************
! bring all ions into the primitive cell
!***********************************************************************
      SUBROUTINE TOPRIM(NIONS,POSION)
      USE prec
      IMPLICIT NONE
      INTEGER NIONS,I
      REAL(q) POSION(3,NIONS)

      DO I=1,NIONS
      POSION(1,I)=MOD(POSION(1,I)+60,1._q)
      POSION(2,I)=MOD(POSION(2,I)+60,1._q)
      POSION(3,I)=MOD(POSION(3,I)+60,1._q)
      ENDDO
      RETURN
      END SUBROUTINE

!***********************************************************************
!                                                                      *
      SUBROUTINE ORDER(N,ARRIN,INDEX)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
!                                                                      *
!***********************************************************************
!                                                                      *
!   Routine ORDER indexes a real array by ascending order. The scheme  *
!   bases on a heapsort algorithm which is described in the book of    *
!   W.H.PRESS et al.: NUMERICAL RECIPIES, CAMBRIDGE UNIV. PR., 1987.   *
!                                                                      *
!                                                                      *
!   Input parameters:                                                  *
!   -----------------                                                  *
!                                                                      *
!      N contains the dimension of array ARRIN (number of data).       *
!      ARRIN contains the input numbers to be sorted (to be indexed).  *
!                                                                      *
!                                                                      *
!   Output parameters:                                                 *
!   ------------------                                                 *
!                                                                      *
!      INDEX contains the permutation table for the ordered numbers.   *
!                                                                      *
!                                                                      *
!***********************************************************************

      DIMENSION ARRIN(N),INDEX(N)

      DO 1 J=1,N
         INDEX(J)=J
    1 CONTINUE
      IF (N==1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
      IF (L>1) THEN
         L=L-1
         INDEXT=INDEX(L)
         QQ=ARRIN(INDEXT)
      ELSE
         INDEXT=INDEX(IR)
         QQ=ARRIN(INDEXT)
         INDEX(IR)=INDEX(1)
         IR=IR-1
         IF (IR==1)THEN
            INDEX(1)=INDEXT
            RETURN
         END IF
      END IF
      I=L
      J=L+L
20    IF (J<=IR) THEN
         IF (J<IR) THEN
            IF(ARRIN(INDEX(J))<ARRIN(INDEX(J+1))) J=J+1
         END IF
         IF (QQ<ARRIN(INDEX(J))) THEN
            INDEX(I)=INDEX(J)
            I=J
            J=J+J
         ELSE
            J=IR+1
         END IF
         GOTO 20
      END IF
      INDEX(I)=INDEXT
      GOTO 10
      END SUBROUTINE ORDER

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

!***********************************************************************
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

!***********************************************************************
! Estimate the maximal rcut due to finite k-grids by scanning six-planes
!***********************************************************************
      SUBROUTINE CALC_MAX_RCUT(N,CR,A,MAXRCUT)
      USE prec
      IMPLICIT NONE
      INTEGER N,CR(3,N)
      REAL(q) A(3,3),MAXRCUT
! LOCAL
      INTEGER NPT,RPT(3,7),I,NDIM(3),J1,J2,J3,IPT
      
      NPT=1; RPT=0
      DO I=1,3
        NDIM(I)=MAXVAL(CR(I,:))-MINVAL(CR(I,:))+1
        IF(NDIM(I).GT.1)NPT=NPT*2
      ENDDO
      NPT=NPT-1
      IF(NPT.EQ.0)NPT=7
      IPT=0
      DO J1=0,NDIM(1),NDIM(1); DO J2=0,NDIM(2),NDIM(2); DO J3=0,NDIM(3),NDIM(3)
      IF(J1.EQ.0.AND.J2.EQ.0.AND.J3.EQ.0)CYCLE
      IF(NPT.NE.7.AND.(J1.EQ.1.OR.J2.EQ.1.OR.J3.EQ.1))CYCLE
      IPT=IPT+1
      RPT(1,IPT)=J1; RPT(2,IPT)=J2; RPT(3,IPT)=J3
      ENDDO; ENDDO; ENDDO
      IF(IPT.NE.NPT)THEN
        WRITE(0,'(" IPT/NPT:",2I2)')IPT,NPT
        STOP' COUNTING ERROR IN CALC_MAX_RCUT!'
      ENDIF

      MAXRCUT=100._q
      DO IPT=1,NPT; MAXRCUT=MIN(MAXRCUT,SQRT(SUM(MATMUL(A,RPT(:,IPT))**2))/2); ENDDO
      RETURN

      END SUBROUTINE CALC_MAX_RCUT
