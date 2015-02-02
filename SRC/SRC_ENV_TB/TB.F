#include "symbol.inc"
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
         LOGICAL :: LROT
         REAL(q) :: MROT(3,3)
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
        INTEGER, POINTER :: ISORT(:)  ! SORT ACCORDING TO LOCAL ENVIRONMENT
        INTEGER, POINTER :: GJ(:)     ! The nationality of each atom
      END TYPE type_info

      TYPE kpoints_struct
!only  KPOINTS,KPT_BAND,KPT_FS
        INTEGER NKPTS,NKPT(3)        ! actual number of k-points
        REAL(q),POINTER:: VKPT(:,:)  ! coordinate of k-point
        REAL(q),POINTER:: WTKPT(:)   ! symmetry weight-factor for each k-point
        INTEGER ISMEAR               ! type of smearing
        REAL(q) SIGMA                ! type of smearing
        INTEGER NEDOS
        REAL(q) EMIN                 ! minimal E for DOS
        REAL(q) EMAX                 ! maximal E for DOS
        REAL(q),POINTER :: DOS(:,:),DOSI(:,:)
        CHARACTER*40  SZNAMK         ! name of k-points file
      END TYPE kpoints_struct

      TYPE wavespin
!only W
        INTEGER NBANDS,ISPIN,RSPIN
        REAL(q) EFERMI,NET
        REAL(q),   POINTER:: FERWE(:,:,:)    ! fermi-weight for each band
        REAL(q),   POINTER:: CELEN(:,:,:)    ! eigenvalues
      END TYPE wavespin

      TYPE CELL_TYPE
        INTEGER DIM,N(3)
        INTEGER,ALLOCATABLE :: R(:,:)
      END TYPE CELL_TYPE

      TYPE QO_TYPE
        INTEGER NSP  ! ISPIN FOR TBH/S
        INTEGER DIM,ID
        REAL(q),POINTER :: SR(:,:,:,:,:,:)
        REAL(q),POINTER :: HR(:,:,:,:,:,:)
        COMPLEX(q),POINTER ::  SK(:,:,:,:)
        COMPLEX(q),POINTER ::  HK(:,:,:,:)
      END TYPE QO_TYPE

      TYPE ENQB_TYPE
        INTEGER NSORT,NGJ ! SORT/NGJ OF ENV-DEPENDENT ATOMS
        REAL(q) NNN,RNN,RCUT ! N.N. CUT OFF RADIUS
        INTEGER,POINTER :: NR12(:,:,:,:) ! NUMBER OF POSITION VECTORS FOR EACH PAIR
        REAL(q),POINTER ::  R12(:,:,:,:,:,:) !  POSITION VECTORS
        TYPE(QO_TYPE) QO
      END TYPE ENQB_TYPE

! VARIABLES
      TYPE (latt)           :: LATT_CUR
      TYPE (type_info)      :: T_INFO
      TYPE (kpoints_struct) :: KPOINTS,KPT_BAND,KPT_FS
      TYPE (wavespin)       :: W
      TYPE (CELL_TYPE)      :: CELL
      TYPE (QO_TYPE)        :: QB
      TYPE (ENQB_TYPE)      :: ENQB
      TYPE (QO_TYPE),ALLOCATABLE :: YZD(:,:) ! Atom pair
      TYPE (QO_TYPE),ALLOCATABLE :: ENYZD(:,:,:,:)

      CONTAINS

! PUBLIC SUBROUTINES
!***********************************************************************
      SUBROUTINE INI_TB(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO

      CALL RD_ENTB(IU,IO)
      CALL SET_ENYZ_PAIR()
      CALL RD_POSCAR(IU,IO)
      CALL RD_ROTCAR(IU)
      CALL SORT_ATOM(IO)
      CALL SET_NET(IO)
      CALL RD_INCAR(IU,IO)
      CALL SET_QB_YZD(IO)
      CALL RD_KPOINTS(KPOINTS,'KPOINTS',IU,IO)
      CALL SET_BAND()
      CALL CALC_FULL_BAND(IU,IO)
      CALL SET_FERMI_WT(IO)
      CALL OUT_TDOS(1)
      RETURN

      END SUBROUTINE INI_TB


! PRIVATE SUBROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_ENTB(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      INTEGER I1,I2,ISP,IR,NBASE,GJ1,GJ2
      INTEGER,ALLOCATABLE::IQ2SORT(:)

      WRITE(*,'(" RD_ENTB...")')
      OPEN(IU,FILE='ENTB.INP',STATUS='OLD')
      READ(IU,*)ENQB%QO%NSP
      READ(IU,*)ENQB%NGJ
      READ(IU,*)ENQB%NSORT
      ALLOCATE(ENYZD(ENQB%NSORT,ENQB%NGJ,ENQB%NSORT,ENQB%NGJ))
      DO I1=1,ENQB%NSORT; READ(IU,*)ENYZD(I1,1,1,1)%ID,ENYZD(I1,1,1,1)%DIM; ENDDO
      ENQB%QO%DIM=SUM(ENYZD(:,1,1,1)%DIM)
      ALLOCATE(IQ2SORT(ENQB%QO%DIM))
      NBASE=0
      DO I1=1,ENQB%NSORT
        IQ2SORT(NBASE+1:NBASE+ENYZD(I1,1,1,1)%DIM)=I1; NBASE=NBASE+ENYZD(I1,1,1,1)%DIM
      ENDDO ! I1
      ALLOCATE(ENQB%NR12(ENQB%NSORT,ENQB%NGJ,ENQB%NSORT,ENQB%NGJ))
      DO GJ1=1,ENQB%NGJ; DO I1=1,ENQB%NSORT; DO GJ2=1,ENQB%NGJ; DO I2=1,ENQB%NSORT
      READ(IU,*)ENQB%NR12(I1,GJ1,I2,GJ2)
      ENDDO; ENDDO; ENDDO; ENDDO
      ALLOCATE(ENQB%R12(0:3,ENQB%NSORT,ENQB%NGJ,ENQB%NSORT,ENQB%NGJ,MAXVAL(ENQB%NR12)))
      DO GJ1=1,ENQB%NGJ; DO I1=1,ENQB%NSORT; DO GJ2=1,ENQB%NGJ; DO I2=1,ENQB%NSORT
      IF(ENQB%NR12(I1,GJ1,I2,GJ2).EQ.0)CYCLE
      READ(IU,*)
      READ(IU,*)ENQB%R12(0:3,I1,GJ1,I2,GJ2,1:ENQB%NR12(I1,GJ1,I2,GJ2))
      ENDDO; ENDDO; ENDDO; ENDDO
      ALLOCATE(ENQB%QO%HR(ENQB%QO%DIM,ENQB%NGJ,ENQB%QO%DIM,ENQB%NGJ,MAXVAL(ENQB%NR12),ENQB%QO%NSP)); ENQB%QO%HR=0
      DO ISP=1,ENQB%QO%NSP
      DO GJ1=1,ENQB%NGJ; DO I1=1,ENQB%QO%DIM; DO GJ2=1,ENQB%NGJ; DO I2=1,ENQB%QO%DIM
      IF(ENQB%NR12(IQ2SORT(I1),GJ1,IQ2SORT(I2),GJ2).EQ.0)CYCLE
      READ(IU,*)
      DO IR=1,ENQB%NR12(IQ2SORT(I1),GJ1,IQ2SORT(I2),GJ2)
        READ(IU,*)ENQB%QO%HR(I1,GJ1,I2,GJ2,IR,ISP)
      ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
      ENDDO ! ISP
      ALLOCATE(ENQB%QO%SR(ENQB%QO%DIM,ENQB%NGJ,ENQB%QO%DIM,ENQB%NGJ,MAXVAL(ENQB%NR12),ENQB%QO%NSP)); ENQB%QO%SR=0
      DO ISP=1,ENQB%QO%NSP
      DO GJ1=1,ENQB%NGJ; DO I1=1,ENQB%QO%DIM; DO GJ2=1,ENQB%NGJ; DO I2=1,ENQB%QO%DIM
      IF(ENQB%NR12(IQ2SORT(I1),GJ1,IQ2SORT(I2),GJ2).EQ.0)CYCLE
      READ(IU,*)
      DO IR=1,ENQB%NR12(IQ2SORT(I1),GJ1,IQ2SORT(I2),GJ2)
        READ(IU,*)ENQB%QO%SR(I1,GJ1,I2,GJ2,IR,ISP)
      ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
      ENDDO ! ISP
      CLOSE(IU)
      DEALLOCATE(IQ2SORT)
      WRITE(*,'(" DONE.")')
      RETURN

      END SUBROUTINE RD_ENTB

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_ROTCAR(IU)
      USE prec
      IMPLICIT NONE
      INTEGER IU
! LOCAL
      INTEGER I

      INQUIRE(FILE='ROTCAR',EXIST=LATT_CUR%LROT)
      IF(LATT_CUR%LROT)THEN
      WRITE(*,'(" ADDITIONAL COOR ROTATION (ROTCAR) WOULD BE TAKEN INTO CONSIDERATION!")')
      OPEN(IU,FILE='ROTCAR',STATUS='OLD')
      READ(IU,*)
      READ(IU,*)(LATT_CUR%MROT(I,:),I=1,3)
      CLOSE(IU)
      ENDIF
      RETURN

      END SUBROUTINE RD_ROTCAR

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SET_ENYZ_PAIR()
      USE prec
      IMPLICIT NONE
! LOCAL
      INTEGER I1,I2,NBASE1,NBASE2,GJ1,GJ2

      WRITE(*,'(" SET_ENYZ_PAIR...")')
      ENQB%NNN=4; ENQB%RNN=2.8_q; ENQB%RCUT=8.5_q
      NBASE1=0
      DO I1=1,ENQB%NSORT
      NBASE2=0
      DO I2=1,ENQB%NSORT
      DO GJ1=1,ENQB%NGJ; DO GJ2=1,ENQB%NGJ
      ENYZD(I1,GJ1,I2,GJ2)%HR=>ENQB%QO%HR(NBASE1+1:NBASE1+ENYZD(I1,1,1,1)%DIM,GJ1:GJ1,NBASE2+1:NBASE2+ENYZD(I2,1,1,1)%DIM,GJ2:GJ2,:,:)
      ENYZD(I1,GJ1,I2,GJ2)%SR=>ENQB%QO%SR(NBASE1+1:NBASE1+ENYZD(I1,1,1,1)%DIM,GJ1:GJ1,NBASE2+1:NBASE2+ENYZD(I2,1,1,1)%DIM,GJ2:GJ2,:,:)
      ENDDO; ENDDO
      NBASE2=NBASE2+ENYZD(I2,1,1,1)%DIM
      ENDDO
      NBASE1=NBASE1+ENYZD(I1,1,1,1)%DIM
      ENDDO
      WRITE(*,'(" DONE.")')
      RETURN

      END SUBROUTINE SET_ENYZ_PAIR

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_POSCAR(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      CHARACTER*255  INPLIN,INPWRK
      INTEGER, EXTERNAL :: NITEMS
      CHARACTER*1  CSEL
      INTEGER I,NT,NI,NSCALE
      REAL(q) SCALEX,SCALEY,SCALEZ

      WRITE(*,'(" RD_POSCAR...")')
      OPEN(IU,FILE='POSCAR',STATUS='OLD')
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
      ALLOCATE(T_INFO%ITYP(T_INFO%NIONS),T_INFO%GJ(T_INFO%NIONS))
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
        READ(IU,*,ERR=400,END=400) T_INFO%POSION(:,NI),T_INFO%GJ(NI)
      ENDDO

      IF (CSEL=='K') THEN
        T_INFO%POSION(1,:)=LATT_CUR%SCALE*T_INFO%POSION(1,:)*SCALEX
        T_INFO%POSION(2,:)=LATT_CUR%SCALE*T_INFO%POSION(2,:)*SCALEY
        T_INFO%POSION(3,:)=LATT_CUR%SCALE*T_INFO%POSION(3,:)*SCALEZ
        CALL KARDIR(T_INFO%NIONS,T_INFO%POSION,LATT_CUR%B)
      ENDIF
!      CALL TOPRIM(T_INFO%NIONS,T_INFO%POSION)
      CLOSE(IU)
      WRITE(*,'(" DONE.")')
      RETURN
 400 CONTINUE
      WRITE(0,*)' No initial positions read in'
      STOP

      END SUBROUTINE RD_POSCAR

!=======================================================================
      SUBROUTINE SET_NET(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER NT
      INTEGER,PARAMETER::NE(3)=(/11,3,6/)

      W%NET=0
      DO NT=1,T_INFO%NTYP
        W%NET=W%NET+NE(NT)*T_INFO%NITYP(NT)
      ENDDO
      WRITE(IO,'(" NET=",F8.1)')W%NET
      WRITE(* ,'(" NET=",F8.1)')W%NET
      RETURN

      END SUBROUTINE SET_NET

!=======================================================================
      SUBROUTINE SORT_ATOM(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER NI,NIP,ICX,ICY,ICZ,INN,IL,ID,NNN,NMAX1,NMAX2
      INTEGER I1,I2,J1,J2,J3,J4,NT,NTP
      REAL(q) DPOS(3),SCORE,DIST,RNN,MINDEV
      REAL(q),PARAMETER::PRIMN(46)=(/2,   3,  5,    7,  11,  13,  17,  19,  23,  29, &
                                  & 31,  37,  41,  43,  47,  53,  59,  61,  67,  71, &
                                  & 73,  79,  83,  89,  97, 101, 103, 107, 109, 113, &
                                  &127, 131, 137, 139, 149, 151, 157, 163, 167, 173, &
                                  &179, 181, 191, 193, 197, 199/)
      INTEGER,ALLOCATABLE::IDX(:)
      REAL(q),ALLOCATABLE::SLIST(:)

      WRITE(*,'(" SORT_ATOM...")')
      NNN=ENQB%NNN; RNN=ENQB%RNN
      ALLOCATE(T_INFO%ISORT(T_INFO%NIONS))
! ! 1ST CYCLE
! Score list to see no duplications.
      NMAX1=(4*3)/2*(2*1  )/2    + & ! 2Cu+2In: 6
           &(4)    *(3*2*1)/(3*2)+ & ! 1Cu+3In: 4
           &(4)    *(3*2  )/2    + & ! 1Cu+2In: 12
           &(4)    *(3*2*1)/(3*2)    !     3IN: 4
      WRITE(IO,'(" MAXIMAL SORT OF Se In CYCLE I:",I6)')NMAX1
      IF(NMAX1.GT.46)STOP'NEED MORE PRIME NUMBERS!'
      ALLOCATE(SLIST(NMAX1),IDX(NMAX1)); IL=0
      DO I1=1,NMAX1; IDX(I1)=I1; ENDDO
! 2Cu+2In
      DO I1=1,4; DO I2=I1+1,4
      DO J1=   1,4; IF(J1.NE.I1.AND.J1.NE.I2)EXIT; ENDDO
      DO J2=J1+1,4; IF(J2.NE.I1.AND.J2.NE.I2)EXIT; ENDDO
      IL=IL+1; SLIST(IL)=PRIMN(I1)*PRIMN(I2)*PRIMN(J1+4)*PRIMN(J2+4)
      ENDDO; ENDDO
! 1Cu+3In
      DO I1=1,4
      DO J1=   1,4; IF(J1.NE.I1)EXIT; ENDDO
      DO J2=J1+1,4; IF(J2.NE.I1)EXIT; ENDDO
      DO J3=J2+1,4; IF(J3.NE.I1)EXIT; ENDDO
      IL=IL+1; SLIST(IL)=PRIMN(I1)*PRIMN(J1+4)*PRIMN(J2+4)*PRIMN(J3+4)
      ENDDO
! 1Cu+2In
      DO I1=1,4
      DO J1=   1,4; IF(J1.EQ.I1)CYCLE
      DO J2=J1+1,4; IF(J2.EQ.I1)CYCLE
      IL=IL+1; SLIST(IL)=PRIMN(I1)*PRIMN(J1+4)*PRIMN(J2+4)
      ENDDO; ENDDO; ENDDO
!    +3In
      DO I1=1,4
      DO J1=   1,4; IF(J1.NE.I1)EXIT; ENDDO
      DO J2=J1+1,4; IF(J2.NE.I1)EXIT; ENDDO
      DO J3=J2+1,4; IF(J3.NE.I1)EXIT; ENDDO
      IL=IL+1; SLIST(IL)=PRIMN(J1+4)*PRIMN(J2+4)*PRIMN(J3+4)
      ENDDO

      IF(IL.NE.NMAX1)STOP'CHECK IL COUNT I!'
      CALL SORT_REDIS_ASC(NMAX1,SLIST,IDX)
      CALL CALC_MIN_DEV(SLIST,NMAX1,MINDEV)
      WRITE(IO,'(" MINDEV IN SCORE_LIST OF CYCLE I:",I3)')INT(MINDEV)
      IF(MINDEV.LT.1.E-3_q)STOP'MINDEV-I TOO SMALL!'
      MINDEV=MINDEV*.1_q
      DO NI=1,T_INFO%NIONS
      IF(T_INFO%ITYP(NI).EQ.1.OR.T_INFO%ITYP(NI).EQ.2)CYCLE ! Always 4 Se with fixed direction, simple.
      INN=0; SCORE=1._q
      DO ICX=-1,1; DO ICY=-1,1; DO ICZ=-1,1
      DO NIP=1,T_INFO%NIONS
      DPOS(1)=ICX; DPOS(2)=ICY; DPOS(3)=ICZ
      DPOS=DPOS+T_INFO%POSION(:,NIP)-T_INFO%POSION(:,NI)
      DPOS=MATMUL(LATT_CUR%A,DPOS)
      DIST=SQRT(SUM(DPOS**2))
      IF(LATT_CUR%LROT)DPOS=MATMUL(LATT_CUR%MROT,DPOS)
      IF(DIST.GT.RNN.OR.DIST.LT.1.E-2_q)CYCLE
      INN=INN+1
      IF(INN.GT.NNN)STOP' FETAL ERROR: INN>NNN!'
      NTP=T_INFO%ITYP(NIP)
      SELECT CASE(NTP)
      CASE(1,2) ! Cu/In
        SCORE=SCORE*PRIMN(ID_VEC(DPOS)+(NTP-1)*4)
      CASE DEFAULT ! Do not expect Se here.
        STOP ' UNEXPECTED ATOMIC TYPE!'
      END SELECT
      ENDDO; ENDDO; ENDDO; ENDDO ! NIP,ICZ,ICY,ICX
      SELECT CASE (T_INFO%ITYP(NI))
      CASE(3)
        CALL SEARCH_SCORE(SLIST,NMAX1,SCORE,ID,MINDEV)
        ID=300000+IDX(ID)
      CASE DEFAULT
        STOP 'UNEXPECTED ATOMIC TYPE!'
      END SELECT
      T_INFO%ISORT(NI)=ID
      ENDDO ! NI
      DEALLOCATE(IDX,SLIST)

! SECOND CYCLE
! Score list to see no duplications.
      NMAX2=NMAX1*(NMAX1-1)*(NMAX1-2)*(NMAX1-3)/4/3/2  + & ! C22_4:   1+1+1+1
           &NMAX1*(NMAX1-1)*(NMAX1-2)/3/2*3            + & ! C22_3*3: 2+1+1
           &NMAX1*(NMAX1-1)/2                          + & ! C22_2:   2+2
           &NMAX1*(NMAX1-1)/2*2                        + & ! C22_2*2: 3+1
           &NMAX1                                          ! C22_1:   4
      WRITE(IO,'(" MAXIMAL SORT OF Se In CYCLE II:",I8)')NMAX2
      ALLOCATE(SLIST(NMAX2),IDX(NMAX2))
      DO J1=1,NMAX2; IDX(J1)=J1; ENDDO
      IL=0
! C22_4:   1+1+1+1
      DO J1=1,NMAX1; DO J2=J1+1,NMAX1; DO J3=J2+1,NMAX1; DO J4=J3+1,NMAX1
        IL=IL+1; SLIST(IL)=PRIMN(J1)*PRIMN(J2)*PRIMN(J3)*PRIMN(J4)
      ENDDO; ENDDO; ENDDO; ENDDO
! C22_3*3: 2+1+1
      DO J1=1,NMAX1; DO J2=J1+1,NMAX1; DO J3=J2+1,NMAX1
        IL=IL+1; SLIST(IL)=PRIMN(J1)*PRIMN(J2)*PRIMN(J3)*PRIMN(J1)
        IL=IL+1; SLIST(IL)=PRIMN(J1)*PRIMN(J2)*PRIMN(J3)*PRIMN(J2)
        IL=IL+1; SLIST(IL)=PRIMN(J1)*PRIMN(J2)*PRIMN(J3)*PRIMN(J3)
      ENDDO; ENDDO; ENDDO
! C22_2:   2+2; 3+1
      DO J1=1,NMAX1; DO J2=J1+1,NMAX1
        IL=IL+1; SLIST(IL)=PRIMN(J1)*PRIMN(J2)*PRIMN(J1)*PRIMN(J2)
        IL=IL+1; SLIST(IL)=PRIMN(J1)*PRIMN(J2)*PRIMN(J2)*PRIMN(J2)
        IL=IL+1; SLIST(IL)=PRIMN(J1)*PRIMN(J1)*PRIMN(J1)*PRIMN(J2)
      ENDDO; ENDDO
! C22_1:   4
      DO J1=1,NMAX1; IL=IL+1; SLIST(IL)=PRIMN(J1)**4; ENDDO
      IF(IL.NE.NMAX2)STOP'CHECK IL COUNT II!'
      CALL SORT_REDIS_ASC(NMAX2,SLIST,IDX)
      CALL CALC_MIN_DEV(SLIST,NMAX2,MINDEV)
      WRITE(IO,'(" MINDEV IN SCORE_LIST OF CYCLE II:",I3)')INT(MINDEV)
      IF(MINDEV.LT.1.E-3_q)STOP'MINDEV-II TOO SMALL!'
      MINDEV=MINDEV*.1_q
! Assign sort
      DO NI=1,T_INFO%NIONS
      IF(T_INFO%ITYP(NI).EQ.3)CYCLE ! Since only one sort for Cu/In, No change for Se in this cycle
      INN=0; SCORE=1._q
      DO ICX=-1,1; DO ICY=-1,1; DO ICZ=-1,1
      DO NIP=1,T_INFO%NIONS
      DPOS(1)=ICX; DPOS(2)=ICY; DPOS(3)=ICZ
      DPOS=DPOS+T_INFO%POSION(:,NIP)-T_INFO%POSION(:,NI)
      DPOS=MATMUL(LATT_CUR%A,DPOS)
      DIST=SQRT(SUM(DPOS**2))
      IF(DIST.GT.RNN.OR.DIST.LT.1.E-2_q)CYCLE
      INN=INN+1
      IF(INN.GT.NNN)STOP' FETAL ERROR: INN>NNN!'
      SELECT CASE(T_INFO%ITYP(NIP))
      CASE(3) ! Se
        SCORE=SCORE*PRIMN(MOD(T_INFO%ISORT(NIP),100000))
      CASE DEFAULT ! Do not expect Cu/In here.
        STOP ' UNEXPECTED ATOMIC TYPE!'
      END SELECT
      ENDDO; ENDDO; ENDDO; ENDDO ! NIP,ICZ,ICY,ICX
      NT=T_INFO%ITYP(NI)
      SELECT CASE (NT)
      CASE(1,2) ! Cu
        CALL SEARCH_SCORE(SLIST,NMAX2,SCORE,ID,MINDEV)
        ID=NT*100000+IDX(ID)
      CASE DEFAULT
        STOP 'UNEXPECTED ATOMIC TYPE!'
      END SELECT
      T_INFO%ISORT(NI)=ID
      ENDDO ! NI
! Classify
      DO NI=1,T_INFO%NIONS
      DO IL=1,ENQB%NSORT
        IF(ENYZD(IL,1,1,1)%ID.NE.T_INFO%ISORT(NI))CYCLE
        GOTO 100
      ENDDO
      WRITE(0,'(" ID=",I6)')T_INFO%ISORT(NI)
      STOP 'FETAL ERROR: ID NOT FOUND IN DATABASE!'
100   CONTINUE
      T_INFO%ISORT(NI)=IL
      ENDDO ! NI
      WRITE(IO,'(" NI:")')
      WRITE(IO,'(10I8)')(NI,NI=1,T_INFO%NIONS)
      WRITE(IO,'(" NT:")')
      WRITE(IO,'(10I8)')(T_INFO%ITYP(NI),NI=1,T_INFO%NIONS)
      WRITE(IO,'(" LABEL:")')
      WRITE(IO,'(10I8)')(ENYZD(T_INFO%ISORT(NI),1,1,1)%ID,NI=1,T_INFO%NIONS)
      WRITE(*,'(" DONE.")')
      RETURN

      END SUBROUTINE SORT_ATOM

!***********************************************************************
      FUNCTION ID_VEC(V)
      USE prec
      IMPLICIT NONE
      INTEGER ID_VEC
      REAL(q) V(3)

      IF    (V(1)<0.AND.V(2)>0.AND.V(3)>0)THEN
        ID_VEC=1
      ELSEIF(V(1)>0.AND.V(2)<0.AND.V(3)>0)THEN
        ID_VEC=2
      ELSEIF(V(1)>0.AND.V(2)>0.AND.V(3)<0)THEN
        ID_VEC=3
      ELSEIF(V(1)<0.AND.V(2)<0.AND.V(3)<0)THEN
        ID_VEC=4
      ELSEIF(V(1)>0.AND.V(2)>0.AND.V(3)>0)THEN
        ID_VEC=5
      ELSEIF(V(1)<0.AND.V(2)<0.AND.V(3)>0)THEN
        ID_VEC=6
      ELSEIF(V(1)>0.AND.V(2)<0.AND.V(3)<0)THEN
        ID_VEC=7
      ELSEIF(V(1)<0.AND.V(2)>0.AND.V(3)<0)THEN
        ID_VEC=8
      ELSE
        WRITE(0,'(" PAIR-V:",3F8.3)')V
        STOP 'CHECK PAIR-V!'
      ENDIF
      RETURN

      END FUNCTION ID_VEC

!=======================================================================
      SUBROUTINE SET_QB_YZD(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER NI1,NBASE1,NI2,NBASE2,IC,IR,I1,I2,ICP,IRM,GJ1,GJ2
      REAL(q) DIST,DPOS(3)
      REAL(q),PARAMETER::SMALL=.1_q

      WRITE(*,'(" SET_QB_YZD...")')
      QB%NSP=ENQB%QO%NSP
      ALLOCATE(YZD(T_INFO%NIONS,T_INFO%NIONS))
      DO NI1=1,T_INFO%NIONS; YZD(NI1,1)%DIM=ENYZD(T_INFO%ISORT(NI1),1,1,1)%DIM; ENDDO
      QB%DIM=SUM(YZD(:,1)%DIM)
      WRITE(IO,'(" QB%DIM=",I8)')QB%DIM
      WRITE(* ,'(" QB%DIM=",I8)')QB%DIM
      ALLOCATE(QB%HR(QB%DIM,1,QB%DIM,1,CELL%DIM,QB%NSP),QB%SR(QB%DIM,1,QB%DIM,1,CELL%DIM,QB%NSP))
      QB%HR=0; QB%SR=0
      NBASE1=0
      DO NI1=1,T_INFO%NIONS
      NBASE2=0
      DO NI2=1,T_INFO%NIONS
      YZD(NI1,NI2)%HR=>QB%HR(NBASE1+1:NBASE1+YZD(NI1,1)%DIM,:,NBASE2+1:NBASE2+YZD(NI2,1)%DIM,:,:,:)
      YZD(NI1,NI2)%SR=>QB%SR(NBASE1+1:NBASE1+YZD(NI1,1)%DIM,:,NBASE2+1:NBASE2+YZD(NI2,1)%DIM,:,:,:)
      NBASE2=NBASE2+YZD(NI2,1)%DIM
      ENDDO
      NBASE1=NBASE1+YZD(NI1,1)%DIM
      ENDDO ! NI1
      DO NI1=1,T_INFO%NIONS; I1=T_INFO%ISORT(NI1); GJ1=T_INFO%GJ(NI1)
      DO NI2=1,T_INFO%NIONS; I2=T_INFO%ISORT(NI2); GJ2=T_INFO%GJ(NI2)
      DO IC=1,CELL%DIM
      DPOS=CELL%R(:,IC)+T_INFO%POSION(:,NI2)-T_INFO%POSION(:,NI1)
      DPOS=MATMUL(LATT_CUR%A,DPOS)
      DIST=SQRT(SUM(DPOS**2))
      IF(DIST.GT.ENQB%RCUT)CYCLE
      ICP=0
      DO IR=1,ENQB%NR12(I1,GJ1,I2,GJ2)
      IF(ABS(DIST   -ENQB%R12(0,I1,GJ1,I2,GJ2,IR)).GT.SMALL)CYCLE
      IF(ABS(DPOS(1)-ENQB%R12(1,I1,GJ1,I2,GJ2,IR)).GT.SMALL)CYCLE
      IF(ABS(DPOS(2)-ENQB%R12(2,I1,GJ1,I2,GJ2,IR)).GT.SMALL)CYCLE
      IF(ABS(DPOS(3)-ENQB%R12(3,I1,GJ1,I2,GJ2,IR)).GT.SMALL)CYCLE
      ICP=ICP+1; IRM=IR
      IF(ICP.GT.1)STOP 'MORE THAN ONE MATCH IDENTIFIED! DECREASE PARAMETER-SMALL!'
      ENDDO ! IR
      IF(ICP.EQ.0)THEN
        WRITE(0,'(" NI1,NI2:",2I4)')NI1,NI2
        WRITE(0,'(" DIST,CELL%R:",F8.2,3I3)')DIST,CELL%R(:,IC)
        WRITE(0,'("FAILED TO LOCATE MATCHED PAIR IN THE DATABASE!")'); read(*,*)
      ENDIF
      YZD(NI1,NI2)%HR(:,1,:,1,IC,:)=ENYZD(I1,GJ1,I2,GJ2)%HR(:,1,:,1,IRM,:)
      YZD(NI1,NI2)%SR(:,1,:,1,IC,:)=ENYZD(I1,GJ1,I2,GJ2)%SR(:,1,:,1,IRM,:)
      ENDDO ! IC
      ENDDO; ENDDO ! NI2,NI1
      WRITE(*,'(" DONE.")')
      RETURN

      END SUBROUTINE SET_QB_YZD

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_INCAR(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      INTEGER I1,I2,I3,IR

      WRITE(*,'(" RD_INCAR...")')
      OPEN(IU,FILE='INCAR',STATUS='OLD')
      READ(IU,*)CELL%N(1:3)
      CLOSE(IU)
      WRITE(*,'(" CELL%N:",3I3)')CELL%N
      CELL%N=CELL%N/2*2+1
      CELL%DIM=CELL%N(1)*CELL%N(2)*CELL%N(3)
      WRITE(IO,'(" CELL%N:",3I4)')CELL%N
      ALLOCATE(CELL%R(3,CELL%DIM))
      IR=0
      DO I1=-CELL%N(1)/2,CELL%N(1)/2
      DO I2=-CELL%N(2)/2,CELL%N(2)/2
      DO I3=-CELL%N(3)/2,CELL%N(3)/2
      IR=IR+1
      CELL%R(1,IR)=I1; CELL%R(2,IR)=I2; CELL%R(3,IR)=I3
      ENDDO; ENDDO; ENDDO
      WRITE(*,'(" NONE.")')
      RETURN

      END SUBROUTINE RD_INCAR

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_KPOINTS(KPOINTS,FNAME,IU,IO)
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

      WRITE(*,'(" RD_KPOINTS...")')
      KPOINTS%NKPT=0
      WRITE(IO,'(" READING KPOINTS")')
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
        IF ((CSEL=='G').OR.(CSEL=='g')) THEN
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
        ELSE
          WRITE(IO,'(" CSEL \= G,g: Not supported!")')
          STOP
        ENDIF
      ENDIF
      CLOSE(IU)
      WRITE(IO,'(" KPOINTS%NKPTS=",I8)')KPOINTS%NKPTS
      WRITE(* ,'(" KPOINTS%NKPTS=",I8)')KPOINTS%NKPTS
      WRITE(IO,'(" KPOINTS (50,Direct):")')
      WRITE(IO,'(" NK   K    WT")')
      DO NKP=1,MIN(KPOINTS%NKPTS,5)
        WRITE(IO,'(I4,3F8.3,F10.5)')NKP,KPOINTS%VKPT(:,NKP),KPOINTS%WTKPT(NKP)
      ENDDO
      WRITE(IO,'(" KPOINTS (50,Cartesian):")')
      WRITE(IO,'(" NK   K    WT")')
      DO NKP=1,MIN(KPOINTS%NKPTS,5)
        WRITE(IO,'(I4,3F8.3,F10.5)')NKP,MATMUL(LATT_CUR%B,KPOINTS%VKPT(:,NKP)),KPOINTS%WTKPT(NKP)
      ENDDO
      WRITE(*,'(" DONE.")')
      RETURN

      END SUBROUTINE RD_KPOINTS

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
      SUBROUTINE SET_BAND()
      USE prec
      IMPLICIT NONE

      W%NBANDS=QB%DIM; W%ISPIN=QB%NSP; W%RSPIN=3-W%ISPIN
! WAVE
      ALLOCATE(W%CELEN (W%NBANDS,KPOINTS%NKPTS,W%ISPIN), &
               W%FERWE (W%NBANDS,KPOINTS%NKPTS,W%ISPIN))
      W%CELEN=0; W%FERWE=0
! DOS
      KPOINTS%ISMEAR=0
      KPOINTS%SIGMA=0.1_q
      KPOINTS%NEDOS=1501
      KPOINTS%EMIN= -5._q
      KPOINTS%EMAX=  10._q
      ALLOCATE(KPOINTS%DOS(KPOINTS%NEDOS,W%ISPIN),KPOINTS%DOSI(KPOINTS%NEDOS,W%ISPIN))
      KPOINTS%DOS=0; KPOINTS%DOSI=0
      RETURN

      END SUBROUTINE SET_BAND

!***********************************************************************
      SUBROUTINE SET_HS1K(VKPT,AK,AR)
      USE prec; USE constant
      IMPLICIT NONE
      REAL(q) VKPT(3),AR(QB%DIM,QB%DIM,CELL%DIM)
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
      SUBROUTINE CALC_FULL_BAND(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      INTEGER ISP,K,NB,NI,NBASE
      REAL(q) KPT(3)
      COMPLEX(q),ALLOCATABLE::EVEC(:,:)
      REAL(q),ALLOCATABLE :: PJ(:,:)

      ALLOCATE(EVEC(QB%DIM,QB%DIM)); EVEC=0
      ALLOCATE(PJ(QB%DIM,QB%DIM)); PJ=0
      WRITE(*,'(" CALC_FULL_BAND...")')
      OPEN(IU,FILE='PROCAR',STATUS='REPLACE')
      WRITE(IU,'(I2,I4,I8," ! ISPIN,NKPTS,DIM")')W%ISPIN,KPOINTS%NKPTS,QB%DIM
      DO ISP=1,W%ISPIN
      DO K=1,KPOINTS%NKPTS
        WRITE(*,'(" K=",I8)')K
        KPT=KPOINTS%VKPT(:,K)
        CALL EIGEN_SOL(KPT,-K,W%CELEN(:,K,ISP),EVEC,QB%DIM,'V',ISP,PJ=PJ)
        DO NB=1,QB%DIM
        WRITE(IU,'(" ISP=",I2," K=",I8," NB=",I8," E/WT=",2F12.6)')ISP,K,NB,W%CELEN(NB,K,ISP),KPOINTS%WTKPT(K)
        NBASE=0
        DO NI=1,T_INFO%NIONS
        WRITE(IU,'(I5,16F8.4)')NI,PJ(NBASE+1:NBASE+YZD(NI,1)%DIM,NB); NBASE=NBASE+YZD(NI,1)%DIM
        ENDDO; ENDDO
      ENDDO; ENDDO ! ISP,K
      CLOSE(IU)
      WRITE(*,'(" DONE.")')
101   FORMAT(<T_INFO%NTYP>I8)
      RETURN

      END SUBROUTINE CALC_FULL_BAND

!***********************************************************************
      SUBROUTINE EIGEN_SOL(KPT,K,W,A,N,JOBZ,ISP,PJ)
      USE prec
      IMPLICIT NONE
      INTEGER K,N,ISP
      REAL(q) W(N),KPT(3)
      COMPLEX(q) A(N,N)
      CHARACTER*1 JOBZ
      REAL(q),OPTIONAL,INTENT(OUT)::PJ(N,N)
! Local
      INTEGER LWORK,INFO,IV,IVP
      REAL(q) RWORK(3*N-2)
      COMPLEX(q) ZERO,ZONE
      COMPLEX(q),ALLOCATABLE :: B(:,:),WORK(:),BP(:,:),SV(:)

      CALL  GET_EFF_HAMILT(A,N,KPT,K,ISP)
      LWORK=32*N
      ALLOCATE(WORK(LWORK)); WORK=0
      ALLOCATE(B(N,N)); B=0
      CALL  GET_EFF_OVERLAP(B,N,KPT,K,ISP)
      IF(PRESENT(PJ))THEN
      ALLOCATE(BP(N,N),SV(N)); BP=B; SV=0
      ENDIF

      CALL ZHEGV(1,JOBZ,'L',N,A,N,B,N,W,WORK,LWORK,RWORK,INFO)
      IF(PRESENT(PJ))THEN
      ZERO=0; ZONE=1._q
      DO IV=1,N
        CALL ZGEMM('N','N',N,1,N,ZONE,BP,N,A(:,IV),N,ZERO,SV,N)
        DO IVP=1,N; PJ(IVP,IV)=CONJG(A(IVP,IV))*SV(IVP); ENDDO
      ENDDO
      DEALLOCATE(BP,SV)
      ENDIF ! PRESENT(PJ)
      DEALLOCATE(WORK,B)
      IF(INFO.NE.0)THEN
        WRITE(0,'(" EIGEN_SOL: INFO=",I4)')INFO
        STOP
      ENDIF
      RETURN

      END SUBROUTINE EIGEN_SOL

!***********************************************************************
      SUBROUTINE GET_EFF_HAMILT(HMAT,NDIM,KPT,K,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(q) KPT(3)
      COMPLEX(q) HMAT(NDIM,NDIM)
! Local
      INTEGER N1,N2

      CALL GET_UCTBH(HMAT,NDIM,KPT,K,ISP)
      RETURN

      END SUBROUTINE GET_EFF_HAMILT

!***********************************************************************
      SUBROUTINE GET_UCTBH(HMAT,NDIM,KPT,K,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NDIM,K,ISP
      REAL(q) KPT(3)
      COMPLEX(q) HMAT(NDIM,NDIM)
! LOCAL
      INTEGER ISP_
      REAL(q),POINTER::HR(:,:,:)

      ISP_=MIN(ISP,QB%NSP)
      IF(K.GT.0)THEN
        HMAT=QB%HK(:,:,K,ISP_)
      ELSE
        HR=>QB%HR(:,1,:,1,:,ISP_)
        CALL SET_HS1K(KPT,HMAT,HR)
        NULLIFY(HR)
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
! Local
      INTEGER N1,N2

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
! LOCAL
      INTEGER ISP_
      REAL(q),POINTER::SR(:,:,:)

      ISP_=MIN(ISP,QB%NSP)
      IF(K.GT.0)THEN
        SMAT=QB%SK(:,:,K,ISP_)
      ELSE
        SR=>QB%SR(:,1,:,1,:,ISP_)
        CALL SET_HS1K(KPT,SMAT,SR)
        NULLIFY(SR)
      ENDIF
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

      WRITE(*,'(" SET_FERMI_WT...")')
      IF(KPOINTS%EMIN.GT.KPOINTS%EMAX)THEN
        DO ISP=1,W%ISPIN
        DO K=1,KPOINTS%NKPTS
        DO N=1,W%NBANDS
          E=W%CELEN(N,K,ISP)
          KPOINTS%EMAX=MAX(KPOINTS%EMAX,E)
          KPOINTS%EMIN=MIN(KPOINTS%EMIN,E)
        ENDDO; ENDDO; ENDDO
        EADD=(KPOINTS%EMAX-KPOINTS%EMIN)*0.05_q
        EADD=MAX(EADD,10*ABS(KPOINTS%SIGMA))
        KPOINTS%EMIN=KPOINTS%EMIN-EADD
        KPOINTS%EMAX=KPOINTS%EMAX+EADD
      ENDIF
! Determine fermi level
      CALL DENMP(W%ISPIN,W%RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,W%NET,W%EFERMI, &
             & KPOINTS%ISMEAR,KPOINTS%SIGMA,W%NBANDS,KPOINTS%NKPTS,W%FERWE,W%CELEN, &
             & KPOINTS%NEDOS,KPOINTS%DOS,KPOINTS%DOSI,KPOINTS%WTKPT)
      WRITE(IO,'(" EFERMI=",F8.3)')W%EFERMI
      WRITE(*,'(" DONE.")')
      RETURN

      END SUBROUTINE SET_FERMI_WT

!********************************************************************
! PRINT OUT TDOS
!********************************************************************
      SUBROUTINE OUT_TDOS(IT)
      USE prec
      IMPLICIT NONE
      INTEGER IT,I
      REAL(q) DELTAE

      DELTAE=(KPOINTS%EMAX-KPOINTS%EMIN)/(KPOINTS%NEDOS-1)
      WRITE(71,'("#TDOS: ITER=",I)')IT
      DO I=2,KPOINTS%NEDOS
        WRITE(71,'(3F12.4)')KPOINTS%EMIN+(I-1)*DELTAE-W%EFERMI,KPOINTS%DOS(I,1),KPOINTS%DOSI(I,1)
      ENDDO
      WRITE(71,*)

      END SUBROUTINE OUT_TDOS


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
             & NEDOS,DOS,DOSI,WTKPT)
      USE prec; USE constant
      IMPLICIT NONE
      INTEGER ISPIN,RSPIN,ISMEAR,NBANDS,NKPTS,NEDOS
      REAL(q) DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN)
      REAL(q) NELECT,EMIN,EMAX,EFERMI,SIGMA
      REAL(q) FERWE(NBANDS,NKPTS,ISPIN),CELEN(NBANDS,NKPTS,ISPIN)
      REAL(q) WTKPT(NKPTS)
! Local
      INTEGER ISP,K,N,NELOW,NEHIG,I,NITER
      REAL(q) DELTAE,EPS,WEIGHT,SFUN_DONE,E,EPSDOS,DFUN,SFUN,DOSTOT
      REAL(q) EF1,EF2,ELECT,X1
      LOGICAL LOWB,HIGHB

      DELTAE=(EMAX-EMIN)/(NEDOS-1)
      DOS =0
      DOSI=0
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
! Search a closest number in array A. (with criteria SMALL)
! A in ascending order
!***********************************************************************
      SUBROUTINE SEARCH_SCORE(A,N,X,IX,SMALL)
      USE prec
      IMPLICIT NONE
      INTEGER N,IX
      REAL(q) A(N),X,SMALL
! LOCAL
      INTEGER JL,JM,JU

      IF(ABS(A(1)-X).LT.SMALL)THEN
        IX=1; RETURN
      ELSEIF(ABS(A(N)-X).LT.SMALL)THEN
        IX=N; RETURN
      ENDIF
      IF(A(1).GT.X.OR.A(N).LT.X)STOP'X OUT OF RANGE!'
      JL=0; JU=N+1
      DO
      IF(JU-JL.EQ.0)EXIT
      JM=(JU+JL)/2
      IF(ABS(A(JM)-X).LT.SMALL)THEN
        IX=JM; RETURN
      ELSEIF(A(JM).GT.X)THEN
        IF(JU.EQ.JM)EXIT
        JU=JM
      ELSE
        IF(JL.EQ.JM)EXIT
        JL=JM
      ENDIF
      ENDDO
      STOP 'FAIL TO LOCATE X IN A!'
      RETURN

      END SUBROUTINE SEARCH_SCORE

!***********************************************************************
! A(N): in ascending order
!***********************************************************************
      SUBROUTINE CALC_MIN_DEV(A,N,MINDEV)
      USE prec
      IMPLICIT NONE
      INTEGER N
      REAL(q) A(N),MINDEV
! LOCAL
      INTEGER I

      MINDEV=1.E4_q
      DO I=1,N-1
        MINDEV=MIN(MINDEV,A(I+1)-A(I))
      ENDDO
      RETURN

      END SUBROUTINE CALC_MIN_DEV

!=======================================================================
! sorts RA in ascending order, and rearanges an index array RB
! seems to be a quicksort, by I am not sure (RA,RB Both changed)
! subroutine writen by Florian Kirchhof
!=======================================================================
      SUBROUTINE SORT_REDIS_ASC(N,RA,RB)
      USE prec
      IMPLICIT NONE
      INTEGER N,RB(N)
      REAL(q) RA(N)
! LOCAL
      REAL(q) RRA
      INTEGER RRB,I,L,J,IR

      IF (N==0) RETURN

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      RETURN

      END SUBROUTINE
