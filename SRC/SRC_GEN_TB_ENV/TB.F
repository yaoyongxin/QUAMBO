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
        INTEGER, POINTER :: ILABEL(:) ! Label for each ion 
        INTEGER, POINTER :: LMDIM(:)  ! QO DIM for each ion
        INTEGER, POINTER :: GJ(:)     ! The nationality of each atom
      END TYPE type_info

      TYPE QO_TYPE
        INTEGER NSP  ! ISPIN FOR TBH/S
        INTEGER DIM
        INTEGER,POINTER :: NI (:)
        INTEGER,POINTER :: EQ (:) ! Corresponding Environment-adapted QO Index
        REAL(q),POINTER :: SR (:,:,:,:)
        REAL(q),POINTER :: HR (:,:,:,:)
        REAL(q),POINTER :: R12(:,:,:,:)
      END TYPE QO_TYPE

      TYPE ENQO_TYPE 
        INTEGER NSP 
        INTEGER DIM,NLABEL,NGJ
        REAL(q) RCUT
        INTEGER,POINTER :: LABEL(:) ! Record each label appeared
        INTEGER,POINTER :: GJ(:)    ! Nationality of each label
        INTEGER,POINTER :: IBASE(:) ! STSART IQ FOR EACH LABEL
        INTEGER,POINTER :: LM(:),ILAB(:) ! For each ENQO
        INTEGER,POINTER :: IR(:,:,:,:),ICP(:,:,:,:,:) ! For distinct pairs of LABELS
        REAL(q),POINTER :: SR (:,:,:,:,:,:,:)
        REAL(q),POINTER :: HR (:,:,:,:,:,:,:)
        REAL(q),POINTER :: R12(:,:,:,:,:,:)
      END TYPE ENQO_TYPE

      TYPE CELL_TYPE
        INTEGER DIM,ORIG
        INTEGER,ALLOCATABLE::R(:,:)
      END TYPE CELL_TYPE

! VARIABLES
      TYPE (latt)          ,SAVE :: LATT_CUR
      TYPE (type_info)     ,SAVE :: T_INFO
      TYPE (CELL_TYPE)     ,SAVE :: CELL
      TYPE (QO_TYPE)       ,SAVE :: QB
      TYPE (ENQO_TYPE)     ,SAVE :: ENQB

      CONTAINS

! PUBLIC SUBROUTINES
!***********************************************************************
      SUBROUTINE INI_TB(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO

      CALL RD_POSCAR(IU,IO)
      CALL READ_TBHS(IU)
      CALL LABEL_ATOM(IO)
      CALL SET_QB_MAP_EQ()
      CALL GEN_ENQB(IO)
      CALL OUT_ENQB(IU,IO)
      RETURN

      END SUBROUTINE INI_TB

! PRIVATE SUBROUTINES
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

      WRITE(IO,'(" READING POSCAR")')
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
      ENDDO
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
!     CALL TOPRIM(T_INFO%NIONS,T_INFO%POSION)

      CLOSE(IU)
      RETURN

 400 CONTINUE
      WRITE(0,*)' No initial positions read in'
      STOP

      END SUBROUTINE RD_POSCAR

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
      SUBROUTINE READ_TBHS(IU)
      USE prec
      IMPLICIT NONE
      INTEGER IU
! LOCAL
      INTEGER IC,IQ1,IQ2,ISP,NT1,NI1,NI2
      REAL(q) RES
      CHARACTER(3) STMP

      OPEN(IU,FILE='TBHC.INP',STATUS='OLD')
! DOS INFO
      READ(IU,*); READ(IU,*)
! BV space CELL info
      READ(IU,*); READ(IU,*)CELL%DIM
      ALLOCATE(CELL%R(3,CELL%DIM))
      DO IC=1,CELL%DIM; READ(IU,*)CELL%R(:,IC); ENDDO
! QO info
      READ(IU,*)
      READ(IU,*)QB%DIM,QB%NSP,ENQB%RCUT
      ALLOCATE(QB%NI(QB%DIM))
      READ(IU,*)
      READ(IU,*)QB%NI
      READ(IU,*)(STMP,IQ1=1,QB%DIM)
      ALLOCATE(T_INFO%LMDIM(T_INFO%NTYP)); T_INFO%LMDIM=0
      NI2=0
      DO IQ1=1,QB%DIM; NI1=QB%NI(IQ1); NT1=T_INFO%ITYP(NI1)
      IF(NI1.NE.NI2)THEN
        IF(T_INFO%LMDIM(NT1).NE.0)THEN; CYCLE
        ELSE; NI2=NI1
        ENDIF
      ENDIF
      T_INFO%LMDIM(NT1)=T_INFO%LMDIM(NT1)+1
      ENDDO
! TBH
      READ(IU,*)
      ALLOCATE(QB%HR(QB%DIM,QB%DIM,CELL%DIM,QB%NSP)); QB%HR=0
      DO
        READ(IU,*,ERR=102,END=102)ISP,IQ1,IQ2,IC,RES
        IF(QB%HR(IQ1,IQ2,IC,ISP).NE.0)THEN
        WRITE(0,'(" ISP,IQ1,IQ2,IC",3I3,I6," ENTRY REPEATED IN TBHC.INP!")')ISP,IQ1,IQ2,IC
        STOP
        ENDIF
        QB%HR(IQ1,IQ2,IC,ISP)=RES
      ENDDO
102   CLOSE(IU)
      OPEN(IU,FILE='TBSC.INP',STATUS='OLD')
      ALLOCATE(QB%SR(QB%DIM,QB%DIM,CELL%DIM,QB%NSP)); QB%SR=0
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
      DO IC=1,CELL%DIM
        IF(CELL%R(1,IC).EQ.0.AND.CELL%R(2,IC).EQ.0.AND.CELL%R(3,IC).EQ.0)THEN
          CELL%ORIG=IC
          GOTO 100
        ENDIF
      ENDDO
      STOP 'NO ORIGIN IN CELL!'
100   CONTINUE
      CALL SET_QB_R12(IU)
      RETURN

      END SUBROUTINE READ_TBHS

!***********************************************************************
      SUBROUTINE SET_QB_R12(IU)
      USE prec
      IMPLICIT NONE
      INTEGER IU
! Local
      INTEGER IQ1,IQ2,NI1,NI2,IC,I
      REAL(q) DPOS(3),MROT(3,3)

      ALLOCATE(QB%R12(0:3,QB%DIM,QB%DIM,CELL%DIM)); QB%R12=0
      DO IQ1=1,QB%DIM; NI1=QB%NI(IQ1); DO IQ2=1,QB%DIM; NI2=QB%NI(IQ2)
      DO IC =1,CELL%DIM
        DPOS=T_INFO%POSION(:,NI2)-T_INFO%POSION(:,NI1)+CELL%R(:,IC)
        DPOS=MATMUL(LATT_CUR%A,DPOS)
        QB%R12(1:3,IQ1,IQ2,IC)=DPOS
        QB%R12(0  ,IQ1,IQ2,IC)=SQRT(SUM(DPOS**2))
      ENDDO; ENDDO; ENDDO
! Possible rotation
      INQUIRE(FILE='ROTCAR',EXIST=LATT_CUR%LROT)
      IF(LATT_CUR%LROT)THEN
      WRITE(*,'(" ADDITIONAL COOR ROTATION (ROTCAR) WOULD BE TAKEN INTO CONSIDERATION!")')
      OPEN(IU,FILE='ROTCAR',STATUS='OLD')
      READ(IU,*)
      READ(IU,*)(LATT_CUR%MROT(I,:),I=1,3)
      CLOSE(IU)
      ENDIF
      RETURN

      END SUBROUTINE SET_QB_R12

!***********************************************************************
      SUBROUTINE LABEL_ATOM(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL 
      INTEGER,PARAMETER::NNN=4     ! MAX number of N.N.s
      REAL(q),PARAMETER::RNN=2.7_q ! N.N. Dist
      INTEGER NI,NIP,NT,NTP,ICX,ICY,ICZ,INN,IL,ID,NMAX1,NMAX2,I1,I2,J1,J2,J3,J4
      REAL(q) DPOS(3),SCORE,DIST,MINDEV
      INTEGER,ALLOCATABLE::IDX(:)
      REAL(q),ALLOCATABLE::SLIST(:)
      REAL(q) LIST1(4,2)
      INTEGER LABEL(200)
      REAL(q),PARAMETER::PRIMN(46)=(/2,   3,  5,    7,  11,  13,  17,  19,  23,  29, &
                                  & 31,  37,  41,  43,  47,  53,  59,  61,  67,  71, &
                                  & 73,  79,  83,  89,  97, 101, 103, 107, 109, 113, &
                                  &127, 131, 137, 139, 149, 151, 157, 163, 167, 173, &
                                  &179, 181, 191, 193, 197, 199/)
      
      ALLOCATE(T_INFO%ILABEL(T_INFO%NIONS))
      DO I1=1,2; DO J1=1,4; LIST1(J1,I1)=PRIMN(J1+(I1-1)*4); ENDDO; ENDDO

! 1ST CYCLE
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
      IL=IL+1; SLIST(IL)=LIST1(I1,1)*LIST1(I2,1)*LIST1(J1,2)*LIST1(J2,2)
      ENDDO; ENDDO
! 1Cu+3In
      DO I1=1,4
      DO J1=   1,4; IF(J1.NE.I1)EXIT; ENDDO
      DO J2=J1+1,4; IF(J2.NE.I1)EXIT; ENDDO
      DO J3=J2+1,4; IF(J3.NE.I1)EXIT; ENDDO
      IL=IL+1; SLIST(IL)=LIST1(I1,1)*LIST1(J1,2)*LIST1(J2,2)*LIST1(J3,2)
      ENDDO
! 1Cu+2In
      DO I1=1,4
      DO J1=   1,4; IF(J1.EQ.I1)CYCLE
      DO J2=J1+1,4; IF(J2.EQ.I1)CYCLE
      IL=IL+1; SLIST(IL)=LIST1(I1,1)*LIST1(J1,2)*LIST1(J2,2)
      ENDDO; ENDDO; ENDDO
!    +3In
      DO I1=1,4
      DO J1=   1,4; IF(J1.NE.I1)EXIT; ENDDO
      DO J2=J1+1,4; IF(J2.NE.I1)EXIT; ENDDO
      DO J3=J2+1,4; IF(J3.NE.I1)EXIT; ENDDO
      IL=IL+1; SLIST(IL)=LIST1(J1,2)*LIST1(J2,2)*LIST1(J3,2)
      ENDDO

      IF(IL.NE.NMAX1)STOP'CHECK IL COUNT I!'
      CALL SORT_REDIS_ASC(NMAX1,SLIST,IDX)
      CALL CALC_MIN_DEV(SLIST,NMAX1,MINDEV)
      WRITE(IO,'(" MINDEV IN SCORE_LIST OF CYCLE I:",I3)')INT(MINDEV)
      IF(MINDEV.LT.1.E-3_q)STOP'MINDEV-I TOO SMALL!'
      MINDEV=MINDEV*.1_q
! Identify sorts
      DO NI=1,T_INFO%NIONS
      IF(T_INFO%ITYP(NI).EQ.1.OR.T_INFO%ITYP(NI).EQ.2)CYCLE ! Always 4 Se with fixed direction, simple.
      INN=0; SCORE=1._q
      DO ICX=-1,1; DO ICY=-1,1; DO ICZ=-1,1
      DO NIP=1,T_INFO%NIONS
      DPOS(1)=ICX; DPOS(2)=ICY; DPOS(3)=ICZ
      DPOS=DPOS+T_INFO%POSION(:,NIP)-T_INFO%POSION(:,NI)
      DPOS=MATMUL(LATT_CUR%A,DPOS)
      IF(LATT_CUR%LROT)DPOS=MATMUL(LATT_CUR%MROT,DPOS)
      DIST=SQRT(SUM(DPOS**2))
      IF(DIST.GT.RNN.OR.DIST.LT.1.E-2_q)CYCLE
      INN=INN+1
      IF(INN.GT.NNN)STOP' FETAL ERROR: INN>NNN!'
      NTP=T_INFO%ITYP(NIP)
      SELECT CASE(NTP)
      CASE(1,2) ! Cu/In
        SCORE=SCORE*LIST1(ID_VEC(DPOS),NTP)
      CASE DEFAULT ! Do not expect Se here.
        STOP ' UNIDENTIFIED ATOMIC TYPE!'
      END SELECT
      ENDDO; ENDDO; ENDDO; ENDDO ! NIP,ICZ,ICY,ICX
      SELECT CASE (T_INFO%ITYP(NI))
      CASE(3)
        CALL SEARCH_SCORE(SLIST,NMAX1,SCORE,ID,MINDEV)
        ID=3E5+IDX(ID)
      CASE DEFAULT
        STOP 'UNEXPECTED ATOMIC TYPE!'
      END SELECT
      T_INFO%ILABEL(NI)=ID
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
        SCORE=SCORE*PRIMN(MOD(T_INFO%ILABEL(NIP),100000))
      CASE DEFAULT ! Do not expect Cu/In here.
        STOP ' UNIDENTIFIED ATOMIC TYPE!'
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
      T_INFO%ILABEL(NI)=ID
      ENDDO ! NI
! Classify
      DO NI=1,T_INFO%NIONS
      IF(NI.EQ.1)THEN
        ENQB%NLABEL=1
        LABEL(ENQB%NLABEL)=T_INFO%ILABEL(NI)
        T_INFO%ILABEL(NI)=ENQB%NLABEL
      ELSE
        DO IL=1,ENQB%NLABEL
        IF(LABEL(IL).NE.T_INFO%ILABEL(NI))CYCLE
        T_INFO%ILABEL(NI)=IL; GOTO 100
        ENDDO
        ENQB%NLABEL=ENQB%NLABEL+1
        LABEL(ENQB%NLABEL)=T_INFO%ILABEL(NI)
        T_INFO%ILABEL(NI)=ENQB%NLABEL
100     CONTINUE
      ENDIF
      ENDDO ! NI
      ALLOCATE(ENQB%LABEL(ENQB%NLABEL)); ENQB%LABEL=LABEL(1:ENQB%NLABEL)
      WRITE(IO,'(" NI:")')
      WRITE(IO,'(10I8)')(NI,NI=1,T_INFO%NIONS)
      WRITE(IO,'(" NT:")')
      WRITE(IO,'(10I8)')(T_INFO%ITYP(NI),NI=1,T_INFO%NIONS)
      WRITE(IO,'(" LABEL:")')
      WRITE(IO,'(10I8)')(ENQB%LABEL(T_INFO%ILABEL(NI)),NI=1,T_INFO%NIONS)
      RETURN

      END SUBROUTINE LABEL_ATOM

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

!***********************************************************************
      SUBROUTINE SET_QB_MAP_EQ()
      USE prec
      IMPLICIT NONE
! LOCAL
      INTEGER ILAB,LMDIM,NBASE,LM,NI,IBASE
      
      ENQB%DIM=0
      DO ILAB=1,ENQB%NLABEL; ENQB%DIM=ENQB%DIM+T_INFO%LMDIM(INT(ENQB%LABEL(ILAB)/100000)); ENDDO
      ENQB%NGJ=MAXVAL(T_INFO%GJ)
      ALLOCATE(ENQB%LM(ENQB%DIM),ENQB%ILAB(ENQB%DIM))
      ALLOCATE(ENQB%IBASE(ENQB%NLABEL))
      NBASE=0
      DO ILAB=1,ENQB%NLABEL
      ENQB%IBASE(ILAB)=NBASE
      LMDIM=T_INFO%LMDIM(INT(ENQB%LABEL(ILAB)/100000))
      DO LM=1,LMDIM; ENQB%LM(NBASE+LM)=LM; ENQB%ILAB(NBASE+LM)=ILAB; ENDDO
      NBASE=NBASE+LMDIM
      ENDDO

      ALLOCATE(QB%EQ(QB%DIM))
      NBASE=0
      DO NI=1,T_INFO%NIONS
      LMDIM=T_INFO%LMDIM(T_INFO%ITYP(NI))
      IBASE=ENQB%IBASE(T_INFO%ILABEL(NI))
      DO LM=1,LMDIM; QB%EQ(NBASE+LM)=IBASE+LM; ENDDO
      NBASE=NBASE+LMDIM
      ENDDO
      IF(NBASE.NE.QB%DIM)STOP ' INCONSISTENCE ERROR: NBASE\=QB%DIM!'
      RETURN

      END SUBROUTINE SET_QB_MAP_EQ

!***********************************************************************
      SUBROUTINE GEN_ENQB(IO)
      USE prec
      IMPLICIT NONE
      INTEGER IO
! LOCAL
      INTEGER ISP,IQ1,IQ2,IC,IEQ1,IEQ2,I,J,GJ1,GJ2
      INTEGER,PARAMETER::NRMAX=100,NCPMAX=100
      INTEGER,POINTER::IR(:,:,:,:),ICP(:,:,:,:,:)
      REAL(q) DEVMAX,DEV
      REAL(q),PARAMETER::SMALL=.1_q

      ENQB%NSP=QB%NSP
      ALLOCATE(ENQB%HR(ENQB%DIM,ENQB%NGJ,ENQB%DIM,ENQB%NGJ,NRMAX,NCPMAX,ENQB%NSP))
      ALLOCATE(ENQB%SR(ENQB%DIM,ENQB%NGJ,ENQB%DIM,ENQB%NGJ,NRMAX,NCPMAX,ENQB%NSP))
      ALLOCATE(ENQB%R12(0:3,ENQB%DIM,ENQB%NGJ,ENQB%DIM,ENQB%NGJ,NRMAX))
      ALLOCATE(ENQB%IR(ENQB%DIM,ENQB%NGJ,ENQB%DIM,ENQB%NGJ),ENQB%ICP(ENQB%DIM,ENQB%NGJ,ENQB%DIM,ENQB%NGJ,NRMAX))
      ENQB%IR=0; ENQB%ICP=0
      IR=>ENQB%IR; ICP=>ENQB%ICP
      DEVMAX=0
      DO ISP=1,ENQB%NSP
      IR=0; ICP=0
      DO IQ1=1,QB%DIM; IEQ1=QB%EQ(IQ1); GJ1=T_INFO%GJ(QB%NI(IQ1))
      DO IQ2=1,QB%DIM; IEQ2=QB%EQ(IQ2); GJ2=T_INFO%GJ(QB%NI(IQ2))
      DO IC=1,CELL%DIM
      IF(QB%R12(0,IQ1,IQ2,IC).GT.ENQB%RCUT)CYCLE
      IF(IR(IEQ1,GJ1,IEQ2,GJ2).GT.0)THEN
        DO I=1,IR(IEQ1,GJ1,IEQ2,GJ2)
        IF(ABS(ENQB%R12(0,IEQ1,GJ1,IEQ2,GJ2,I)-QB%R12(0,IQ1,IQ2,IC)).GT.SMALL)CYCLE
        IF(ABS(ENQB%R12(1,IEQ1,GJ1,IEQ2,GJ2,I)-QB%R12(1,IQ1,IQ2,IC)).GT.SMALL)CYCLE
        IF(ABS(ENQB%R12(2,IEQ1,GJ1,IEQ2,GJ2,I)-QB%R12(2,IQ1,IQ2,IC)).GT.SMALL)CYCLE
        IF(ABS(ENQB%R12(3,IEQ1,GJ1,IEQ2,GJ2,I)-QB%R12(3,IQ1,IQ2,IC)).GT.SMALL)CYCLE
        GOTO 101
        ENDDO ! I
        IR(IEQ1,GJ1,IEQ2,GJ2)=IR(IEQ1,GJ1,IEQ2,GJ2)+1
        ICP(IEQ1,GJ1,IEQ2,GJ2,IR(IEQ1,GJ1,IEQ2,GJ2))=1
        ENQB%R12(0:3,IEQ1,GJ1,IEQ2,GJ2,IR(IEQ1,GJ1,IEQ2,GJ2))=QB%R12(0:3,IQ1,IQ2,IC)
        ENQB%HR(IEQ1,GJ1,IEQ2,GJ2,IR(IEQ1,GJ1,IEQ2,GJ2),1,ISP)=QB%HR(IQ1,IQ2,IC,ISP)
        ENQB%SR(IEQ1,GJ1,IEQ2,GJ2,IR(IEQ1,GJ1,IEQ2,GJ2),1,ISP)=QB%SR(IQ1,IQ2,IC,ISP)
        CYCLE
101     CONTINUE
! Check error of HR
        DEV=0
        DO J=1,ICP(IEQ1,GJ1,IEQ2,GJ2,I)
          DEV=MAX(DEV,ABS(ENQB%HR(IEQ1,GJ1,IEQ2,GJ2,I,J,ISP)-QB%HR(IQ1,IQ2,IC,ISP)))
        ENDDO ! J
        DEVMAX=MAX(DEVMAX,DEV)
        IF(DEV.GT..05_q)THEN
          WRITE(IO,'(" WARNING: SIGNIFISCANT DEVIATION FOUND! MAY CONSIDER NARROW DOWN ENV.")')
          WRITE(IO,'(" HR:")')
          WRITE(IO,'(10F8.3)')ENQB%HR(IEQ1,GJ1,IEQ2,GJ2,I,1:ICP(IEQ1,GJ1,IEQ2,GJ2,I),ISP),QB%HR(IQ1,IQ2,IC,ISP)
          WRITE(IO,'(" SR:")')
          WRITE(IO,'(10F8.3)')ENQB%SR(IEQ1,GJ1,IEQ2,GJ2,I,1:ICP(IEQ1,GJ1,IEQ2,GJ2,I),ISP),QB%SR(IQ1,IQ2,IC,ISP)
          WRITE(IO,'(" R12:")')
          WRITE(IO,'( 4F8.3)')ENQB%R12(0:3,IEQ1,GJ1,IEQ2,GJ2,I)
          WRITE(IO,'(" LM:",2I5)')ENQB%LM(IEQ1),ENQB%LM(IEQ2)
          WRITE(IO,'(" LABEL:",2I8)')ENQB%LABEL(ENQB%ILAB(IEQ1)),ENQB%LABEL(ENQB%ILAB(IEQ2))
          WRITE(IO,'(" GJ:",2I3," NI:",2I6)')GJ1,GJ2,QB%NI(IQ1),QB%NI(IQ2)
          WRITE(IO,'(" ISP,IQ1,IQ2,IC",4I8)')ISP,IQ1,IQ2,IC
        ENDIF
        ICP(IEQ1,GJ1,IEQ2,GJ2,I)=ICP(IEQ1,GJ1,IEQ2,GJ2,I)+1
!        ENQB%R12(0:3,IEQ1,GJ1,IEQ2,GJ2,I)=(ENQB%R12(0:3,IEQ1,GJ1,IEQ2,GJ2,I)*(ICP(IEQ1,IEQ2,I)-1)+QB%R12(0:3,IQ1,IQ2,IC))/ICP(IEQ1,IEQ2,I)
! Not necessary since I use a common rigid lattice to identify them.
        ENQB%HR(IEQ1,GJ1,IEQ2,GJ2,I,ICP(IEQ1,GJ1,IEQ2,GJ2,I),ISP)=QB%HR(IQ1,IQ2,IC,ISP)
        ENQB%SR(IEQ1,GJ1,IEQ2,GJ2,I,ICP(IEQ1,GJ1,IEQ2,GJ2,I),ISP)=QB%SR(IQ1,IQ2,IC,ISP)
      ELSE
        IR (IEQ1,GJ1,IEQ2,GJ2)  =1
        ICP(IEQ1,GJ1,IEQ2,GJ2,1)=1
        ENQB%R12(0:3,IEQ1,GJ1,IEQ2,GJ2,1)=QB%R12(0:3,IQ1,IQ2,IC)
        ENQB%HR(IEQ1,GJ1,IEQ2,GJ2,1,1,ISP)=QB%HR(IQ1,IQ2,IC,ISP)
        ENQB%SR(IEQ1,GJ1,IEQ2,GJ2,1,1,ISP)=QB%SR(IQ1,IQ2,IC,ISP)
      ENDIF
      ENDDO; ENDDO; ENDDO 
      ENDDO
      WRITE(IO,'(" MAX HR DEVIATION:",F8.3)')DEVMAX
      WRITE(IO,'(" MAX ICP:",I4)')MAXVAL(ICP)
      RETURN

      END SUBROUTINE GEN_ENQB

!***********************************************************************
      SUBROUTINE OUT_ENQB(IU,IO)
      USE prec
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL 
      INTEGER IQ1,IQ2,IR,ISP,ILAB,ILABP,NR12,GJ1,GJ2
      REAL(q) MEAN,DEVMAX

      WRITE(* ,'(" TOTAL NUMBER OF SORTS:",I3)')ENQB%NLABEL
      WRITE(IO,'(" TOTAL NUMBER OF SORTS:",I3)')ENQB%NLABEL
      WRITE(* ,'(" DIM_GJ:",I3)')ENQB%NGJ
      WRITE(IO,'(" DIM_GJ:",I3)')ENQB%NGJ
      OPEN(IU,FILE='ENTB.OUT',STATUS='REPLACE')
      WRITE(IU,'(I8," ! NSP")')ENQB%NSP
      WRITE(IU,'(I8," ! NGJ")')ENQB%NGJ
      WRITE(IU,'(I8," ! NSORT")')ENQB%NLABEL
      DO ILAB=1,ENQB%NLABEL
        WRITE(IU,'(I8,I5," ! ID,DIM,ISORT=",I3)')ENQB%LABEL(ILAB),T_INFO%LMDIM(ENQB%LABEL(ILAB)/100000),ILAB
      ENDDO
      NR12=0
      DO GJ1=1,ENQB%NGJ; DO ILAB =1,ENQB%NLABEL; IQ1=ENQB%IBASE(ILAB )+1
      DO GJ2=1,ENQB%NGJ; DO ILABP=1,ENQB%NLABEL; IQ2=ENQB%IBASE(ILABP)+1
      WRITE(IU,'(I5," ! NR12,ISORT/GJ PAIR=",4I4)')ENQB%IR(IQ1,GJ1,IQ2,GJ2),ILAB,GJ1,ILABP,GJ2
      NR12=NR12+ENQB%IR(IQ1,GJ1,IQ2,GJ2)
      ENDDO; ENDDO; ENDDO; ENDDO
      WRITE(* ,'(" TOTAL NUMBER OF P_VECTORS:",I8)')NR12
      WRITE(IO,'(" TOTAL NUMBER OF P_VECTORS:",I8)')NR12
      DO GJ1=1,ENQB%NGJ; DO ILAB =1,ENQB%NLABEL; IQ1=ENQB%IBASE(ILAB )+1
      DO GJ2=1,ENQB%NGJ; DO ILABP=1,ENQB%NLABEL; IQ2=ENQB%IBASE(ILABP)+1
      IF(ENQB%IR(IQ1,GJ1,IQ2,GJ2).EQ.0)CYCLE
      WRITE(IU,'(" R12 LIST: ISORT/GJ PAIR=",4I4)')ILAB,GJ1,ILABP,GJ2
      WRITE(IU,'(4F10.4)')ENQB%R12(0:3,IQ1,GJ1,IQ2,GJ2,1:ENQB%IR(IQ1,GJ1,IQ2,GJ2))
      ENDDO; ENDDO; ENDDO; ENDDO
      DO ISP=1,ENQB%NSP
      DO GJ1=1,ENQB%NGJ; DO IQ1=1,ENQB%DIM
      DO GJ2=1,ENQB%NGJ; DO IQ2=1,ENQB%DIM
      IF(ENQB%IR(IQ1,GJ1,IQ2,GJ2).EQ.0)CYCLE
      WRITE(IU,'(" HR/DEVMAX/ERR%/DIST LIST: IQ/GJ PAIR=",4I4," ISP=",I2)')IQ1,GJ1,IQ2,GJ2,ISP
      DO IR=1,ENQB%IR(IQ1,GJ1,IQ2,GJ2)
        MEAN=SUM(ENQB%HR(IQ1,GJ1,IQ2,GJ2,IR,1:ENQB%ICP(IQ1,GJ1,IQ2,GJ2,IR),ISP))/ENQB%ICP(IQ1,GJ1,IQ2,GJ2,IR)
        DEVMAX=MAXVAL(ABS(ENQB%HR(IQ1,GJ1,IQ2,GJ2,IR,1:ENQB%ICP(IQ1,GJ1,IQ2,GJ2,IR),ISP)-MEAN))
        WRITE(IU,'(2F12.6,I8,"%",F8.2)')MEAN,DEVMAX,INT(DEVMAX/MAX(ABS(MEAN),1.E-8_q)*100),ENQB%R12(0,IQ1,GJ1,IQ2,GJ2,IR)
      ENDDO; ENDDO; ENDDO ; ENDDO; ENDDO
      ENDDO
      DO ISP=1,ENQB%NSP
      DO GJ1=1,ENQB%NGJ; DO IQ1=1,ENQB%DIM
      DO GJ2=1,ENQB%NGJ; DO IQ2=1,ENQB%DIM
      IF(ENQB%IR(IQ1,GJ1,IQ2,GJ2).EQ.0)CYCLE
      WRITE(IU,'(" SR/DEVMAX/ERR%/DIST LIST: IQ/GJ PAIR=",4I4," ISP=",I2)')IQ1,GJ1,IQ2,GJ2,ISP
      DO IR=1,ENQB%IR(IQ1,GJ1,IQ2,GJ2)
        MEAN=SUM(ENQB%SR(IQ1,GJ1,IQ2,GJ2,IR,1:ENQB%ICP(IQ1,GJ1,IQ2,GJ2,IR),ISP))/ENQB%ICP(IQ1,GJ1,IQ2,GJ2,IR)
        DEVMAX=MAXVAL(ABS(ENQB%SR(IQ1,GJ1,IQ2,GJ2,IR,1:ENQB%ICP(IQ1,GJ1,IQ2,GJ2,IR),ISP)-MEAN))
        WRITE(IU,'(2F12.6,I8,"%",F8.2)')MEAN,DEVMAX,INT(DEVMAX/MAX(ABS(MEAN),1.E-8_q)*100),ENQB%R12(0,IQ1,GJ1,IQ2,GJ2,IR)
      ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
      ENDDO
      CLOSE(IU)
      RETURN

      END SUBROUTINE OUT_ENQB


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
