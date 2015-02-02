#include "symbol.inc"
      PROGRAM TIGHTBING
      USE prec; USE TB_MODULE
      IMPLICIT NONE

      OPEN(7,FILE='log',STATUS='REPLACE')
      CALL VERSION()
      CALL INI_TB(14,7)
      CALL CALC_BAND_STRU(14,7)
!      CALL CALC_FS(14,7)
      CLOSE(7)

      END PROGRAM TIGHTBING

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE VERSION()
      USE prec
      IMPLICIT NONE
      CHARACTER DATE*8,TIME*10

      WRITE(7,'("************************************")')
      WRITE(7,'(" Tight-binding code 2009.9 by Y.Yao ")')
      WRITE(7,'(" Mainly following VASP terminology. ")')
      WRITE(7,'("************************************")')
      CALL DATE_AND_TIME(date, time)
      WRITE(7,'(" Current date and time:",A8,":",A10)')DATE,TIME

      END SUBROUTINE VERSION

!***********************************************************************
      SUBROUTINE CALC_BAND_STRU(IU,IO)
      USE prec; USE TB_MODULE
      IMPLICIT NONE
      INTEGER IU,IO

      CALL RD_KPOINTS(KPT_BAND,'KPOINTS_BAND',IU,IO)
      CALL CALC_BAND_STRU_(IU)
      RETURN

      END SUBROUTINE CALC_BAND_STRU

!***********************************************************************
      SUBROUTINE CALC_BAND_STRU_PJ(IU,IO)
      USE prec; USE TB_MODULE
      IMPLICIT NONE
      INTEGER IU,IO

      CALL RD_KPOINTS(KPT_BAND,'KPOINTS_BAND',IU,IO)
      CALL CALC_BAND_STRU_PJ_(IU)
      RETURN

      END SUBROUTINE CALC_BAND_STRU_PJ

!***********************************************************************
      SUBROUTINE CALC_FS(IU,IO)
      USE prec; USE TB_MODULE
      IMPLICIT NONE
      INTEGER IU,IO

      CALL RD_KPOINTS(KPT_FS,'KPOINTS_FS',IU,IO)
      CALL CALC_FS_(IU)
      RETURN

      END SUBROUTINE CALC_FS

!***********************************************************************
      SUBROUTINE CALC_BAND_STRU_(IU)
      USE prec; USE TB_MODULE
      IMPLICIT NONE
      INTEGER IU
! Local
      INTEGER ISP,NK
      REAL(q) EVAL(QB%DIM),DKLEN,KLEN
      COMPLEX(q),ALLOCATABLE::EVEC(:,:)

      WRITE(*,'(" CALC_BAND_STRU...")')
      ALLOCATE(EVEC(QB%DIM,QB%DIM)); EVEC=0
      DO ISP=1,W%ISPIN
      IF(ISP.EQ.1)THEN
        OPEN(IU,FILE='TBBAND.1',STATUS='REPLACE')
      ELSE
        OPEN(IU,FILE='TBBAND.2',STATUS='REPLACE')
      ENDIF
      KLEN=0
      DO NK=1,KPT_BAND%NKPTS
        CALL EIGEN_SOL(KPT_BAND%VKPT(:,NK),-NK,EVAL,EVEC,QB%DIM,'N',ISP)
        IF(NK.EQ.1)THEN
          DKLEN=0
        ELSE
          CALL CALC_DLEN(DKLEN,KPT_BAND%VKPT(:,NK-1),KPT_BAND%VKPT(:,NK),LATT_CUR%B,LATT_CUR%SCALE)
        ENDIF
        KLEN=KLEN+DKLEN
        WRITE(IU,100)KLEN,EVAL-W%EFERMI
      ENDDO
      CLOSE(IU)
      ENDDO ! ISPIN
100   FORMAT(<QB%DIM+1>F10.4)
      WRITE(*,'(" DONE.")')
      RETURN

      END SUBROUTINE CALC_BAND_STRU_

!***********************************************************************
      SUBROUTINE CALC_BAND_STRU_PJ_(IU)
      USE prec; USE TB_MODULE
      IMPLICIT NONE
      INTEGER IU
! Local
      INTEGER ISP,NK,NB
      REAL(q) DKLEN,KLEN
      COMPLEX(q),ALLOCATABLE:: EVEC(:,:)
      REAL(q),ALLOCATABLE :: EVAL(:,:),PJ(:,:)

      ALLOCATE(EVAL(QB%DIM,KPT_BAND%NKPTS),PJ(QB%DIM,QB%DIM)); EVAL=0; PJ=0
      ALLOCATE(EVEC(QB%DIM,QB%DIM)); EVEC=0
      DO ISP=1,W%ISPIN
      IF(ISP.EQ.1)THEN
        OPEN(IU,FILE='TBBAND.1',STATUS='REPLACE')
      ELSE
        OPEN(IU,FILE='TBBAND.2',STATUS='REPLACE')
      ENDIF
      EVAL=0; PJ=0
      DO NK=1,KPT_BAND%NKPTS
        CALL EIGEN_SOL(KPT_BAND%VKPT(:,NK),-NK,EVAL(:,NK),EVEC,QB%DIM,'V',ISP,PJ=PJ)
      ENDDO
      DO NB=1,QB%DIM
      WRITE(IU,'("# NB=",I5)')NB
      KLEN=0
      DO NK=1,KPT_BAND%NKPTS
      IF(NK.EQ.1)THEN
          DKLEN=0
        ELSE
          CALL CALC_DLEN(DKLEN,KPT_BAND%VKPT(:,NK-1),KPT_BAND%VKPT(:,NK),LATT_CUR%B,LATT_CUR%SCALE)
        ENDIF
        KLEN=KLEN+DKLEN
        WRITE(IU,100)KLEN,EVAL(NB,NK)-W%EFERMI,PJ(1,NB),SUM(PJ(2:4,NB)),SUM(PJ(5:9,NB)),SUM(PJ(10:16,NB))
      ENDDO
      WRITE(IU,*)
      ENDDO
      CLOSE(IU)
      ENDDO ! ISPIN
      DEALLOCATE(EVAL,PJ)
100   FORMAT(<QB%DIM+4>F10.4)
      RETURN

      END SUBROUTINE CALC_BAND_STRU_PJ_

!***********************************************************************
      SUBROUTINE CALC_FS_(IU)
      USE prec; USE TB_MODULE
      IMPLICIT NONE
      INTEGER IU
! Local
      INTEGER ISP,NK,NB1,NB2,I
      COMPLEX(q),ALLOCATABLE :: EVEC(:,:)
      REAL(q),ALLOCATABLE :: EVAL(:,:)

      ALLOCATE(EVAL(QB%DIM,KPT_FS%NKPTS)); EVAL=0
      ALLOCATE(EVEC(QB%DIM,QB%DIM)); EVEC=0
      DO ISP=1,W%ISPIN
      IF(ISP.EQ.1)THEN
        OPEN(IU,FILE='fermi.1.bxsf',STATUS='REPLACE')
      ELSE
        OPEN(IU,FILE='fermi.2.bxsf',STATUS='REPLACE')
      ENDIF
      EVAL=0
      DO NK=1,KPT_FS%NKPTS
        CALL EIGEN_SOL(KPT_FS%VKPT(:,NK),-NK,EVAL(:,NK),EVEC,QB%DIM,'N',ISP)
      ENDDO
      WRITE(IU,'("BEGIN_INFO")')
      WRITE(IU,'(" #Launch as: xcrysden --bxsf fermi.bxsf")')
      WRITE(IU,'("  Fermi Energy: ",F10.4)')W%EFERMI
      WRITE(IU,'("END_INFO")')
      WRITE(IU,*)
! Determine band range for plot
      NB1=0; NB2=0
      DO I=1,W%NBANDS; DO NK=1,KPT_FS%NKPTS
      IF(NB1.EQ.0)THEN
        IF(ABS(EVAL(I,NK)-W%EFERMI).LT..2_q)THEN
        NB1=I; EXIT
        ENDIF
      ELSE
        IF(ABS(EVAL(I,NK)-W%EFERMI).LT..2_q)THEN
        NB2=I; EXIT
        ENDIF
      ENDIF
      ENDDO; ENDDO
      WRITE(IU,'("BEGIN_BLOCK_BANDGRID_3D")')
      WRITE(IU,'("  X1")')
      WRITE(IU,'("  BEGIN_BANDGRID_3D_X1")')
      WRITE(IU,'("    ",I14)')NB2-NB1+1
      WRITE(IU,'("    ",3I14)')KPT_FS%NKPT+1
      WRITE(IU,'("    ",3F14.7)')0._q,0._q,0._q
      DO I=1,3
      WRITE(IU,'("    ",3F14.7)')LATT_CUR%B(:,I)
      ENDDO
      DO I=NB1,NB2
      WRITE(IU,'("    BAND:  ",I2)')I
      WRITE(IU,'("    ",5F10.4)')EVAL(I,:)
      ENDDO
      WRITE(IU,'("  END_BANDGRID_3D_1")')
      WRITE(IU,'("END_BLOCK_BANDGRID_3D")')
      CLOSE(IU)
      ENDDO ! ISPIN
      DEALLOCATE(EVAL)
      RETURN

      END SUBROUTINE CALC_FS_

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
