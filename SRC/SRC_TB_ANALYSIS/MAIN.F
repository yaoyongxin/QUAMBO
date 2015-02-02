      PROGRAM TBANALYSIS
      USE prec; USE TB_MODULE
      IMPLICIT NONE

      OPEN(7,FILE='TB_ANALYSIS.log',STATUS='REPLACE')
      CALL VERSION()
      CALL TB_ANALYSIS(14,7)
      CLOSE(7)

      END PROGRAM TBANALYSIS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE VERSION()
      USE prec
      IMPLICIT NONE
      CHARACTER DATE*8,TIME*10

      WRITE(7,'("************************************")')
      WRITE(7,'(" Tight-binding Post Processing code 2009.9 by Y.Yao ")')
      WRITE(7,'(" Mainly following VASP terminology. ")')
      WRITE(7,'("************************************")')
      CALL DATE_AND_TIME(date, time)
      WRITE(7,'(" Current date and time:",A8,":",A10)')DATE,TIME
      RETURN

      END SUBROUTINE VERSION
