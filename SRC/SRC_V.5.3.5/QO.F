#include "symbol.inc"
!***********************************************************************
! QO Module, Yongxin Yao 10-16-2008
! For serial version only.
!***********************************************************************

      MODULE QO_MODULE

      USE prec; USE radial; USE LATTICE
      IMPLICIT NONE

! Structure for each type of atom (Yuan Su) in QO
      TYPE QO_YS
        INTEGER         :: LMAX,LMMAX,LMQD
        CHARACTER(3),POINTER::LMCHAR(:)  ! s,p,d character
        LOGICAL,POINTER :: LIN(:)        ! Tag for L to be used for QO 
        INTEGER,POINTER :: LPS(:)        ! L-character
        REAL(q)         :: PSMAXN
        TYPE (rgrid)    :: R             ! Target orbital radial grid
        REAL(q),POINTER :: RW(:,:)       ! Target orbital r*wave
        REAL(q),POINTER :: QW(:,:)      ! Rec. rep. of target orbital r*wave
        REAL(q),POINTER :: ALPHA(:)      ! Scaling factors target orbitals.
      END TYPE QO_YS

      TYPE WS_VEC ! Wigner_Seitz Vector
        INTEGER DIM,NGRID(3)
        INTEGER,POINTER :: IR(:,:),NDEG(:)
      END TYPE WS_VEC
! Structure for description of QO
      TYPE QO_DES
        LOGICAL           :: LORTHLOC,LORTHALL,LSKIPDAV,LWRTQO,LWSCEL
        INTEGER           :: KEY ! 1: XFQ's optimised virtual method; 2: YYX's
        INTEGER           :: RSPIN,ISPIN,NKPTS,NQBND,NRPLWV,NPROD,NIONS,NBANDS,NTYP
        INTEGER           :: NPSNL
        INTEGER           :: LDIM,LMDIM
        INTEGER           :: NBMIN,NBMAX,NBSUB
        REAL(q)           :: RCUT,DFERMIUP,DFERMIDN,EFERMI,CCUT,RBV,NET ! RBV: dimension of BV space
        REAL(q)           :: EWINDN,EWINUP ! Outer energy window for original quambo construction, relative to fermi level
        INTEGER,POINTER   :: ITYP(:),NITYP(:)
        INTEGER,POINTER   :: HOMO(:,:)
        INTEGER,POINTER   :: IONIQ(:,:)    ! Mapping ION_LM to IQ
        INTEGER,POINTER   :: IQION(:)      ! Mapping IQ TO ION
        CHARACTER(3),POINTER::IQLMCHAR(:)  ! s,p,d character
        REAL(q),POINTER   :: PAIRDIST(:,:,:) ! DIST of <NI|NIP(R)>
        REAL(q),POINTER   :: NINIPR(:,:,:,:) ! VECTOR NIP-->NI
        REAL(q),POINTER   :: CELEN(:,:,:)    ! eigenvalues, safe here for real. for basis.
        REAL(q),POINTER   :: TBH(:,:,:,:),TBS(:,:,:,:)
        COMPLEX(q),POINTER:: WATOM(:,:)    ! LCAO of Atoms
        COMPLEX(q),POINTER:: WAPRJ(:,:)    ! projector (complex) of LCAO
        COMPLEX(q),POINTER:: CPTWFP(:,:) ! chosen wavefunctions as basis.
        COMPLEX(q),POINTER:: CPROJ (:,:) ! projector (complex)  for basis
        COMPLEX(q),POINTER:: PSICOEF(:,:,:,:)! Coefficients of QO for each Psi basis.
        REAL(q),POINTER   :: POSREF(:,:)
        TYPE (LATT)          LAT
        TYPE (WS_VEC)     :: WS
      END TYPE QO_DES

      TYPE (QO_YS ),ALLOCATABLE,SAVE:: QOYS(:)
      TYPE (QO_DES),SAVE:: QO
      INTEGER,PARAMETER:: QOIO=77,QOIU6=76,QOIU5=75,QOIU4=74,QOIU3=73,QONKDIM=500

      END MODULE QO_MODULE
