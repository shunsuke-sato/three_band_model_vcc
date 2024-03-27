!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
! Three band model (3D)
module global_variables
! mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zI = (0d0,1d0)
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0

! Material parameter
  integer :: NKr,NKz,NKrz,NKrz_s,NKrz_e,NKrz_ave,NKrz_remainder
  integer,allocatable :: ikr_table(:),ikz_table(:)
  real(8) :: kr_max,kz_max,dkr,dkz
  real(8),allocatable :: kr(:),kz0(:),kz(:)
  real(8),parameter :: fact_intra = 1d0, fact_inter = 1d0
! d: (semi) core , v: valence, c: conduction
  real(8),parameter :: eps_d = -(0.1440380881432848d+001 - 0.2586486109145065d+00)
  real(8),parameter :: eps_c1 = 0d0
  real(8),parameter :: eps_c2 = eps_c1 + 0.2d0/(2d0*Ry) 
  real(8),parameter :: mass_c1 = 2.5d0, mass_c2 = mass_c1
  real(8),parameter :: piz_dc1 = 1d0,piz_dc2 = 0d0
  real(8),parameter :: piz_dcc = 0d0
  real(8),parameter :: deps12_2 = (0.5d0/(2d0*Ry) )**2

! Time-propagation
  integer :: Nt
  real(8) :: dt
  real(8),allocatable :: eps(:,:)
  complex(8),allocatable :: zCt(:,:,:)
  real(8),allocatable :: Act(:),Act_dt2(:),jtz(:),jtz_intra(:),jtz_inter(:)
  character(20) :: envelope_1,envelope_2

! Laser parameters
  real(8) :: E0_1,omega_1,tpulse_1,omega_ev_1,tpulse_fs_1,Iwcm2_1
  real(8) :: E0_2,omega_2,tpulse_2,omega_ev_2,tpulse_fs_2,Iwcm2_2
  real(8) :: Tdelay_fs,Tdelay

! MPI
  include 'mpif.h'
  integer :: Myrank,Nprocs,ierr

! I/O
  integer,parameter :: nf_current = 21
  character(99),parameter :: cf_current = "current.out"
end module global_variables
