!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine calc_nex(nex_c1,nex_c2)
  use global_variables
  implicit none
  real(8),intent(out) :: nex_c1,nex_c2
  real(8) :: nex_c1_l,nex_c2_l
  integer :: ik,ikr,ikz


  nex_c1_l = 0d0
  nex_c2_l = 0d0
  do ik = NKrz_s,NKrz_e
    ikr = ikr_table(ik); ikz = ikz_table(ik)

! state 1
    nex_c1_l = nex_c1_l + abs(zCt(2,1,ik))**2*kr(ikr)
    nex_c2_l = nex_c2_l + abs(zCt(3,1,ik))**2*kr(ikr)

  end do
  nex_c1_l=nex_c1_l*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  nex_c2_l=nex_c2_l*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

  call MPI_ALLREDUCE(nex_c1_l,nex_c1,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(nex_c2_l,nex_c2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  return
end subroutine calc_nex
