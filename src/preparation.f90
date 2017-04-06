!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation
  use global_variables
  implicit none
  integer :: ikr,ikz,ik


  NKrz = NKr*(2*NKz+1)
  NKrz_ave=NKrz/Nprocs; NKrz_remainder=mod(NKrz,Nprocs)
  if (NKrz_remainder == 0) then
    NKrz_s=NKrz_ave*Myrank+1
    NKrz_e=NKrz_ave*(Myrank+1)
  else
    if (Myrank +1 <= NKrz_remainder) then
      NKrz_s=(NKrz_ave+1)*Myrank+1
      NKrz_e=NKrz_s + (NKrz_ave+1)-1 !(NKrz_ave+1)*(Myrank+1)
    else
      NKrz_s=(NKrz_ave+1)*NKrz_remainder + NKrz_ave*(Myrank-NKrz_remainder)+1
      NKrz_e=NKrz_s + NKrz_ave-1
    end if
  end if

  if(myrank==0)then
    write(*,"(A,2x,I7)")"NKrz=",NKrz
    write(*,"(A,2x,I7)")"NKrz_ave",NKrz_ave
    write(*,"(A,2x,I7)")"NKrz_remainder",NKrz_remainder
  end if
!  write(*,"(9I9)")myrank,NKrz_s,NKrz_e

  allocate(zCt(3,1,NKrz_s:NKrz_e),eps(3,NKrz_s:NKrz_e))
  allocate(kz0(-NKz:NKz),kz(-NKz:NKz),kr(NKr))
  allocate(ikr_table(NKrz),ikz_table(NKrz))
  zCt = 0d0; zCt(1,1,:) = 1d0
  eps = 0d0

! table
  ik = 0
  do ikr = 1,NKr
    do ikz = -NKz,NKz
      ik = ik + 1
      ikr_table(ik) = ikr
      ikz_table(ik) = ikz
    end do
  end do
  

  dkr = kr_max/dble(NKr)
  dkz = kz_max/dble(NKz)

  do ikz = -NKz,NKz
    kz0(ikz) = dkz*dble(ikz)
  end do
  kz = kz0

  do ikr = 1,NKr
    kr(ikr) = dkr*dble(ikr)
  end do

  do ik = NKrz_s,NKrz_e
    ikr = ikr_table(ik)
    ikz = ikz_table(ik)
    kz(ikz) = kz0(ikz)
    eps(1,ik) = eps_d
    eps(2,ik) = eps_c1 +0.5d0/mass_c1*(kr(ikr)**2+kz(ikz)**2)
    eps(3,ik) = eps_c2 +0.5d0/mass_c2*(kr(ikr)**2+kz(ikz)**2)
  end do

end subroutine preparation
