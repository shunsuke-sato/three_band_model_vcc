!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dt_evolve(it) ! Now coding
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ik,ikr,ikz
  complex(8) :: zHmat(3,3),zEig(3,3)
  complex(8) :: zc1,zc2,zc3
  real(8) :: eps_t(3),de12
  real(8) :: Act_old,Act_new,Et_old,Et_new

  Act_old = Act(it); Act_new = Act(it+1)
  Et_old = -0.5d0*(Act(it+1)-Act(it-1))/dt
  Et_new = -0.5d0*(Act(it+1+1)-Act(it+1-1))/dt

  do ik = NKrz_s,NKrz_e
    ikr = ikr_table(ik); ikz = ikz_table(ik)

!== Enforced time-reversal symmetry scheme

!== First half dt
    kz(ikz) = kz0(ikz) + Act_old
    eps_t(1) = eps_d
    eps_t(2) = -0.5d0/mass_v*(kr(ikr)**2+kz(ikz)**2)
    eps_t(3) = eps_g +0.5d0/mass_c*(kr(ikr)**2+kz(ikz)**2)
    de12 = eps_t(1) - eps_t(2)  
    zHmat(1,1) = eps_t(1); zHmat(2,2) = eps_t(2); zHmat(3,3) = eps_t(3)
    zHmat(1,2) = -zI*piz_dv*Et_old*de12/(de12**2+deps12_2); zHmat(2,1)=conjg(zHmat(1,2))
    zHmat(1,3) = -zI*piz_dc*Et_old/(eps_t(1)-eps_t(3)); zHmat(3,1)=conjg(zHmat(1,3))

    zHmat(2,3) = -zI*piz_vc*Et_old/(eps_t(2)-eps_t(3)); zHmat(3,2)=conjg(zHmat(2,3))

!    call diag3x3(zHmat,zEig,eps_t)
    call zheevh3(zHmat,zEig,eps_t)

! state 1
    zc1=sum(conjg(zEig(:,1))*zCt(:,1,ik))*exp(-0.5d0*zI*dt*eps_t(1))
    zc2=sum(conjg(zEig(:,2))*zCt(:,1,ik))*exp(-0.5d0*zI*dt*eps_t(2))
    zc3=sum(conjg(zEig(:,3))*zCt(:,1,ik))*exp(-0.5d0*zI*dt*eps_t(3))
    zCt(:,1,ik)=zc1*zEig(:,1)+zc2*zEig(:,2)+zc3*zEig(:,3)
! state 2
    zc1=sum(conjg(zEig(:,1))*zCt(:,2,ik))*exp(-0.5d0*zI*dt*eps_t(1))
    zc2=sum(conjg(zEig(:,2))*zCt(:,2,ik))*exp(-0.5d0*zI*dt*eps_t(2))
    zc3=sum(conjg(zEig(:,3))*zCt(:,2,ik))*exp(-0.5d0*zI*dt*eps_t(3))
    zCt(:,2,ik)=zc1*zEig(:,1)+zc2*zEig(:,2)+zc3*zEig(:,3)

!== Second half dt
    kz(ikz) = kz0(ikz) + Act_new
    eps_t(1) = eps_d
    eps_t(2) = -0.5d0/mass_v*(kr(ikr)**2+kz(ikz)**2)
    eps_t(3) = eps_g +0.5d0/mass_c*(kr(ikr)**2+kz(ikz)**2)
    de12 = eps_t(1) - eps_t(2)  

    zHmat(1,1) = eps_t(1); zHmat(2,2) = eps_t(2); zHmat(3,3) = eps_t(3)
    zHmat(1,2) = -zI*piz_dv*Et_new*de12/(de12**2+deps12_2); zHmat(2,1)=conjg(zHmat(1,2))
    zHmat(1,3) = -zI*piz_dc*Et_new/(eps_t(1)-eps_t(3)); zHmat(3,1)=conjg(zHmat(1,3))
    zHmat(2,3) = -zI*piz_vc*Et_new/(eps_t(2)-eps_t(3)); zHmat(3,2)=conjg(zHmat(2,3))

!    call diag3x3(zHmat,zEig,eps_t)
    call zheevh3(zHmat,zEig,eps_t)

! state 1
    zc1=sum(conjg(zEig(:,1))*zCt(:,1,ik))*exp(-0.5d0*zI*dt*eps_t(1))
    zc2=sum(conjg(zEig(:,2))*zCt(:,1,ik))*exp(-0.5d0*zI*dt*eps_t(2))
    zc3=sum(conjg(zEig(:,3))*zCt(:,1,ik))*exp(-0.5d0*zI*dt*eps_t(3))
    zCt(:,1,ik)=zc1*zEig(:,1)+zc2*zEig(:,2)+zc3*zEig(:,3)
! state 2
    zc1=sum(conjg(zEig(:,1))*zCt(:,2,ik))*exp(-0.5d0*zI*dt*eps_t(1))
    zc2=sum(conjg(zEig(:,2))*zCt(:,2,ik))*exp(-0.5d0*zI*dt*eps_t(2))
    zc3=sum(conjg(zEig(:,3))*zCt(:,2,ik))*exp(-0.5d0*zI*dt*eps_t(3))
    zCt(:,2,ik)=zc1*zEig(:,1)+zc2*zEig(:,2)+zc3*zEig(:,3)
    
  end do

  return
end subroutine dt_evolve
