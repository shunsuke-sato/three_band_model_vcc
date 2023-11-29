!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it,it2
  real(8) :: jav_t
  real(8) :: ss,ss_l

  call current(0,jav_t)
  jtz(0) = jav_t
  do it = 0,Nt

    ss_l = sum(abs(zCt(:,:,:))**2)/(2d0*dble(NKrz))
    call MPI_ALLREDUCE(ss_l,ss,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myrank==0)write(*,*)"norm",ss

    if(myrank == 0)write(*,*)'it=',it,'/',Nt
    call dt_evolve(it)
    call current(it+1,jav_t)
    jtz(it+1) = jav_t

  end do

  if(myrank == 0)then
    open(nf_current,file=cf_current)
    write(nf_current,"(A)")"# tt(au), Act(au), current(au)"
    do it2=0,Nt
      write(nf_current,"(999e26.16e3)")dt*it2,Act(it2),jtz(it2)
    end do
    close(nf_current)
  end if

end subroutine time_propagation
