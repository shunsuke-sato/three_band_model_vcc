!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine input
  use global_variables
  implicit none

  if(myrank == 0)then
    read(*,*)kr_max,kz_max
    read(*,*)NKr, NKz
    read(*,*)Nt,dt
    read(*,*)envelope_1
    read(*,*)Iwcm2_1,omega_ev_1,tpulse_fs_1
    read(*,*)envelope_2
    read(*,*)Iwcm2_2,omega_ev_2,tpulse_fs_2
    read(*,*)Tdelay_fs
  end if

  call MPI_BCAST(kr_max,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(kz_max,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(NKr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(NKz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Nt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(envelope_1,20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Iwcm2_1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(omega_ev_1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tpulse_fs_1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(envelope_2,20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Iwcm2_2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(omega_ev_2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tpulse_fs_2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Tdelay_fs,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

end subroutine input
