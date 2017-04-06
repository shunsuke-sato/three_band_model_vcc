!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
! Reference Joachim Kopp, Int. J. Mod. Phys. C 19 (2008) 523-548
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine diag3x3(zA,zB,lambda)
  implicit none
  real(8),parameter :: pi=4d0*atan(1d0)
  complex(8), intent(in) :: zA(3,3)
  complex(8), intent(out) :: zB(3,3)
  real(8), intent(out) :: lambda(3)

  real(8) :: c0,c1,c2
  real(8) :: p,q,t,ss,phi,xx(3),xx_t(3)
  complex(8) :: zC(3,3),zx,zy
  integer :: ivec,imin

  c2 = -real(zA(1,1)+zA(2,2)+zA(3,3))
  c1 = real( zA(1,1)*zA(2,2) + zA(2,2)*zA(3,3) + zA(3,3)*zA(1,1) &
       -abs(zA(1,2))**2 -abs(zA(1,3))**2 -abs(zA(2,3))**2 )
  c0 = real( zA(1,1)*abs(zA(2,3))**2 + zA(2,2)*abs(zA(1,3))**2 &
       + zA(3,3)*abs(zA(1,2))**2 - zA(1,1)*zA(2,2)*zA(3,3) &
       - 2d0*conjg(zA(1,3))*zA(1,2)*zA(2,3) )

  p = c2**2 - 3d0*c1
  q = -27d0*c0/2d0 -c2**3 + 9d0*c2*c1/2d0
  ss = c1**2*(p-c1)/4d0+c0*(q+27d0*c0/4d0)
  ss = sqrt(27d0*ss)/q
  phi = atan(ss)/3d0
  if(phi < 0)phi = phi + pi

  xx(1) = 2d0*cos(phi)
  xx(2) = 2d0*cos(phi+2d0*pi/3d0)
  xx(3) = 2d0*cos(phi-2d0*pi/3d0)

  lambda(:) = (sqrt(p)*xx(:)-c2)/3d0

  do ivec = 1,3

    zC = zA
    zC(1,1) = zC(1,1) -lambda(ivec)
    zC(2,2) = zC(2,2) -lambda(ivec)
    zC(3,3) = zC(3,3) -lambda(ivec)

    xx(1) = abs(zC(1,1))
    xx(2) = abs(zC(2,2))
    xx(3) = abs(zC(3,3))
    if(xx(1)<xx(2))then
      if(xx(1)<xx(3))then
        imin = 1
      else
        imin = 3
      end if
    else
      if(xx(2)<xx(3))then
        imin=2
      else
        imin=3
      end if
    end if

    select case(imin)
    case(1)
! First vector
      zx = (conjg(zC(1,3))*zC(2,3) - zC(3,3)*conjg(zC(1,2))) &
        /(zC(2,2)*zC(3,3)-abs(zC(2,3))**2)
      zy = -(conjg(zC(1,3)) + conjg(zC(2,3))*zx)/zC(3,3)
  
      ss = 1d0 + abs(zx)**2 + abs(zy)**2; ss = 1d0/sqrt(ss)
      zB(1,ivec) = ss; zB(2,ivec) = zx*ss; zB(3,ivec) = zy*ss
    case(2)
! Second vector
      zx = (conjg(zC(1,3))*zC(1,2) - zC(1,1)*conjg(zC(2,3))) &
        /(zC(1,1)*zC(3,3)-abs(zC(1,3))**2)
      zy = -(zC(1,2) + zC(1,3)*zx)/zC(1,1)
      ss = 1d0 + abs(zx)**2 + abs(zy)**2; ss = 1d0/sqrt(ss)
      zB(1,ivec) = zy*ss; zB(2,ivec) = ss; zB(3,ivec) = zx*ss
    case(3)
! Third vector
      zx = (zC(1,2)*zC(2,3) - zC(2,2)*zC(1,3)) &
        /(zC(1,1)*zC(2,2)-abs(zC(1,2))**2)
      zy = -(zC(2,3) + conjg(zC(1,2))*zx)/zC(2,2)
      ss = 1d0 + abs(zx)**2 + abs(zy)**2; ss = 1d0/sqrt(ss)
      zB(1,ivec) = zx*ss; zB(2,ivec) = zy*ss; zB(3,ivec) = ss
    case default
      stop "Failed in diag3x3"
    end select
  end do
end subroutine diag3x3
