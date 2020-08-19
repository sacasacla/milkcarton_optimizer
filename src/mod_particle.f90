module mod_particle
use,intrinsic :: ISO_FORTRAN_ENV, only : REAL64
implicit none
private
public :: particle
public ::  assignment (=)
!
  type particle
    real(REAL64) :: q(3)    = 0D0
    real(REAL64) :: f(3)    = 0D0
  end type particle
!
! type tension
!   real(REAL64) :: core    = 1D0
!   real(REAL64) :: r       = 0D0
! contains
!   procedure    :: force => tension_force
! end type tension
!
! type pressure
!   real(REAL64) :: h       = 0D0
!   real(REAL64) :: vert(3) = [0D0,0D0,1D0] / 4D0
! contains
!   procedure    :: force => pressure_force
! end type pressure
!
  interface assignment (=)
    module procedure particle_assign
  end interface assignment (=)
!
contains
!
  pure elemental subroutine particle_assign(l,r)
  class(particle),intent(inout) :: l
  class(particle),intent(in)    :: r
    l%q(:) = r%q(:)
    l%f(:) = r%f(:)
  end subroutine particle_assign
!
! pure elemental subroutine tension_force(this,p0,p1,f)
! class(tension),intent(in)              :: this
! class(particle),intent(inout)          :: p0,p1
! class(particle),intent(inout),optional :: f
! real(REAL64)                           :: v(3),r,r1,g
!
!   v(:) = p1%q(:) - p0%q(:)
!   g    = DOT_PRODUCT( v, v ) ; if(g<=0D0) RETURN
!   r    = SQRT( g )
!   r1   = 1D0 / r
!   g    = ( r - this%r ) * r1
!   r1   = this%r * r1
!   if(r1>1D0)then
!     r1   = r1 * r1
!     r1   = r1 * r1
!     r1   = r1 * r1
!     r1   = r1 * r1
!     r1   = 1D0 - r1
!   else
!     r1   = 0D0
!   endif
!
!   g    = g + this%core * r1
!   !g    = g + this%core * ( r1 - r1 * r1 )
!   v(:) = g * v(:)
!
!   if(PRESENT(f))then
!     f%f(:)  = v(:)
!   else
!     p0%f(:) = p0%f(:) + v(:)
!     p1%f(:) = p1%f(:) - v(:)
!   endif
!
! end subroutine tension_force
!
! pure elemental subroutine pressure_force(this,p0,p1,p2,p3,f)
! class(pressure),intent(in)             :: this
! class(particle),intent(inout)          :: p0,p1,p2,p3
! class(particle),intent(inout),optional :: f
! real(REAL64),parameter                 :: rsix = 1D0 / 6D0
! real(REAL64)                           :: a(3),b(3),n(3),g,s
!
!   a(:) = p0%q(:) + p1%q(:) + p2%q(:) + p3%q(:)
!   g    = this%h - DOT_PRODUCT( this%vert(:), a(:) )
!   if(g<0D0) g = 0D0
!
!   a(:) = p2%q(:) - p0%q(:)
!
!   b(:) = p1%q(:) - p0%q(:)
!   s    = SQRT( 1D0 - DOT_PRODUCT(a,b)**2 / ( DOT_PRODUCT(a,a) * DOT_PRODUCT(b,b) ) )
!
!   b(:) = p3%q(:) - p0%q(:)
!   s    = SQRT( 1D0 - DOT_PRODUCT(a,b)**2 / ( DOT_PRODUCT(a,a) * DOT_PRODUCT(b,b) ) ) + s
!
!   b(:) = p3%q(:) - p1%q(:)
!
!   n(1) = a(2) * b(3) - a(3) * b(2)
!   n(2) = a(3) * b(1) - a(1) * b(3)
!   n(3) = a(1) * b(2) - a(2) * b(1)
!
!   n(:) = rsix * g * s * n(:)
!
!   if(PRESENT(f))then
!     f%f(:) = - n(:)
!   else
!     p0%f(:) = p0%f(:) - n(:)
!     p1%f(:) = p1%f(:) - n(:)
!     p2%f(:) = p2%f(:) - n(:)
!     p3%f(:) = p3%f(:) - n(:)
!   endif
!
! end subroutine pressure_force
!
end module mod_particle
