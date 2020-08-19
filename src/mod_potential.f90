module mod_potential
use,intrinsic :: ISO_FORTRAN_ENV, only : REAL64
use mod_particle
implicit none
private
public ::  assignment (=)
public :: tension,pressure
!
  type potential
    real(REAL64) :: u     = 0D0
    real(REAL64) :: f(3)  = 0D0
  end type potential
!
  type,extends(potential) :: tension
    real(REAL64) :: c       = 1D0
    real(REAL64) :: core    = 1D0
    real(REAL64) :: r       = 0D0
  contains
    procedure    :: calc  => tension_calc
  end type tension
!
  type,extends(potential) :: pressure
    real(REAL64) :: rho     = 1D0
    real(REAL64) :: gra     = 1D0
    real(REAL64) :: h       = 0D0
    real(REAL64) :: vert(3) = [0D0,0D0,1D0]
  contains
    procedure    :: calc  => pressure_calc
  end type pressure
!
  interface assignment (=)
    module procedure tension_assign,pressure_assign
  end interface assignment (=)
!
contains
!
  pure elemental subroutine tension_assign(l,r)
  class(tension),intent(inout) :: l
  class(tension),intent(in)    :: r
    l%u    = r%u
    l%f    = r%f
    l%c    = r%c
    l%core = r%core
    l%r    = r%r
  end subroutine tension_assign
!
  pure elemental subroutine pressure_assign(l,r)
  class(pressure),intent(inout) :: l
  class(pressure),intent(in)    :: r
    l%u    = r%u
    l%f    = r%f
    l%rho  = r%rho
    l%gra  = r%gra
    l%h    = r%h
    l%vert = r%vert
  end subroutine pressure_assign
!
  pure elemental subroutine tension_calc(this,p0,p1)
  class(tension),intent(inout) :: this
  class(particle),intent(in)   :: p0,p1
  real(REAL64)                 :: r1
!
    this%f(:) = p1%q(:) - p0%q(:)
    this%u    = SQRT( DOT_PRODUCT( this%f, this%f ) )
    r1        = 1D0 / this%u
    this%f(:) = this%f(:) * r1
    this%u    = this%u - this%r
    r1        = this%r * r1
    if(r1>1D0)then
      r1   = r1 * r1
      r1   = r1 * r1
      r1   = r1 * r1
      r1   = r1 * r1
      r1   = 1D0 - r1
    else
      r1   = 0D0
    endif
!
    this%u    = this%u + this%core * r1
    this%f(:) = this%u * this%f(:)
!
  end subroutine tension_calc
!
  pure elemental subroutine pressure_calc(this,p0,p1,p2,p3)
  class(pressure),intent(inout) :: this
  class(particle),intent(in)    :: p0,p1,p2,p3
  real(REAL64),parameter        :: r24 = 1D0 / (6D0*4D0)
  real(REAL64)                  :: a(3),b(3),s
!
!   3-2
!   | |
!   0-1
!
    this%u = this%h + this%h
    this%u = this%u + this%u
!
    a(:)   = p0%q(:) + p1%q(:) + p2%q(:) + p3%q(:)
    this%u = this%u - DOT_PRODUCT( this%vert(:), a(:) )
!
    this%u = this%u * r24
!
    if(this%u<0D0)then
      this%u = 0D0 ; this%f(:) = 0D0 ; RETURN
    endif
!
    a(:) = p2%q(:) - p0%q(:)
    b(:) = p1%q(:) - p0%q(:)
    s    = SQRT( 1D0 - DOT_PRODUCT(a,b)**2 / ( DOT_PRODUCT(a,a) * DOT_PRODUCT(b,b) ) )
    b(:) = p3%q(:) - p0%q(:)
    s    = SQRT( 1D0 - DOT_PRODUCT(a,b)**2 / ( DOT_PRODUCT(a,a) * DOT_PRODUCT(b,b) ) ) + s
!
    b(:) = p1%q(:) - p3%q(:)
!
    this%f(1) = a(2) * b(3) - a(3) * b(2)
    this%f(2) = a(3) * b(1) - a(1) * b(3)
    this%f(3) = a(1) * b(2) - a(2) * b(1)
!
    this%f(:) = this%f(:) / SQRT( DOT_PRODUCT(this%f(:),this%f(:)) )
!
    this%u    = this%rho * this%gra * this%u * s         ! rho*g*h*s
    this%f(:) = this%u   * this%f(:)
!
  end subroutine pressure_calc
!
end module mod_potential
