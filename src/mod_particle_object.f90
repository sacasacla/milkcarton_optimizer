module mod_particle_object
use,intrinsic :: ISO_FORTRAN_ENV, only : REAL64
use mod_particle
implicit none
private
public :: particle_object
public ::  assignment (=)
!
  type particle_object
    type(particle),allocatable :: p(:)
  contains
    procedure :: malloc      => particle_object_malloc
    procedure :: merge_force => particle_object_merge_force
    procedure :: reset_force => particle_object_reset_force
    procedure :: cost        => particle_object_cost
    procedure :: move        => particle_object_move
    procedure :: np          => particle_object_np
    procedure :: clear       => particle_object_clear
    final     :: particle_object_destroy
  end type particle_object
!
  interface assignment (=)
    module procedure particle_object_assign
  end interface assignment (=)
!
contains
!
  pure elemental subroutine particle_object_malloc(this,n)
  class(particle_object),intent(inout) :: this
  integer,intent(in)                   :: n
    call this%clear() ; ALLOCATE( this%p(n) )
  end subroutine particle_object_malloc
!
  pure elemental subroutine particle_object_assign(l,r)
  class(particle_object),intent(inout) :: l
  class(particle_object),intent(in)    :: r
    if(ALLOCATED(r%p)) l%p = r%p
  end subroutine particle_object_assign
!
  pure elemental subroutine particle_object_merge_force(l,r)
  class(particle_object),intent(inout) :: l
  class(particle_object),intent(in)    :: r
  integer                              :: i
    do i=1,l%np()
      l%p(i)%f(:) = l%p(i)%f(:) + r%p(i)%f(:)
    enddo
  end subroutine particle_object_merge_force
!
  pure elemental subroutine particle_object_reset_force(this)
  class(particle_object),intent(inout) :: this
  integer                              :: i
    do i=1,this%np()
      this%p(i)%f(:) = 0D0
    enddo
  end subroutine particle_object_reset_force
!
  pure elemental function particle_object_np(this) result(res)
  class(particle_object),intent(in) :: this
  integer                           :: res
    if(ALLOCATED(this%p))then
      res = SIZE(this%p)
    else
      res = 0
    endif
  end function particle_object_np
!
  function particle_object_cost(this) result(res)
  class(particle_object),intent(inout) :: this
  real(REAL64)                         :: res
  integer                              :: i
!
    res = 0D0
!
    if(this%np()<1) RETURN
!
    do i=1,this%np()
      res = res + DOT_PRODUCT( this%p(i)%f, this%p(i)%f )
    enddo
!
    res = res / REAL(this%np(),REAL64)
!
    if(res/=res) res = HUGE(0D0)
!
  end function particle_object_cost
!
  pure subroutine particle_object_move(this,s)
  class(particle_object),intent(inout) :: this
  real(REAL64),intent(in)              :: s
  integer                              :: i
!
    do i=1,this%np()
      this%p(i)%q(:) = this%p(i)%q(:) + s * this%p(i)%f(:)
    enddo
!
  end subroutine particle_object_move
!
  pure elemental subroutine particle_object_clear(this)
  class(particle_object),intent(inout) :: this
    if(ALLOCATED(this%p)) DEALLOCATE(this%p)
  end subroutine particle_object_clear
!
  pure elemental subroutine particle_object_destroy(this)
  type(particle_object),intent(inout) :: this
    call this%clear()
  end subroutine particle_object_destroy
!
end module mod_particle_object
