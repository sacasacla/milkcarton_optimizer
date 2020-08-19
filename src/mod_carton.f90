module mod_carton
use,intrinsic :: ISO_FORTRAN_ENV, only : REAL64
use mod_particle
use mod_particle_object
use mod_potential
implicit none
private
public :: carton
!
  real(REAL64),parameter     :: DEF_LX     =  7.00D-2   ! m
  real(REAL64),parameter     :: DEF_LY     =  7.00D-2   ! m
  real(REAL64),parameter     :: DEF_LZ     = 28.00D-2   ! m
  real(REAL64),parameter     :: DEF_DX     =  1.75D-2   ! m
!
  type,extends(particle_object) :: carton
    integer,private       :: nx = 0
    integer,private       :: ny = 0
    integer,private       :: nz = 0
    integer,private       :: n1 = 0      ! nx+ny
    real(REAL64)          :: lx = DEF_LX
    real(REAL64)          :: ly = DEF_LY
    real(REAL64)          :: lz = DEF_LZ
    real(REAL64)          :: dx = DEF_DX
    real(REAL64)          :: h0 = 0D0
    real(REAL64)          :: v0 = 1D0    ! fluid volume ( L )
    type(tension)         :: tf
    type(pressure)        :: pf
  contains
    procedure :: init       => carton_init
    procedure :: clone      => carton_clone
    procedure :: calc_force => carton_calc_force
    procedure :: volume     => carton_volume
    procedure :: surface    => carton_surface
    procedure :: dump       => carton_dump
    procedure :: load       => carton_load
    procedure :: export     => carton_export
  end type carton
!
contains
!
  subroutine carton_init(this)
  class(carton),intent(inout) :: this
  real(REAL64)                :: q(3)
  integer                     :: i,j,k
!
    this%nx   = int( this%lx / this%dx ) / 2
    this%ny   = int( this%ly / this%dx ) / 2 + 1
    this%nz   = int( this%lz / this%dx ) + 1
    this%n1   = this%nx + this%ny
    this%tf%r = this%dx
!
    call this%malloc( this%n1 * this%nz )
!
    k    = 0
    q(:) = 0D0
!
    do j=1,this%nz
!
      q(1)  = 0D0
      q(2)  = (this%ny-1) * this%dx
!
      do i=1,this%nx
        k = k + 1
        this%p(k)%q(:) = q(:)
        q(1)           = q(1) + this%dx
      enddo
!
      do i=this%nx+1,this%nx+this%ny
        k = k + 1
        this%p(k)%q(:) = q(:)
        q(2)           = q(2) - this%dx
      enddo
!
      q(3) = q(3) + this%dx
!
    enddo
!
  end subroutine carton_init
!
  pure elemental subroutine carton_clone(l,r)
  class(carton),intent(inout) :: l
  class(carton),intent(in)    :: r
    call l%clear()
    if(.not.ALLOCATED(r%p)) RETURN
    l%nx = r%nx
    l%ny = r%ny
    l%nz = r%nz
    l%n1 = r%n1
    l%lx = r%lx
    l%ly = r%ly
    l%lz = r%lz
    l%dx = r%dx
    l%h0 = r%h0
    l%v0 = r%v0
    l%tf = r%tf
    l%pf = r%pf
    call l%malloc( r%np() )
    l%p  = r%p
  end subroutine carton_clone
!
  subroutine carton_calc_force(this)
  class(carton),intent(inout) :: this
  integer                     :: i,j,k
!
    call this%reset_force()
!
    k = 0
!
    do j=1,this%nz
!
      k = k + 1
!
      call this%tf%calc( this%p(k), this%p(k+1) )
      this%p(k)%f(2:3) = this%p(k)%f(2:3) + this%tf%f(2:3) + this%tf%f(2:3)
      this%p(k+1)%f(:) = this%p(k+1)%f(:) - this%tf%f(:)
!
      do i=2,this%n1-2
        k = k + 1
        call this%tf%calc( this%p(k), this%p(k+1) )
        this%p(k)%f(:)   = this%p(k)%f(:)   + this%tf%f(:)
        this%p(k+1)%f(:) = this%p(k+1)%f(:) - this%tf%f(:)
      enddo
!
      k = k + 1
!
      call this%tf%calc( this%p(k), this%p(k+1) )
      this%p(k)%f(:)   = this%p(k)%f(:)   + this%tf%f(:)
      this%p(k+1)%f(1) = this%p(k+1)%f(1) - this%tf%f(1) - this%tf%f(1)
      this%p(k+1)%f(3) = this%p(k+1)%f(3) - this%tf%f(3) - this%tf%f(3)
!
      k = k + 1
!
    enddo
!
    k = 0
!
    do j=1,this%nz-1
      do i=1,this%n1
!
        k = k + 1
!
        call this%tf%calc( this%p(k), this%p(k+this%n1) )
!
        this%p(k)%f(:)    = this%p(k)%f(:)    + this%tf%f(:)
        this%p(k+this%n1)%f(:) = this%p(k+this%n1)%f(:) - this%tf%f(:)
!
      enddo
    enddo
!
    do i=1,this%np()
      this%p(i)%f(:) = this%p(i)%f(:) * this%tf%c
    enddo
!
    this%pf%h = this%surface()
    this%h0   = this%pf%h
    k         = 0
!
    do j=1,this%nz-1
!
      k = k + 1
!
      call this%pf%calc( this%p(k), this%p(k+1), this%p(k+this%n1+1), this%p(k+this%n1) )
!
      this%p(k)%f(2)           = this%p(k)%f(2)           + this%pf%f(2) + this%pf%f(2)
      this%p(k)%f(3)           = this%p(k)%f(3)           + this%pf%f(3) + this%pf%f(3)
      this%p(k+1)%f(:)         = this%p(k+1)%f(:)         + this%pf%f(:)
      this%p(k+this%n1)%f(2)   = this%p(k+this%n1)%f(2)   + this%pf%f(2) + this%pf%f(2)
      this%p(k+this%n1)%f(3)   = this%p(k+this%n1)%f(3)   + this%pf%f(3) + this%pf%f(3)
      this%p(k+this%n1+1)%f(:) = this%p(k+this%n1+1)%f(:) + this%pf%f(:)
!
      do i=2,this%n1-2
        k = k + 1
        call this%pf%calc( this%p(k), this%p(k+1), this%p(k+this%n1+1), this%p(k+this%n1) )
        this%p(k)%f(:)           = this%p(k)%f(:)           + this%pf%f(:)
        this%p(k+1)%f(:)         = this%p(k+1)%f(:)         + this%pf%f(:)
        this%p(k+this%n1)%f(:)   = this%p(k+this%n1)%f(:)   + this%pf%f(:)
        this%p(k+this%n1+1)%f(:) = this%p(k+this%n1+1)%f(:) + this%pf%f(:)
      enddo
!
      k = k + 1
      call this%pf%calc( this%p(k), this%p(k+1), this%p(k+this%n1+1), this%p(k+this%n1) )
!
      this%p(k)%f(:)           = this%p(k)%f(:)           + this%pf%f(:)
      this%p(k+1)%f(1)         = this%p(k+1)%f(1)         + this%pf%f(1) + this%pf%f(1)
      this%p(k+1)%f(3)         = this%p(k+1)%f(3)         + this%pf%f(3) + this%pf%f(3)
      this%p(k+this%n1)%f(:)   = this%p(k+this%n1)%f(:)   + this%pf%f(:)
      this%p(k+this%n1+1)%f(1) = this%p(k+this%n1+1)%f(1) + this%pf%f(1) + this%pf%f(1)
      this%p(k+this%n1+1)%f(3) = this%p(k+this%n1+1)%f(3) + this%pf%f(3) + this%pf%f(3)
!
      k = k + 1
!
    enddo
!
    do i=1,this%n1
      this%p(i)%f(:)           = 0D0
    enddo
!
    do i=this%n1*(this%nz-1)+1,this%n1*this%nz
      this%p(i)%f(:)           = 0D0
    enddo
!
  end subroutine carton_calc_force
!
  function carton_surface(this) result(res)
  class(carton),intent(in) :: this
  real(REAL64),parameter   :: rsix = 4d3 / 6d0 ! l
  real(REAL64)             :: res
  real(REAL64)             :: v,v0,v1
  integer                  :: i
!
    res = 0D0 ; if(this%np()<1) RETURN
!
    v0  = this%v0 / rsix
!
    do i=1,this%np()-this%n1,this%n1
!
      v  = volume_side( this%n1, this%p(i:), this%p(i+this%n1:) )
      v1 = res + v + volume_top( this%n1, this%p(i+this%n1:) )
!
      if( v1 > v0 )then
!
        v1  = v + res + volume_top( this%n1, this%p(i+this%n1:) ) - v0
        res =   - res + volume_top( this%n1, this%p(i:) )         + v0
!
        res = ( v1  * SUM( this%p(i:i+this%n1-1)%q(3) )             &
       &    +   res * SUM( this%p(i+this%n1:i+this%n1*2-1)%q(3) ) ) &
       &    / ( ( v1 + res ) * REAL( this%n1, REAL64 ) )
!
        RETURN
!
      endif
!
      res = res + v
!
    enddo
!
    res = this%lz
!
  end function carton_surface
!
  pure function carton_volume(this) result(res)
  class(carton),intent(in) :: this
  real(real64),parameter   :: rsix = 4d3 / 6d0 ! l
  real(real64)             :: res
  integer                  :: i
!
    res = 0D0
!
    do i=1,this%np()-this%n1,this%n1
      res = res + volume_side( this%n1, this%p(i:), this%p(i+this%n1:) )
    enddo
    res = res + volume_top(  this%n1, this%p(this%np()-this%n1+1:) )
    res = res * rsix
!
  end function carton_volume
!
  pure function volume_side(n1,p0,p1) result(res)
  integer,intent(in)         :: n1
  class(particle),intent(in) :: p0(n1),p1(n1)
  real(real64)               :: res,cross(3)
  integer                    :: i
!
    res = 0D0
!
    do i=1,n1-1
!
      cross(1) = p0(i+1)%q(2) * p1(i)%q(3) - p0(i+1)%q(3) * p1(i)%q(2)
      cross(2) = p0(i+1)%q(3) * p1(i)%q(1) - p0(i+1)%q(1) * p1(i)%q(3)
      cross(3) = p0(i+1)%q(1) * p1(i)%q(2) - p0(i+1)%q(2) * p1(i)%q(1)
!
      res      = res - DOT_PRODUCT( cross, p0(i)%q(:) )
      res      = res + DOT_PRODUCT( cross, p1(i+1)%q(:) )
!
    enddo
!
  end function volume_side
!
  pure function volume_top(n1,p) result(res)
  integer,intent(in)        :: n1
  type(particle),intent(in) :: p(n1)
  real(real64)              :: res
  integer                   :: i
!
    res = 0D0
!
    do i=1,n1-1
      res = res + p(i+1)%q(1) * p(i)%q(2) - p(i+1)%q(2) * p(i)%q(1)
    enddo
!
    res = res * SUM( p(:)%q(3) ) / REAL( n1, REAL64 )
!
  end function volume_top
!
  subroutine carton_export(this,dev)
  class(carton),intent(in) :: this
  integer,intent(in)       :: dev
  real(real64)             :: h
  integer                  :: i,j,k,l
!
    k = 0 ; l = 0
!
    do j=1,this%nz
!
      do i=1,this%n1
        k = k + 1 ; l = l + 1
        write( dev,'(i8,3(",",f16.9))') l,  this%p(k)%q(1),  this%p(k)%q(2),  this%p(k)%q(3)
      enddo
      do i=1,this%n1-1
        k = k - 1 ; l = l + 1
        write( dev,'(i8,3(",",f16.9))') l,  this%p(k)%q(1), -this%p(k)%q(2),  this%p(k)%q(3)
      enddo
      do i=1,this%n1-1
        k = k + 1 ; l = l + 1
        write( dev,'(i8,3(",",f16.9))') l, -this%p(k)%q(1), -this%p(k)%q(2),  this%p(k)%q(3)
      enddo
      do i=1,this%n1-2
        k = k - 1 ; l = l + 1
        write( dev,'(i8,3(",",f16.9))') l, -this%p(k)%q(1),  this%p(k)%q(2),  this%p(k)%q(3)
      enddo
      k = k + this%n1 - 2
!
    enddo
!
    k = 0
!
    do i=1,this%n1
      k = k + 1 ; l = l + 1
      write( dev,'(i8,3(",",f16.9))') l,  this%p(k)%q(1),  this%p(k)%q(2),  this%h0
    enddo
    do i=1,this%n1-1
      k = k - 1 ; l = l + 1
      write( dev,'(i8,3(",",f16.9))') l,  this%p(k)%q(1), -this%p(k)%q(2),  this%h0
    enddo
    do i=1,this%n1-1
      k = k + 1 ; l = l + 1
      write( dev,'(i8,3(",",f16.9))') l, -this%p(k)%q(1), -this%p(k)%q(2),  this%h0
    enddo
    do i=1,this%n1-2
      k = k - 1 ; l = l + 1
      write( dev,'(i8,3(",",f16.9))') l, -this%p(k)%q(1),  this%p(k)%q(2),  this%h0
    enddo
!
  end subroutine carton_export
!
  subroutine carton_dump(this,dev)
  class(carton),intent(in) :: this
  integer,intent(in)       :: dev
  integer                  :: ios
!
    write( dev, '(3I8)',    IOSTAT=ios, ERR=100 ) this%nx,this%ny,this%nz
    write( dev, '(6G20.10)', IOSTAT=ios, ERR=100 ) this%lx,this%ly,this%lz,this%dx,this%h0,this%v0
    write( dev, '(6G20.10)', IOSTAT=ios, ERR=100 ) this%tf%c,this%tf%r,this%tf%core,this%pf%h,this%pf%rho,this%pf%gra
!
    if(ALLOCATED(this%p))then
      write( dev, '(6G20.10)', IOSTAT=ios, ERR=100 ) this%p
    endif
!
100 RETURN
  end subroutine carton_dump
!
  subroutine carton_load(this,dev)
  class(carton),intent(inout) :: this
  integer,intent(in)          :: dev
  integer                     :: ios
!
    call this%clear()
!
    read( dev, *, IOSTAT=ios, ERR=100 ) this%nx,this%ny,this%nz
    read( dev, *, IOSTAT=ios, ERR=100 ) this%lx,this%ly,this%lz,this%dx,this%h0,this%v0
    read( dev, *, IOSTAT=ios, ERR=100 ) this%tf%c,this%tf%r,this%tf%core,this%pf%h,this%pf%rho,this%pf%gra
!
    this%n1   = this%nx + this%ny
!
    call this%malloc( this%n1 * this%nz )
!
    if(ALLOCATED(this%p))then
      read( dev, *, IOSTAT=ios, ERR=100 ) this%p
    endif
!
100 RETURN
  end subroutine carton_load
end module mod_carton
