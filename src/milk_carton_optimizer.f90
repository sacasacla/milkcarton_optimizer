program main
use,intrinsic :: ISO_FORTRAN_ENV, only : REAL64,STDOUT=>OUTPUT_UNIT
use mod_carton
implicit none
integer,parameter          :: NITER      = 1000000
integer,parameter          :: SNAP       = 10000
type(carton)               :: car1,car2
integer                    :: iter
!
  call car1%load(5)
! call car1%init()
! call car1%dump(10)
  call car2%clone(car1)
!
  do iter=0,NITER
!
    call car1%calc_force()
!
    if(MODULO(iter,SNAP)==0)then
      call car1%dump(iter/SNAP+11)
      print'(i8,G24.9,*(f16.9))',iter/SNAP,car1%cost(),car1%surface(),car1%volume()
      FLUSH(STDOUT)
    endif
!
    call linear_search( )
!
  enddo
!
  open(10,FILE='carton.csv')
  call car1%export(10)
  close(10)
!
contains
!
  subroutine linear_search()
  real(REAL64),parameter  :: thre    = 1D-4
  integer,parameter       :: maxiter = 100
  real(REAL64)            :: s(0:2),f(2)
  real(REAL64)            :: sr,fr ! reflect
  real(REAL64)            :: se,fe ! expand
  real(REAL64)            :: si,fi ! inner contract
  real(REAL64)            :: so,fo ! outer contract
  real(REAL64)            :: conv
  integer                 :: best,worst
  integer                 :: i
!
! nelder_mead parameters
!
!    alpha = 1d0
!    beta  = 2d0
!    gam   = 5d-1
!
    s(1)   = 0D0  ; f(1) = cost()
    s(2)   = 1D0  ; f(2) = cost(s(2))
    s(0)   = s(MINLOC(f,1))
!
    conv   = 0D0
!
    do i=1,maxiter
!
      best  = MINLOC(f,1) ; worst = MAXLOC(f,1)
!
! X_c = X_best in 1D simplex
! sr is refrected points : S = s(best) + alpha * (s(best) - s(worst))
!
      sr     = s(best) + s(best) - s(worst)
      fr     = cost( sr )
!
      if(f(best) <= fr .and. fr < f(worst))then
        s(worst) = sr ; f(worst) = fr
      elseif(fr < f(best))then
!! se is expanded point : S = s(best) + beta * (sr - s(best))
        se     = sr + sr - s(best)
        fe     = cost( se )
        if(fe < fr)then ; s(worst) = se ; f(worst) = fe
        else            ; s(worst) = sr ; f(worst) = fr
        endif
      else
!! se is inner or outer contraction point : S = s(best) + gam * (sr - s(best))
        si   = 0.5D0 * (s(best) + sr)
        fi   = cost( si )
!!  as S = s(best) - gam * (sr - s(best))
        so   = 1.5D0 * s(best) - 0.5D0 * sr
        fo   = cost( so )
        if(fi<fo)then ; s(worst) = si ; f(worst) = fi
        else          ; s(worst) = so ; f(worst) = fo
        endif
      endif
!
      s(0) = s(1)*s(1) + s(2)*s(2)
      conv = abs(conv - s(0))/s(0)
!
      if(conv < thre) EXIT
!
      conv = s(0)
!
    enddo
!
    call car1%move( s(best) )
!
  end subroutine linear_search
!
  function cost(s) result(res)
  real(REAL64),intent(in),optional :: s
  real(REAL64)                     :: res
!
    res = 0D0
!
    car2 = car1
    if(PRESENT(s)) call car2%move(s)
!
    call car2%calc_force()
    res = car2%cost()
!
  end function cost
!
! pure subroutine move_particle(p,s)
! class(particle),intent(inout) :: p(:,:)
! real(REAL64),intent(in)       :: s
! real(REAL64)                  :: s2
! integer                       :: i,j,n1
!
!   n1 = SIZE(p,1)
!
!   do j=2,NZ
!
!     p(1,j)%q(2)  = p(1,j)%q(2)  + s * p(1,j)%f(2)
!     p(1,j)%q(3)  = p(1,j)%q(3)  + s * p(1,j)%f(3)
!
!     do i=2,n1-1
!       p(i,j)%q(:) = p(i,j)%q(:) + s * p(i,j)%f(:)
!     enddo
!
!     p(n1,j)%q(1) = p(n1,j)%q(1) + s * p(n1,j)%f(1)
!     p(n1,j)%q(3) = p(n1,j)%q(3) + s * p(n1,j)%f(3)
!
!   enddo
!
! end subroutine move_particle
!
end program main

