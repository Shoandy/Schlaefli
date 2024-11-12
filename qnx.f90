! GNU compiler: gfortran-9.2 -O3 -fdefault-real-16 -o qnx qnx.f90
! Run: ./qnx

!===============================================================================
! This program computes function q_n(x) for chosen n and x from [n-1,n+1]
! and saves the values x and q_n(x) in the file qnx.dat
!===============================================================================

program schlaefli

	implicit none

	integer :: m,i,j,n,p,k
	real (kind=16) :: y,gq,gn,x,chebev
	real (kind=16), dimension (:), allocatable :: c,a,d
	real (kind=16), dimension (:,:), allocatable :: cosi
	real (kind=16), parameter :: PI=3.1415926535897932384626433832795029_16
	external gq,gn,chebev

	real (kind=4) :: ts,tf

	call cpu_time(ts) ! Measure computation time.

	open(unit=10,file='qnx.dat',status='replace',&
	action='write',form='formatted', position='asis') ! Record x and q_n(x)

! Dimension of Euclidean space n: n = 4,5,...,10000
  do
		print '(/,a/)', 'Enter integer n: 4,5...,10000'
		read *, n
		if (n >= 4) exit
	end do

	write(10,'(i5)') n

! Number of Chebyshev coefficients:
	m = 54 ! N = 54 gives minimal Chebyshev truncation error for n = 4 and 5,
	! which is of order |a_N| which is just above the machine accuracy ~ 1.9e-34

! Compute Chebyshev approximation matrix:
	allocate(cosi(1:m,1:m))

	do i=1,m
		do j=1,m
			cosi(i,j) = cos((PI*(i-1))*((j-0.5_16)/m))
		end do
	end do

	if (mod(n,2) == 0) then
		p = 4
	else
		p = 5
	end if

	! Compute Chebyshev coefficients a(i) of Q_(p-2)(y):
	allocate(a(1:m))
	allocate(d(1:m))

  call chebgq(cosi,d,m,p,gq)
	call cheba(d,a,m,p)

	! If n>p, compute Chebyshev coefficients a(i) of Q_n(y):
  if (n>p) then
		allocate(c(1:m))
		do k=p,n-2,2
			call chebgn(cosi,c,m,k+2,gn) ! compute c(i) of G_n(y)
			call dca(d,c,a,m) ! compute d(i) of GQ_n(y)
			call cheba(d,a,m,k+2) ! compute a(i) of Q_n(y)
		end do
	end if

! Compute and save x and q_n(x) repeatedly, for any x in [n-1,n+1]:
  do i=0,20000
		x = n - 1.0_16 + i*0.0001_16
		y = x - n
		write(10,'(f11.3,f39.34)') x, chebev(a,m,y)
	end do

	call cpu_time(tf)
	write(*,'(/,a,f10.7,1x,a,/)') 'Computation time:', tf-ts, 'seconds.'

	close (unit=10,status='keep')

end program schlaefli


! Subroutines:

! Function GQ to be Chebyshev approximated:
function gq(y,p)

	integer :: p
	real (kind=16) :: y,gq
	real (kind=16), parameter :: PI=3.1415926535897932384626433832795029_16

  if (p == 4) then
		gq = 18.0_16*acos(1.0_16/(y+2.0_16))/(y+4.0_16) &
		   / sqrt((y+4.0_16)**2-1.0_16)/sqrt(y+1.0_16)
	else
		gq = 96.0_16*sqrt(5.0_16)*(acos(1.0_16/(y+3.0_16))-PI/3.0_16) &
		   / (y+1.0_16)/(y+5.0_16)/sqrt((y+5.0_16)**2-1.0_16)
	end if

	return
end function gq

! Subroutine Chebyshev fit: chebft
! Given a function func = f(y) on the interval [-1,1],
! and a maximum degree m, this routine computes the m coefficients c_{k}
! such that $f(y) \approx -c_{1}/2 + \sum_{k=1}^{m}c_{k}T_{k-1}(y)$.
! This routine is to be used with moderately large m (e.g., 30 or 50).
! Parameters: maximum expected value of m and $\pi$.

subroutine chebgq(cosi,c,m,p,func)

	integer :: m,p
	integer, parameter :: MAX=70
	real (kind=16) :: c(m),func
	real (kind=16) :: cosi(m,m)
	external func
	integer :: j,k
	real (kind=16) :: fac,y,f(MAX),sum

! We evaluate the function at m points,
! as required by Chebyshev approximation.
	do k=1,m
		f(k) = func(cosi(2,k),p)
	end do

	fac = 2.0_16/m
	do j=1,m
		sum = 0.0_16
		do k=1,m
			sum = sum+f(k)*cosi(j,k)
		end do
		c(j) = fac*sum
	end do
return
end subroutine chebgq

! Subroutine evaluatinging Chebyshev coefficients a(i) of function f(y),
! which solves the differential equation
! 2*(1+y)*Q_n'+(n-1)*Q_n = GQ_n(y), -y<=x<=y,
! from the Chebyshev coefficients d(i) of the function GQ_n(y).

subroutine cheba(d,a,m,n)

	integer :: m,n
	real (kind=16) :: d(m),a(m)
	integer :: i

  a(m) = d(m)/(2*m+n-3)
	a(m-1) = (d(m-1)-4*(m-1)*a(m))/(2*m+n-5)

	do i=m-2,1,-1
		a(i) = (d(i)-d(i+2)+(n-2*i-3)*a(i+2)-4*i*a(i+1))/(2*i+n-3)
	end do

	return
end subroutine cheba

! Function G_n(y) to be Chebyshev approximated:
function gn(n,y)

	integer :: n
	real (kind=16) :: y,gn

	gn = (n-1.0_16)**2*sqrt(n*(n-2.0_16))/(y+n)/sqrt((y+n)**2-1.0_16)

	return
end function gn

! Compute m Chebyshev coefficients c(i) of G_n(y):
subroutine chebgn(cosi,c,m,n,func)

	integer :: n,m
	integer, parameter :: MAX=70
	real (kind=16) :: c(m),func
	real (kind=16) :: cosi(m,m)
	external func
	integer :: j,k
	real (kind=16) :: fac,y,f(MAX),sum

! We evaluate the function at m points,
! as required by Chebyshev approximation.
	do k=1,m
		y = cosi(2,k)
		f(k) = func(n,y)
	end do

	fac = 2.0_16/m
	do j=1,m
		sum = 0.0_16
		do k=1,m
			sum = sum+f(k)*cosi(j,k)
		end do
		c(j) = fac*sum
	end do

return
end subroutine chebgn

! Compute Chebyshev coefficients d(i) of GQ_n(y) from c(i) and a(i):

subroutine dca(d,c,a,m)

	integer :: m
	real (kind=16) :: a(m),c(m),d(m)
	integer i,k

	d(1) = c(1)*a(1)/2
	d(2) = 0.0_16

	do i=2,m
		d(1) = d(1)+c(i)*a(i)
		d(2) = d(2)+(c(i)*a(i-1)+c(i-1)*a(i))/2
	end do

	do k=3,m
		d(k) = 0.0_16
		do i=2,k-1
			d(k) = d(k)+c(i)*a(k+1-i)/2
		end do
		do i=k,m
			d(k) = d(k)+(c(i)*a(i-k+1)+c(i-k+1)*a(i))/2
		end do
	end do

	return
end subroutine dca

! Chebyshev evaluation: chebev
! All arguments are input. c(1:m) is an array of Chebyshev coefficients,
! the first m elements of c output from chebft or chint.
! The Chebyshev polynomial $-c_{1}/2 + \sum_{k=1}^{m}c_{k}T_{k-1}(y)$,
! is evaluated at y and the result is returned as the function value.

function chebev(c,m,y)
	integer :: m
	real (kind=16) :: chebev,y,c(m)
	integer :: j
	real (kind=16) :: d,dd,sv,y2

	d = 0.0_16
	dd = 0.0_16
	y2 = 2.0_16*y
	do j=m,2,-1 ! Clenshaw's recurrence.
		sv = d
		d = y2*d-dd+c(j)
		dd = sv
	end do

	chebev = y*d-dd+0.5_16*c(1)

return
end function chebev
