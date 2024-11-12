! GNU compiler: gfortran-9.2 -O3 -fdefault-real-16 -o cheb_cfs cheb_cfs.f90
! Run: ./cheb_cfs

!===============================================================================
! This program computes Chebyshev coefficients for functions q_n(x), where
! n = 4,5,6,7,...,10000 and x from the interval [n-1,n+1].
! The coefficients are stored in the file cfs.dat
! There are m rows for each n; the first column contains n and the second column
! contains the corresponding Chebyshev coefficients.
!===============================================================================


program Chebyshev

	implicit none

	integer :: m,mm,i,j,n,p,k
	real (kind=16) :: gq,gn
	real (kind=16), dimension (:), allocatable :: a1,d1,c1,a2,d2,c2
	real (kind=16), dimension (:,:), allocatable :: cosi
	real (kind=16), parameter :: PI=3.1415926535897932384626433832795029_16
	external gq,gn

	real (kind=4) :: ts,tf

	open(unit=10,file='cfs.dat',status='replace',action='write',&
	     form='formatted',position='asis') ! Open file for data record.

	call cpu_time(ts) ! Measure computation time.

! Enter the maximal Euclidean space dimension n <= 10000.
  do
		print '(/,a/)', 'Enter even integer from 4 to 10000'
		read *, n
		if (n >= 4) exit
	end do

! Define number of Chebyshev coefficients:
	m = 54 ! m = 54 gives minimal Chebyshev truncation error for n = 4 and 5,
	! which is of order |a_m|, which is just above machine accuracy ~ 1.9e-34

! NOTE: We shall use computed Chebyshev coefficients in Python program with
! double precision computation. Thus, we'll save only the first 24 coefficients
! and restrict the number of significant digits corresponding to np.float128()

mm = 24

! Compute Chebyshev expansion matrix:
	allocate(cosi(1:m,1:m))

	do i=1,m
		do j=1,m
			cosi(i,j) = cos((PI*(i-1))*((j-0.5_16)/m))
		end do
	end do

	! Initialise arrays:
	allocate(a1(1:m))
	allocate(d1(1:m))
	allocate(c1(1:m))
	allocate(a2(1:m))
	allocate(d2(1:m))
	allocate(c2(1:m))

  ! Compute Chebyshev coefficients a(i) for n = 4:
	p = 4
	call chebgq(cosi,d1,m,p,gq)
	call cheba(d1,a1,m,p)

	do i=1,mm
		write (10,'(i5,f29.20)') p, a1(i)
	end do

	! Compute Chebyshev coefficients a(i) for n = 5:
	call chebgq(cosi,d2,m,p+1,gq)
	call cheba(d2,a2,m,p+1)
	do i=1,mm
		write (10,'(i5,f29.20)') p+1, a2(i)
	end do

	! Compute Chebyshev coefficients a(i) for n = p+2,p+3,...:

	do k=p+2,n,2
		call chebgn(cosi,c1,m,k,gn) ! compute c(i)
		call dca(d1,c1,a1,m) ! compute d(i)
		call cheba(d1,a1,m,k) ! compute a(i)
		! Record Chebyshev coefficients a(i) for n = k:
		do i=1,mm
			write (10,'(i5,f29.20)') k, a1(i)
		end do
		if (k < n) then
			call chebgn(cosi,c2,m,k+1,gn) ! compute c(i)
			call dca(d2,c2,a2,m) ! compute d(i)
			call cheba(d2,a2,m,k+1) ! compute a(i)
			! Record Chebyshev coefficients a(i) for n = k:
			do i=1,mm
				write (10,'(i5,f29.20)') k+1, a2(i)
			end do
		end if
	end do

	call cpu_time(tf)
	write(*,'(/,a,f10.7,1x,a,/)') 'Chebyshev coefficients computation time:', &
	tf-ts, 'seconds.'

	close (unit=10,status='keep')

end program Chebyshev


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
