
program serial
	implicit none
	include "mpif.h"

	integer n, i, j, k, step
	parameter (n=8)

	real*8 u(n+2,n+2), up(n+2,n+2)
	real*8 a(n), b(n), c(n), d(n), x(n)
	real*8 dx, dt, dt2
	real*8 t0, t1, t2

	! dx might actually be 1/(n+1)
	dx = 1d0/n
	dt = 0.01d0

	! initialize things to 0
	u = 0.0
	u(1,:) = 1d0
	up = 0.0 ! u-prime
	up(1,:) = 1d0

	! these never change, precompute
	a = -1d0
	b = 2d0*(dx*dx/dt+1d0)
	c = -1d0

	! start time-stepping
	t0 = MPI_Wtime()
	do step=1,1000
		! start in one direction
		do i=1,n
			! load up d
			do k=1,n
				d(k) = u(i+2,k+1) + 2d0*(dx*dx/dt-1d0)*u(i+1,k+1) + u(i,k+1)
			enddo
			d(1) = d(1) + u(i+1,1)
			d(n) = d(n) + u(i+1,n+2)
			! call tridi
			call solve_tridiag(a,b,c,d,x,n)
			! load result into other matrix
			up(i+1,2:n+1) = x
		enddo

		! now do other direction
		do j=1,n
			do k=1,n
				d(k) = up(k+1,j+2) + 2d0*(dx*dx/dt-1d0)*up(k+1,j+1) + up(k+1,j)
			enddo
			d(1) = d(1) + up(1,j+1)
			d(n) = d(n) + up(n+2,j+1)
			call solve_tridiag(a,b,c,d,x,n)
			! back to original matrix
			u(2:n+1,j+1) = x
		enddo
	enddo
	t1 = MPI_Wtime()

	t2 = t1-t0
	print*, t2

	! output to screen
	open(25,file='out_serial')
	write(25,*), n, n
	do i=2,n+1
		write(25,*),u(2:n+1,i)
	enddo
	close(25)

end program


! SHAMELESSLY stolen from wikipedia
subroutine solve_tridiag(a,b,c,v,x,n)
	implicit none
	  !        a - sub-diagonal (means it is the diagonal below the main diagonal)
	  !        b - the main diagonal
	  !        c - sup-diagonal (means it is the diagonal above the main diagonal)
	  !        v - right part
	  !        x - the answer
	  !        n - number of equations
	integer,intent(in) :: n
    real(8),dimension(n),intent(in) :: a,b,c,v
    real(8),dimension(n),intent(out) :: x
    real(8),dimension(n) :: bp,vp
    real(8) :: m
    integer i

! Make copies of the b and v variables so that they are unaltered by this sub
    bp(1) = b(1)
    vp(1) = v(1)

!The first pass (setting coefficients):
firstpass: do i = 2,n
	m = a(i)/bp(i-1)
	bp(i) = b(i) - m*c(i-1)
	vp(i) = v(i) - m*vp(i-1)
	end do firstpass
	
	x(n) = vp(n)/bp(n)
	!The second pass (back-substition)
backsub:do i = n-1, 1, -1
    x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
    end do backsub

end subroutine solve_tridiag
