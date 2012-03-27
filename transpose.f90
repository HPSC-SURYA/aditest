
program trans
	implicit none
	include "mpif.h"

	integer myid, numproc, ierr
	integer n, i, j, k, step, startindex, indexspan, bufflen
	parameter (n=8)

	! probably could use allocatable arrays to reduce memory duplication
	! oh well
	real*8 u(0:n+1,0:n+1), up(0:n+1,0:n+1)
	real*8, dimension(:), allocatable :: sendbuffer, recvbuffer
	real*8 a(n), b(n), c(n), d(n), x(n)
	real*8 dx, dt, dt2

	! MPI goodness
	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, numproc, ierr)
	! slab decomposition
	startindex = myid*n/numproc+1
	indexspan = n/numproc
	bufflen = n*n/numproc
	allocate(sendbuffer(bufflen))
	allocate(recvbuffer(bufflen))

	! dx might actually be 1/(n+1)
	dx = 1d0/n
	dt = 0.01d0

	! initialize things to 0
	u = 0.0
	u(0,:) = 1d0
	up = 0.0 ! u-prime
	up(0,:) = 1d0

	! these never change, precompute
	a = -1d0
	b = 2d0*(dx*dx/dt+1d0)
	c = -1d0

	call MPI_Barrier(MPI_COMM_WORLD, ierr)

	! start time-stepping
	do step=1,1
		! start in one direction
		do i=startindex,indexspan+startindex-1
			write(*,*),myid,i
			! load up d
			do k=1,n
				d(k) = u(i+1,k) + 2d0*(dx*dx/dt-1d0)*u(i,k) + u(i-1,k)
			enddo
			! call tridi
			call solve_tridiag(a,b,c,d,x,n)
			! load result into other matrix
			up(i,1:n) = x
		enddo

		! START TRANSPOSE
		
		! linearize data to send
	!	write(*,*), myid, 'start linearize'
	!	j = 1
	!	do i=startindex, indexspan
	!		do k=1,n
	!			sendbuffer(j+k-1) = up(i,k)
	!		enddo
	!		j = j + n
	!	enddo
	!	write(*,*), myid, 'end linearize'

		! END TRANSPOSE

		! now do other direction
		do j=startindex,indexspan+startindex-1
			do k=1,n
				d(k) = up(k,j+1) + 2d0*(dx*dx/dt-1d0)*up(k,j) + u(k,j-1)
			enddo
			call solve_tridiag(a,b,c,d,x,n)
			! back to original matrix
			u(1:n,j) = x
		enddo

		! transpose back

	enddo

	! output to screen
	if (myid .eq. 0) then
		do i=0,n+1
			write(*,'(10f5.1)'),up(i,0:n+1)
		enddo
	endif
	write(*,*), 'end'
	deallocate(sendbuffer)
	deallocate(recvbuffer)

	call MPI_Finalize(ierr)
	
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