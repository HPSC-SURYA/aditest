
! seeing how easy it is to shove all of the wakatani stuff
!  into a single subroutine

program tridi_wktn
	implicit none
	include "mpif.h"
	
	integer n, myid, numproc, ierr, stat(MPI_STATUS_SIZE)
	real*8, dimension(:), allocatable :: a, b, c, d, x
	integer indexspan, startindex, ii
	parameter (n=8)

	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, numproc, ierr)

	indexspan = n/numproc
	startindex = myid*indexspan+1

	allocate(a(indexspan))
	allocate(b(indexspan))
	allocate(c(indexspan))
	allocate(d(indexspan))
	allocate(x(indexspan))

	! initialize things
	a = 1d0
	b = -2d0
	c = 1d0
	d = -5d-1

	if (numproc .eq. 1) then
		call solve_tridiag(a,b,c,d,x,n)
	endif


	! output
	do ii=0,numproc-1
		if (myid .eq. ii) then
			write(*,*), "proc", myid
			write(*,'(8f8.1)'), x
		endif
		call MPI_Barrier(MPI_COMM_WORLD, ierr)
	enddo

	! too lazy to deallocate
	call MPI_Finalize(ierr)

end program tridi_wktn


! the magic happens here
subroutine tridi(a,b,c,d,x,n)
	implicit none
	include "mpif.h"
	integer n, i
	real(8), dimension(n) ::  a,b,c,d ! apparently needs the ::
	real(8), dimension(n) :: x
	real(8), dimension(n) :: bp, vp
	real(8) :: m
	
end subroutine tridi

! for reference
subroutine solve_tridiag(a,b,c,v,x,n)
	implicit none
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
