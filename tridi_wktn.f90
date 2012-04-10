
! seeing how easy it is to shove all of the wakatani stuff
!  into a single subroutine

program tridi_wktn
	implicit none
	include "mpif.h"
	
	integer n, myid, numproc, ierr, stat(MPI_STATUS_SIZE)
	real*8, dimension(:), allocatable :: a, b, c, d, x, bp, acorr
	integer indexspan, startindex, ii
	parameter (n=8)

	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, numproc, ierr)

	indexspan = n/numproc
	startindex = myid*indexspan+1

	allocate(a(n))
	allocate(b(n))
	allocate(c(n))
	allocate(d(n))
	allocate(x(indexspan))
	allocate(bp(n))
	allocate(acorr(indexspan))

	! initialize things
	a = 1d0
	b = -2d0
	c = 1d0
	d = -5d-1

	! we want to precompute bp
	bp(1) = b(1)
	do ii=2,n
		bp(ii) = b(ii) - (a(ii)/bp(ii-1))*c(ii-1)
	enddo

	! we also want to pre-compute the correction coefs for the entire vector
	acorr = 0d0
	if (myid .ne. 0) then
		acorr(1) = -a(startindex)/bp(startindex-1)
		do ii=2,indexspan
			acorr(ii) = acorr(ii-1) * (-a(startindex+ii-1)/bp(startindex+ii-1))
		enddo
	endif

	! for testing on 1 processor
!	if (numproc .eq. 1) then
!		call solve_tridiag(a,b,c,d,x,n)
!	endif

	call tridi_condition(a,b,c,d,n,bp)


	! output
	do ii=0,numproc-1
		if (myid .eq. ii) then
			write(*,*), "proc", myid
			write(*,'(4f8.1)'), x
		endif
		call MPI_Barrier(MPI_COMM_WORLD, ierr)
	enddo

	! too lazy to deallocate
	call MPI_Finalize(ierr)

end program tridi_wktn


subroutine tridi_condition(a,b,c,d,n,bp)
	implicit none
	integer n, ii
	real(8), dimension(n) :: a,b,c,d,bp

	bp(1) = b(1)
	do ii=2,n
		bp(ii) = b(ii) - (a(ii)/bp(ii-1))*c(ii-1)
	enddo

end subroutine tridi_condition

! the magic happens here
subroutine tridi(a,b,c,d,x,n,myid,numproc)
	implicit none
	include "mpif.h"
	integer n, myid, numproc, ii, ierr, stat(MPI_STATUS_SIZE), req
	real(8), dimension(n) ::  a,b,c,d ! apparently needs the ::
	real(8), dimension(n) :: x, acorr
	real(8), dimension(n) :: bp, vp
	real(8) :: m, corr, sendcorr

	! do first non-wakatani pass
	bp(1) = b(1)
	do ii=2,n
		m = a(ii)/bp(ii-1)
        bp(ii) = b(ii) - m*c(ii-1)
	enddo

! Start First Pass
	! Guess Phase
	if (myid .eq. 0) then
		vp(1) = d(1)
	else
		vp(1) = 0d0
	endif
	acorr(1) = 1d0
	do ii=2,n
		m = a(ii)/bp(ii-1)
		vp(ii) = d(ii) - m*vp(ii-1)
		acorr(ii) = acorr(ii-1) * (-m) ! prepare correction factor
	enddo

	! Propagation Phase
	if (myid .ne. 0) then
        call MPI_Recv(corr, 1, MPI_REAL8, myid-1, 0, MPI_COMM_WORLD, stat, ierr)
		vp(n) = vp(n) + acorr(n)*corr
	endif
	if (myid .ne. numproc-1) then
		sendcorr = vp(n)
	    call MPI_Isend(sendcorr, 1, MPI_REAL8, myid+1, 0, MPI_COMM_WORLD, req, ierr)
	endif

	! Correction Phase
!	acorr = 1d0
!	do ii=1,indexspan-1
!		acorr = acorr * a(ii)
!		w(ii) = w(ii) + acorr*corr
!	enddo
!	call MPI_Wait(req, stat, ierr)



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
