
! seeing how easy it is to shove all of the wakatani stuff
!  into a single subroutine

program tridi_wktn
	implicit none
	include "mpif.h"
	
	integer n, myid, numproc, ierr, stat(MPI_STATUS_SIZE)
	real*8, dimension(:), allocatable :: a, b, c, d, x, bp, acorr, acorr2
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
	allocate(acorr2(indexspan))

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
			acorr(ii) = (-a(startindex+ii-1)/bp(startindex+ii-2)) * acorr(ii-1)
		enddo
	endif
	acorr2 = 0d0
	if (myid .ne. numproc-1) then
		acorr2(indexspan) = -c(startindex+indexspan-1)/bp(startindex+indexspan-1)
		do ii=indexspan-1,1,-1
			acorr2(ii) = (-c(startindex+ii-1)/bp(startindex+ii-1)) * acorr2(ii+1)
		enddo
	endif

	! for testing on 1 processor
	if (numproc .eq. 1) then
		call solve_tridiag(a,b,c,d,x,n)
	else
		call tridi(a,b,c,d,x,n,myid,numproc,bp,acorr,acorr2,indexspan,startindex)
	endif

	! output
	do ii=0,numproc-1
		if (myid .eq. ii) then
			write(*,*), "proc", myid
			write(*,'(8f8.2)'), x
		endif
		call MPI_Barrier(MPI_COMM_WORLD, ierr)
	enddo

	! too lazy to deallocate
	call MPI_Finalize(ierr)

end program tridi_wktn


! the magic happens here
subroutine tridi(a,b,c,d,x,n,myid,numproc,bp,acorr,acorr2,indexspan,startindex)
	implicit none
	include "mpif.h"
	integer n, myid, numproc, ii, ierr, req
	integer stat(MPI_STATUS_SIZE)
	integer indexspan, startindex
	real(8), dimension(n) :: a,b,c,d,bp ! apparently needs the ::
	real(8), dimension(indexspan) :: vp, x, acorr, acorr2
	real(8) :: m, corr, sendcorr

! Start First Pass
	! Guess Phase
	vp(1) = d(startindex)
	do ii=2,indexspan
		m = a(startindex+ii-1)/bp(startindex+ii-2)
		vp(ii) = d(startindex+ii-1) - m*vp(ii-1) ! fuck this line of code.
	enddo
	! Propagation Phase
	if (myid .ne. 0) then
        call MPI_Recv(corr, 1, MPI_REAL8, myid-1, 0, MPI_COMM_WORLD, stat, ierr)
		vp(indexspan) = vp(indexspan) + acorr(indexspan)*corr
	endif
	sendcorr = vp(indexspan) ! could have sworn this was different. whatevs.
	if (myid .ne. numproc-1) then
	    call MPI_Isend(sendcorr, 1, MPI_REAL8, myid+1, 0, MPI_COMM_WORLD, req, ierr)
	endif

	! Correction Phase
	do ii=1,indexspan-1
		vp(ii) = vp(ii) + acorr(ii)*corr
	enddo
	if (myid .ne. numproc-1) then
		call MPI_Wait(req, stat, ierr)
	endif

! Start Second Phase
	! Guess Phase
	x(indexspan) = vp(indexspan)/bp(indexspan+startindex-1)
	do ii=indexspan-1,1,-1
		x(ii) = vp(ii)/bp(startindex+ii-1) - (c(startindex+ii-1)/bp(startindex+ii-1))*x(ii+1)
	enddo
	! Propagation Phase
	if (myid .ne. numproc-1) then
        call MPI_Recv(corr, 1, MPI_REAL8, myid+1, 0, MPI_COMM_WORLD, stat, ierr)
		x(1) = x(1) + acorr2(1)*corr
	endif
	sendcorr = x(1) ! could have sworn this was different. whatevs.
	if (myid .ne. 0) then
	    call MPI_Isend(sendcorr, 1, MPI_REAL8, myid-1, 0, MPI_COMM_WORLD, req, ierr)
	endif

	! Correction Phase
	do ii=2,indexspan
		x(ii) = x(ii) + acorr2(ii)*corr
	enddo
	if (myid .ne. 0) then
		call MPI_Wait(req, stat, ierr)
	endif

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
