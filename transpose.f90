
program trans
	implicit none
	include "mpif.h"

	integer myid, numproc, ierr, stat(MPI_STATUS_SIZE), dt_slab, req
	integer n, i, j, k, step, startindex, indexspan, stage, slabsize, maxstep
	integer id1, id2
!	parameter (n=8)

	! probably could use allocatable arrays to reduce memory duplication
	! oh well
	real*8, dimension(:,:), allocatable :: slab1, slab2, slab3, slab4
	real*8, dimension(:), allocatable :: a,b,c,d,x
!	real*8 a(n), b(n), c(n), d(n), x(n)
	real*8 dx, dt, dt2

	! timing stuff
	real*8 t0, t1, t2, avgtime
	character *100 buffer

	call getarg(1,buffer)
	read(buffer,*) n
	call getarg(2,buffer)
	read(buffer,*) maxstep

	! MPI goodness
	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, numproc, ierr)
	! slab decomposition
	startindex = myid*n/numproc+2
	indexspan = n/numproc
	slabsize = (n+2)*(indexspan+2)
	! both slabs are columnar to speed things up? sure.
	allocate(slab1(n+2,indexspan+2))
	allocate(slab2(n+2,indexspan+2))
	allocate(slab3(n+2,indexspan+2))
	allocate(slab4(n+2,indexspan+2))
	allocate(a(n))
	allocate(b(n))
	allocate(c(n))
	allocate(d(n))
	allocate(x(n))

	call MPI_TYPE_CONTIGUOUS(slabsize, MPI_DOUBLE_PRECISION, dt_slab, ierr)
	call MPI_TYPE_COMMIT(dt_slab, ierr)

	! dx might actually be 1/(n+1)
	dx = 1d0/n
	dt = 0.01d0

	! initialize things to 0
	slab1 = 0.0
	slab2 = 0.0
	slab3 = 0.0
	slab4 = 0.0
	call boundary(myid, numproc, slab1, n, indexspan)
	call boundary(myid, numproc, slab2, n, indexspan)
	call boundary2(myid, numproc, slab3, n, indexspan)
	call boundary2(myid, numproc, slab4, n, indexspan)

	! these never change, precompute
	a = -1d0
	b = 2d0*(dx*dx/dt+1d0)
	c = -1d0

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t0 = MPI_WTIME()
	! start time-stepping
	do step=1,maxstep
		! start in one direction
		do i=2,indexspan+1
			!i2 = i+startindex-1
			! load up d
			do k=2,n+1
				d(k-1) = slab1(k,i+1) + 2d0*(dx*dx/dt-1d0)*slab1(k,i) + slab1(k,i-1)
			enddo
			! deal with boundary conditions
			d(1) = d(1) + slab1(1,i)
			d(n) = d(n) + slab1(n+2,i)
			! call tridi
			call solve_tridiag(a,b,c,d,x,n)
			! load result into other matrix
			slab2(2:n+1,i) = x
		enddo
		! START TRANSPOSE
		call stransposeMPI(myid, numproc, n, indexspan, slab2, slab3)
		! END TRANSPOSE
		! now do other direction
		do j=2,indexspan+1
			do k=2,n+1
				d(k-1) = slab3(k,j+1) + 2d0*(dx*dx/dt-1d0)*slab3(k,j) + slab3(k,j-1)
			enddo
			d(1) = d(1) + slab3(1,j)
			d(n) = d(n) + slab3(n+2,j)
			call solve_tridiag(a,b,c,d,x,n)
			! back to original matrix
			slab4(2:n+1,j) = x
		enddo

		! transpose back
		call stransposeMPI(myid, numproc, n, indexspan, slab4, slab1)
	enddo
	t1 = MPI_WTIME();

	t2 = t1-t0

	call MPI_ALLREDUCE(t2, avgtime, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
	if (myid .eq. 0) then
		write(*,*), avgtime/numproc
	endif

	! output to file for legit comparing
	if (myid .eq. 0) then
		open(25,file='out_trans')
		write(25,*), n, n
		do i=2,indexspan+1
			write(25,*),slab1(2:n+1,i)
		enddo
		do i=1,numproc-1
			call MPI_RECV(slab1, 1, dt_slab, i, 5, MPI_COMM_WORLD, stat, ierr)
			do j=2,indexspan+1
				write(25,*),slab1(2:n+1,j)
			enddo
		enddo
		close(25)
	else
		call MPI_SEND(slab1, 1, dt_slab, 0, 5, MPI_COMM_WORLD, ierr)
	endif

	deallocate(slab1)
	deallocate(slab2)
	deallocate(slab3)
	deallocate(slab4)
	call MPI_TYPE_FREE(dt_slab, ierr)
	call MPI_Finalize(ierr)
end program



! simple (crappy) transpose
subroutine stransposeMPI(myid, numproc, slablength, blocklength, sendslab, recvslab)
	include "mpif.h"
	integer myid, numproc, slablength, blocklength, offset
	real*8 sendslab(slablength+2, blocklength+2), recvslab(slablength+2, blocklength+2)
	real*8 block(blocklength+2, blocklength)
	integer ii, ij, ii2, ik, ierror, stat(MPI_STATUS_SIZE)
	integer blocksize

	blocksize = blocklength*(blocklength+2)

	do ii2=0,numproc-1
		ii = ii2 ! something terrible is happening
		if (ii .eq. myid) then
			! send everything
			do ij=0,numproc-1
				if (ij .ne. myid) then
	!				write(*,*), "proc ", myid, "sending to ", ij
					offset = ij*blocklength
					block = sendslab(1+offset:2+blocklength+offset, 2:1+blocklength)
					call MPI_SEND(block, blocksize, MPI_REAL8, ij, 0, MPI_COMM_WORLD, ierror)
	!				write(*,*), "proc ", myid, "sent to ", ij
				else
					offset = ij*blocklength
					block = sendslab(1+offset:2+blocklength+offset, 2:1+blocklength)
					recvslab(2+offset:1+blocklength+offset, 1:2+blocklength) = transpose(block)
				endif
			enddo
		else
			! recv from ii
	!		write(*,*), "proc ", myid, "post recv from ", ii
			call MPI_RECV(block, blocksize, MPI_REAL8, ii, 0, MPI_COMM_WORLD, stat, ierror)
			offset = ii*blocklength
			recvslab(2+offset:1+blocklength+offset, 1:2+blocklength) = transpose(block)
	!		write(*,*), "proc ", myid, "recv from ", ii
		endif
		call MPI_BARRIER(MPI_COMM_WORLD, ierror)
	enddo
end subroutine stransposeMPI






! this is the advanced (faster) version that doesnt work at all
subroutine transposeMPI(myid, numproc)
	include "mpif.h" ! crashes if i dont have this again. why?
	! params
	integer myid, numproc
	! local variables
	integer ii, tradewith, ierror, stat
	integer sendbuf, recvbuf

	ii = 1
	i2 = 2**(ii-1)
	sendbuf = 10
	do while (ii .le. numproc/2)
		if (ii .eq. numproc/2) then
			! final stage, no need for 2 sub-stages
		else
			! non-final stage, needs 2 sub-stages
			if (iand(myid,i2) .eq. 0) then
				tradewith = mod(myid+ii, numproc)
				write(*,*), "proc ", myid, " trade up with ", tradewith, " at stage ", i2
				call MPI_SEND(sendbuf, 1, MPI_INT, tradewith, 0, MPI_COMM_WORLD, ierror)
				call MPI_RECV(recvbuf, 1, MPI_INT, tradewith, 0, MPI_COMM_WORLD, stat, ierror)
			else
				tradewith = mod(myid-ii+numproc, numproc)
				write(*,*), "proc ", myid, " trade down with ", tradewith, " at stage ", i2
				call MPI_RECV(recvbuf, 1, MPI_INT, tradewith, 0, MPI_COMM_WORLD, stat, ierror)
				call MPI_SEND(sendbuf, 1, MPI_INT, tradewith, 0, MPI_COMM_WORLD, ierror)
			endif
		endif
		write(*,*), "proc ", myid, " end stage ", i2
		ii = ii+1
		i2 = 2**(ii-1)
		call MPI_Barrier(MPI_COMM_WORLD, ierror)
	enddo
	
end subroutine transposeMPI

subroutine boundary(myid, numproc, slab, slablength, blocklength)
	implicit none
	integer myid, numproc, slablength, blocklength
	real*8 slab(slablength+2, blocklength+2)

	! set boundary conditions
	slab(1,:) = 1d0
!	slab(slablength+2,:) = 1d0
!	if (myid .eq. 0) then
!		slab(:,1) = 1d0
!	endif
!	if (myid .eq. numproc-1) then
!		slab(:,blocklength+2) = 1d0
!	endif
end subroutine boundary

! this is the transpose of the other one
subroutine boundary2(myid, numproc, slab, slablength, blocklength)
	implicit none
	integer myid, numproc, slablength, blocklength
	real*8 slab(slablength+2, blocklength+2)

	! set boundary conditions
!	slab(:,1) = 1d0
!	slab(:,slablength+2) = 1d0
	if (myid .eq. 0) then
		slab(:,1) = 1d0
	endif
!	if (myid .eq. numproc-1) then
!		slab(blocklength+2,:) = 1d0
!	endif
end subroutine boundary2

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
