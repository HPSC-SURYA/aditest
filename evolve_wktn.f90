
! wakatani code to compare with transpose method
! evolves a single diff eq on a grid of arbitrary size

program evolve_wktn
	implicit none
	include "mpif.h"

	integer myid, numproc, ierr, stat(MPI_STATUS_SIZE), dt_slab, req
	real*8, dimension(:,:), allocatable :: a, b, c, d, x, bp, acorr, acorr2
	real*8, dimension(:,:), allocatable :: slab1, slab2
	real*8, dx, dt, dt2, t0, t1, t2, avgtime
	integer indexspan, startindex, ii, ij, m, slabsize
	parameter (n=8)
	parameter (m=1) ! message length

	! MPI goodness
	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, numproc, ierr)
	! slab decomposition
	startindex = myid*n/numproc+2
	indexspan = n/numproc
	slabsize = (n+2)*(indexspan+2)
	allocate(slab1(n+2,indexspan+2)) ! includes boundary layers
	allocate(slab2(n+2,indexspan+2))
	! things for the tridis
	allocate(a(n,n))
	allocate(b(n,n))
	allocate(c(n,n))
	allocate(d(n,n))
	allocate(x(n,indexspan)) ! is this needed?
	! precomputation things
	allocate(bp(n,n))
	allocate(acorr(n,indexspan))
	allocate(acorr2(n,indexspan))

	! dx might actually be 1/(n+1)
	dx = 1d0/n
	dt = 0.01d0

	! initialize things to 0
	slab1 = 0.0
	slab2 = 0.0
	call boundary(myid, numproc, slab1, n, indexspan)
	call boundary(myid, numproc, slab2, n, indexspan)

	! these never change, precompute
	a = -1d0
	b = 2d0*(dx*dx/dt+1d0)
	c = -1d0
	! precompute bp
	bp(:,1) = b(:,1)
	do ii=2,n
		bp(:,ii) = b(:,ii) - (a(:,ii)/bp(:,ii-1))*c(:,ii-1)
	enddo
	! precompute correction factors
	acorr = 0d0
	if (myid .ne. 0) then
		acorr(:,1) = -a(:,startindex)/bp(:,startindex-1)
		do ii=2,indexspan
			acorr(:,ii) = (-a(:,startindex+ii-1)/bp(:,startindex+ii-2)) * acorr(:,ii-1)
		enddo
	endif
	acorr2 = 0d0
	if (myid .ne. numproc-1) then
		acorr2(:,indexspan) = -c(:,startindex+indexspan-1)/bp(:,startindex+indexspan-1)
		do ii=indexspan-1,1,-1
			acorr2(:,ii) = (-c(:,startindex+ii-1)/bp(:,startindex+ii-1)) * acorr2(:,ii+1)
		enddo
	endif


	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t0 = MPI_WTIME()
	! start time-stepping
	do step=1,1000
		! start in one direction
		! this is cache-preferred direction
		do i=2,indexspan+1
			!i2 = i+startindex-1
			! load up d
			do k=2,n+1
				d(k-1,i) = slab1(k,i+1) + 2d0*(dx*dx/dt-1d0)*slab1(k,i) + slab1(k,i-1)
			enddo
			! deal with boundary conditions
			d(1,i) = d(1,i) + slab1(1,i)
			d(n,i) = d(n,i) + slab1(n+2,i)
			! call tridi
			call solve_tridiag(a(:,i),b(:,i),c(:,i),d(:,i),x(:,i),n)
			! load result into other matrix
			slab2(2:n+1,i) = x(:,i)
		enddo
	
		! communicate ghost points?

		! now do other direction
		! set up d(,)
		do k=2,n+1
			d(k-1) = slab2(j+1,k) + 2d0*(dx*dx/dt-1d0)*slab2(j,k) + slab2(j-1,k)
		enddo
		d(1) = d(1) + slabt(j,1)
		d(n) = d(n) + slabt(j,n+2)
		call tridi_block(a,b,c,d,x,n,myid,numproc,bp,acorr,acorr2,indexspan,startindex,m,n/m)
		slab1(2:n+1,2:indexspan+1) = x
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
	deallocate(slabt)
	deallocate(slabt2)
	call MPI_TYPE_FREE(dt_slab, ierr)
	call MPI_Finalize(ierr)

end program



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
