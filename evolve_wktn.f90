
! wakatani code to compare with transpose method
! evolves a single diff eq on a grid of arbitrary size

program evolve_wktn
	implicit none
	include "mpif.h"

	integer myid, numproc, ierr, stat(MPI_STATUS_SIZE), dt_slab, req
	real*8, dimension(:,:), allocatable :: a, b, c, d, x, bp, acorr, acorr2
	real*8, dimension(:,:), allocatable :: slab1, slab2
	real*8, dimension(:), allocatable :: x2, ghost1, ghost2
	real*8 dx, dt, dt2, t0, t1, t2, avgtime
	integer step, n, indexspan, startindex, ii, ij, ik, ms, slabsize, asdf, maxstep
!	parameter (n=1024)
!	parameter (ms=16) ! message length
	character *100 buffer
	
	call getarg(1,buffer)
	read(buffer,*) n
	call getarg(2,buffer)
	read(buffer,*) maxstep
	call getarg(3,buffer)
	read(buffer,*) ms

	! MPI goodness
	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, numproc, ierr)
	! slab decomposition
	startindex = myid*n/numproc+1
	indexspan = n/numproc
	slabsize = (n+2)*(indexspan+2)
	asdf = n/ms
	allocate(slab1(n+2,indexspan+2)) ! includes boundary layers
	allocate(slab2(n+2,indexspan+2))
	! things for the tridis
	allocate(a(n,n))
	allocate(b(n,n))
	allocate(c(n,n))
	allocate(d(n,n))
	allocate(x(n,indexspan)) ! is this needed?
	allocate(x2(n))
	allocate(ghost1(n))
	allocate(ghost2(n))
	! precomputation things
	allocate(bp(n,n))
	allocate(acorr(n,indexspan))
	allocate(acorr2(n,indexspan))

	call MPI_TYPE_CONTIGUOUS(slabsize, MPI_DOUBLE_PRECISION, dt_slab, ierr)
	call MPI_TYPE_COMMIT(dt_slab, ierr)

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
	do step=1,maxstep
		! start in one direction

		! wakatani call
		do ij=1,n
			do ik=2,indexspan+1
				d(ij,startindex+ik-2) = slab1(ij+2,ik) + 2d0*(dx*dx/dt-1d0)*slab1(ij+1,ik) + slab1(ij,ik)
			enddo
			d(ij,startindex) = d(ij,startindex) + slab1(ij+1,1)
			d(ij,startindex+indexspan-1) = d(ij,startindex+indexspan-1) + slab1(ij+1,indexspan+2)
			! problem probably has to do with d() and whether or not it gets ghost points in it
			! also where they go. does the mpi part use d() outside of a slab?
		enddo
		call tridi_block(a,b,c,d,x,n,myid,numproc,bp,acorr,acorr2,indexspan,startindex,ms,asdf)
		slab2(2:n+1,2:indexspan+1) = x

		call ghosts(myid, numproc, slab2, n, indexspan)

		! this is cache-preferred direction
		do ii=2,indexspan+1
			!i2 = i+startindex-1
			! load up d
			do ik=2,n+1
				d(ik-1,ii) = slab2(ik,ii+1) + 2d0*(dx*dx/dt-1d0)*slab2(ik,ii) + slab2(ik,ii-1)
			enddo
			! deal with boundary conditions
			d(1,ii) = d(1,ii) + slab2(1,ii)
			d(n,ii) = d(n,ii) + slab2(n+2,ii)
			! call tridi
			call solve_tridiag(a(:,ii),b(:,ii),c(:,ii),d(:,ii),x2,n)
			! load result into other matrix
			slab1(2:n+1,ii) = x2
		enddo
	
	enddo
	t1 = MPI_WTIME()
	t2 = t1-t0
	
	call MPI_ALLREDUCE(t2, avgtime, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
	if (myid .eq. 0) then
		write(*,*), avgtime/numproc
	endif

	! output to file for legit comparing
	if (myid .eq. 0) then
		open(25,file="out_wktn")
		write(25,*),n,n
		do ii=2,indexspan+1
			write(25,*),slab1(2:n+1,ii)
		enddo
		do ii=1,numproc-1
			call MPI_RECV(slab1, 1, dt_slab, ii, 5, MPI_COMM_WORLD, stat, ierr)
			do ij=2,indexspan+1
				write(25,*),slab1(2:n+1,ij)
			enddo
		enddo
		close(25)
	else
		call MPI_SEND(slab1, 1, dt_slab, 0, 5, MPI_COMM_WORLD, ierr)
	endif

	deallocate(slab1)
	deallocate(slab2)
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


subroutine ghosts(myid, numproc, slab, n, indexspan)
	implicit none
	include "mpif.h"
	integer myid, numproc, n, indexspan, req, ierr
	integer stat(MPI_STATUS_SIZE)
	real*8 slab(n+2,indexspan+2)
	real*8 sendbuff(n), recvbuff(n)

	if (myid .ne. 0) then
		sendbuff = slab(2:n+1,2)
		call MPI_Isend(sendbuff, n, MPI_REAL8, myid-1, 7, MPI_COMM_WORLD, req, ierr)
	endif
	if (myid .ne. numproc-1) then
		call MPI_Recv(recvbuff, n, MPI_REAL8, myid+1, 7, MPI_COMM_WORLD, stat, ierr)
		slab(2:n+1,indexspan+2) = recvbuff
	endif
	if (myid .ne. 0) then
		call MPI_Wait(req, stat, ierr)
	endif
	if (myid .ne. numproc-1) then
		sendbuff = slab(2:n+1,indexspan+1)
		call MPI_Isend(sendbuff, n, MPI_REAL8, myid+1, 7, MPI_COMM_WORLD, req, ierr)
	endif
	if (myid .ne. 0) then
		call MPI_Recv(recvbuff, n, MPI_REAL8, myid-1, 7, MPI_COMM_WORLD, stat, ierr)
		slab(2:n+1,1) = recvbuff
	endif
	if (myid .ne. numproc-1) then
		call MPI_Wait(req, stat, ierr)
	endif
end subroutine ghosts



! acorr, acorr2, and bp can all be computed once by a subroutine on program start
! for now, assume n mod m = 0
! SO MANY PARAMETERS
subroutine tridi_block(a,b,c,d,x,n,myid,numproc,bp,acorr,acorr2,indexspan,startindex,m,numblocks)
	implicit none
	include "mpif.h"
	integer n, myid, numproc, ii, ij, ierr
	integer stat(MPI_STATUS_SIZE)
	integer indexspan, startindex, m
	integer numblocks, s, e, crow, checkrecv, flag, done, z, span
	integer, dimension(numblocks) :: sendreq, recvreq
	real(8), dimension(n,n) :: a,b,c,d,bp ! apparently needs the ::
	real(8), dimension(n,indexspan) :: vp, x, acorr, acorr2
	real(8), dimension(n) :: recvbuff
	real(8) :: r, corr, sendcorr
	parameter (z=1) ! number of rows to calc while waiting

	sendreq = 0
	recvreq = 0

! Tridi pass 1
	if (myid .eq. 0) then
		! all root needs to do is calc and send off
		do ii=1,numblocks
			s=1+(ii-1)*m
			e=ii*m
!			print*, "root calc from ",s," to ",e
			vp(s:e,1) = d(s:e,startindex)
			do ij=2,indexspan
				vp(s:e,ij) = d(s:e,startindex+ij-1) - a(s:e,startindex+ij-1)*vp(s:e,ij-1)/bp(s:e,startindex+ij-2)
			enddo
!			print*, "root isend block ",ii," to proc 1"
			call MPI_Isend(vp(s:e,indexspan), m, MPI_REAL8, myid+1, 0, MPI_COMM_WORLD, sendreq(ii), ierr)
		enddo
	else
		! post a bunch of recvs for checking later
		do ii=1,numblocks
			s = 1+(ii-1)*m
			e = ii*m
!			print*, myid, "posting irecv of block ",ii," for later"
			call MPI_Irecv(recvbuff(s:e), m, MPI_REAL8, myid-1, 0, MPI_COMM_WORLD, recvreq(ii), ierr)
		enddo

		! cant do anything without the first block
!		print*, myid, "doing first block of guessing"
		vp(1:m,1) = d(1:m,startindex)
		do ij=2,indexspan
			vp(1:m,ij) = d(1:m,startindex+ij-1) - a(1:m,startindex+ij-1)*vp(1:m,ij-1)/bp(1:m,startindex+ij-2)
		enddo
		crow = 1
		checkrecv = 1
		! all other procs are waiting on recvs, so work on things
		! assume blocks will come in order
		done = 0
		do while (done .eq. 0)
			
			! check for recv
			call MPI_Test(recvreq(checkrecv), flag, stat, ierr)
			! if recv'd, look for next one
			if (flag .eq. 1) then
!				print*, myid, "recvd block ", checkrecv
				! quickly correct final row and send it off
				s = 1+(checkrecv-1)*m
				e = checkrecv*m
				vp(s:e,indexspan) = vp(s:e,indexspan) + acorr(s:e,indexspan)*recvbuff(s:e)
				if (myid .ne. numproc-1) then
!					print*, myid, "propagating to next proc"
					call MPI_Isend(vp(s:e,indexspan), m, MPI_REAL8, myid+1, 0, MPI_COMM_WORLD, sendreq(checkrecv), ierr)
				endif
				checkrecv = checkrecv + 1
				! if recvd all, exit while loop
				if (checkrecv .gt. numblocks) then
!					print*, myid, "all blocks recvd"
					done = 1
				else
					! still have more blocks to come
					! make sure guess is prepped
					s = e+1 ! next block
					e = e+m
					vp(s:e,1) = d(s:e,startindex)
					do ij=2,indexspan
						vp(s:e,ij) = d(s:e,startindex+ij-1) - a(s:e,startindex+ij-1)*vp(s:e,ij-1)/bp(s:e,startindex+ij-2)
					enddo
				endif
			else

			! correct a few rows
			if (crow .le. n) then
				span = z
!				print*,myid,"waste time?",span,crow,n,checkrecv
				if (crow+span-1 .gt. n) then
					span = n-crow+1
				endif
				! dont correct if havent recvd a block yet
				if (crow+span-1 .gt. m*(checkrecv-1)) then
					span = m*(checkrecv-1)-crow+1
!					print*,myid,"trunc2 to",span
				endif
				! correct from row (crow : crow+span-1)
				if (span .gt. 0) then
					s = crow
					e = crow+span-1
!					print*, myid, "wasting time correcting from ",s," to ",e
					do ij=1,indexspan-1
						vp(s:e,ij) = vp(s:e,ij) + acorr(s:e,ij)*recvbuff(s:e)
					enddo
					crow = crow + span
				endif
			endif

			endif

		end do
		! correct remaining things
		if (crow .le. n) then
			s = crow
			e = n
!			print*, myid, "computing rest of corrections from ",s," to ",e
			do ij=1,indexspan-1
				vp(s:e,ij) = vp(s:e,ij) + acorr(s:e,ij)*recvbuff(s:e)
			enddo
		endif
	endif

	! reset recv req array
	recvreq = 0
	! wait for all sends to complete?
	! not necessary if clever
	do ij=1,numblocks ! mpi_waitall looks annoying
		call MPI_Wait(sendreq(ij), stat, ierr)
	enddo
	sendreq = 0

! Start Second Phase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (myid .eq. numproc-1) then
		! all root needs to do is calc and send off
		do ii=1,numblocks
			s=1+(ii-1)*m
			e=ii*m
!			print*, "PHASE2 root calc from ",s," to ",e
			x(s:e,indexspan) = vp(s:e,indexspan)/bp(s:e,indexspan+startindex-1)
			do ij=indexspan-1,1,-1
				x(s:e,ij) = vp(s:e,ij)/bp(s:e,startindex+ij-1) - (c(s:e,startindex+ij-1)/bp(s:e,startindex+ij-1))*x(s:e,ij+1)
			enddo
!			print*, "PHASE2 root isend block ",ii," to proc n-1"
			call MPI_Isend(x(s:e,1), m, MPI_REAL8, myid-1, 0, MPI_COMM_WORLD, sendreq(ii), ierr)
		enddo
	else
		! post a bunch of recvs for checking later
		do ii=1,numblocks
			s = 1+(ii-1)*m
			e = ii*m
!			print*, myid, "PHASE2 posting irecv of block ",ii," for later from",s,"to",e
			call MPI_Irecv(recvbuff(s:e), m, MPI_REAL8, myid+1, 0, MPI_COMM_WORLD, recvreq(ii), ierr)
		enddo

		! cant do anything without the first block
!		print*, myid, "PHASE2 doing first block of guessing"
		x(1:m,indexspan) = vp(1:m,indexspan)/bp(1:m,indexspan+startindex-1)
		do ij=indexspan-1,1,-1
			x(1:m,ij) = vp(1:m,ij)/bp(1:m,startindex+ij-1) - (c(1:m,startindex+ij-1)/bp(1:m,startindex+ij-1))*x(1:m,ij+1)
		enddo
		crow = 1
		checkrecv = 1
		! all other procs are waiting on recvs, so work on things
		! assume blocks will come in order
		done = 0
		do while (done .eq. 0)
			
			! check for recv
			call MPI_Test(recvreq(checkrecv), flag, stat, ierr)
			! if recv'd, look for next one
			if (flag .eq. 1) then
!				print*, myid, "PHASE2 recvd block ", checkrecv
				! quickly correct final row and send it off
				s = 1+(checkrecv-1)*m
				e = checkrecv*m
				x(s:e,1) = x(s:e,1) + acorr2(s:e,1)*recvbuff(s:e)
				if (myid .ne. 0) then
!					print*, myid, "PHASE2 propagating to next proc"
					call MPI_Isend(x(s:e,1), m, MPI_REAL8, myid-1, 0, MPI_COMM_WORLD, sendreq(checkrecv), ierr)
				endif
				checkrecv = checkrecv + 1
				! if recvd all, exit while loop
				if (checkrecv .gt. numblocks) then
!					print*, myid, "PHASE2 all blocks recvd"
					done = 1
				else
					! still have more blocks to come
					! make sure guess is prepped
					s = e+1 ! next block
					e = e+m
					x(s:e,indexspan) = vp(s:e,indexspan)/bp(s:e,indexspan+startindex-1)
					do ij=indexspan-1,1,-1
						x(s:e,ij) = vp(s:e,ij)/bp(s:e,startindex+ij-1) - (c(s:e,startindex+ij-1)/bp(s:e,startindex+ij-1))*x(s:e,ij+1)
					enddo
				endif
			else

			! correct a few rows
			if (crow .le. n) then
				span = z
!				print*,myid,"waste time?",span,crow,n,checkrecv
				if (crow+span-1 .gt. n) then
					span = n-crow+1
				endif
				! dont correct if havent recvd a block yet
				if (crow+span-1 .gt. m*(checkrecv-1)) then
					span = m*(checkrecv-1)-crow+1
!					print*,myid,"trunc2 to",span
				endif
				! correct from row (crow : crow+span-1)
				if (span .gt. 0) then
					s = crow
					e = crow+span-1
!					print*, myid, "PHASE2 wasting time correcting from ",s," to ",e
					do ij=2,indexspan
						x(s:e,ij) = x(s:e,ij) + acorr2(s:e,ij)*recvbuff(s:e)
					enddo
					crow = crow + span
				endif
			endif

			endif

		end do
		! correct remaining things
		if (crow .le. n) then
			s = crow
			e = n
!			print*, myid, "PHASE2 computing rest of corrections from ",s," to ",e
			do ij=2,indexspan
				x(s:e,ij) = x(s:e,ij) + acorr2(s:e,ij)*recvbuff(s:e)
			enddo
		endif
	endif

	! reset recv req array
	recvreq = 0
	! wait for all sends to complete?
	! not necessary if clever
	do ij=1,numblocks ! mpi_waitall looks annoying
		call MPI_Wait(sendreq(ij), stat, ierr)
	enddo

	! DOOONE
end subroutine tridi_block


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
