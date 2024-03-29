
! trying for the message vectorization version

program tridi_wktn_block
	implicit none
	include "mpif.h"
	
	integer n, myid, numproc, ierr, stat(MPI_STATUS_SIZE)
	real*8, dimension(:,:), allocatable :: a, b, c, d, x, bp, acorr, acorr2
	integer indexspan, startindex, ii, ij, m
	parameter (n=8)
	parameter (m=2) ! message length

	real*8 t0, t1, t2, avgtime

	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, numproc, ierr)

	indexspan = n/numproc
	startindex = myid*indexspan+1

	allocate(a(n,n))
	allocate(b(n,n))
	allocate(c(n,n))
	allocate(d(n,n))
	allocate(x(n,indexspan))
	allocate(bp(n,n))
	allocate(acorr(n,indexspan))
	allocate(acorr2(n,indexspan))

	! initialize things
	a = 1d0
	b = -2d0
	c = 1d0
	d = -5d-1

	if (numproc .gt. 1) then
		! we want to precompute bp
		bp(:,1) = b(:,1)
		do ii=2,n
			bp(:,ii) = b(:,ii) - (a(:,ii)/bp(:,ii-1))*c(:,ii-1)
		enddo
		
		! we also want to pre-compute the correction coefs for the entire vector
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
	endif

	t0 = MPI_Wtime()
	! for testing on 1 processor
	if (numproc .eq. 1) then
		do ii=1,n
			call solve_tridiag(a(ii,:),b(ii,:),c(ii,:),d(ii,:),x(ii,:),n)
		enddo
	else
		call tridi_block(a,b,c,d,x,n,myid,numproc,bp,acorr,acorr2,indexspan,startindex,m,n/m)
	endif
	t1 = MPI_Wtime(ierr) - t0

	! output
	do ii=0,numproc-1
		if (myid .eq. ii) then
			write(*,*), "proc", myid
			do ij=1,indexspan
				write(*,'(8f8.2)'), x(:,ij)
			enddo
		endif
		call MPI_Barrier(MPI_COMM_WORLD, ierr)
	enddo

	call MPI_ALLREDUCE(t1, avgtime, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
	if (myid .eq. 0) then
		write(*,*) (avgtime/numproc)
	endif

	! too lazy to deallocate
	call MPI_Finalize(ierr)

end program tridi_wktn_block


! the magic happens here
! acorr, acorr2, and bp can all be computed once by a subroutine on program start
! for now, assume n mod m = 0
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
!			print*, "PHASE2 root isend block ",ii," to proc 1"
			call MPI_Isend(x(s:e,1), m, MPI_REAL8, myid-1, 0, MPI_COMM_WORLD, sendreq(ii), ierr)
		enddo
	else
		! post a bunch of recvs for checking later
		do ii=1,numblocks
			s = 1+(ii-1)*m
			e = ii*m
!			print*, myid, "PHASE2 posting irecv of block ",ii," for later"
			call MPI_Irecv(recvbuff(s:e), m, MPI_REAL8, myid+1, 0, MPI_COMM_WORLD, recvreq(ii), ierr)
		enddo

		! cant do anything without the first block
!		print*, myid, "PHASE 2doing first block of guessing"
		x(1:m,indexspan) = vp(1:m,indexspan)/bp(1:m,indexspan+startindex-1)
		do ij=indexspan-1,1,-1
			x(1:m,ij) = vp(1:m,ij)/bp(1:m,startindex+ij-1) - (c(1:m,startindex+ij-1)/bp(1:m,startindex+ij-1))*x(1:m,ij+1)
!			x(1:m,ij) = x(1:)
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
!						x(s:e,ij) = x(s:e,ij+1)
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
