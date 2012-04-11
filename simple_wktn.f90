
! the problem with the wakatani method is that it is very abstracted from our problem
! you have to take adi, isolate the tridi solve, then isolate a single recursion relation
! (of which there are two), then apply this algorithm

! the purpose of this code is to just code up the recursion relation
! then we can build on this and have something that is useful to us


! w(i) = a(i) * w(i-1) + b(i)

program simple_wktn
	implicit none
	include "mpif.h"

	integer myid, numproc, ierr, req, stat(MPI_STATUS_SIZE)
	! if you forget that (MPI_..) may Zeus strike you down
	real*8, dimension(:), allocatable :: w, a, b
	real*8 corr, acorr, sendcorr
	integer startindex, indexspan, n, ii, ij
	parameter (n=16)

	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, numproc, ierr)

	! determine problem size
	indexspan = n/numproc
	startindex = myid*indexspan+1
	! allocate things
	allocate(w(indexspan))
	allocate(a(indexspan))
	allocate(b(indexspan))
	! pre-compute coefs
	do ii=1,indexspan
		a(ii) = 1.5d0
		b(ii) = 1d0
		w(ii) = 0d0
	enddo
	if (myid .eq. 0) then
		b(1) = 1d0
	endif
	! for correcting the final data point first
	acorr = 1d0
	do ii=1,indexspan
		acorr = acorr*a(ii)
	enddo

	call MPI_Barrier(MPI_COMM_WORLD, ierr)

	! Guess Phase
	w(1) = a(1)*0d0 + b(1) ! just to illustrate
	do ii=2,indexspan
		w(ii) = a(ii)*w(ii-1) + b(ii)
	enddo

	! Propagation Phase
	if (myid .ne. 0) then
		call MPI_Recv(corr, 1, MPI_REAL8, myid-1, 0, MPI_COMM_WORLD, stat, ierr)
		w(indexspan) = w(indexspan) + acorr*corr
	endif
	if (myid .ne. numproc-1) then
		sendcorr = w(indexspan)
		call MPI_Isend(sendcorr, 1, MPI_REAL8, myid+1, 0, MPI_COMM_WORLD, req, ierr)
	endif

	! Correction Phase
	acorr = 1d0
	do ii=1,indexspan-1
		acorr = acorr * a(ii)
		w(ii) = w(ii) + acorr*corr
	enddo

	call MPI_Wait(req, stat, ierr)

	! output
	do ii=0,numproc-1
		if (myid .eq. ii) then
			write(*,*), "proc", myid
			write(*,'(4f8.1)'), w
		endif
		call MPI_Barrier(MPI_COMM_WORLD, ierr)
	enddo

	deallocate(w)
	deallocate(a)
	deallocate(b)

	call MPI_Finalize(ierr)

end program simple_wktn
