
program compare
   ! This sample program creates two arrays, writes each array out to a file,
   ! reads both files back into arrays, and compares the two arrays element
   ! by element, looking for a difference greater than 0.01d0.
   integer rows, columns
   double precision, allocatable :: expected_array(:,:), actual_array(:,:)
   rows = 8
   columns = 8
   allocate(expected_array(rows, columns))
   allocate(actual_array(rows, columns))

   ! fill both arrays with the same thing
   do i=1,rows
     do j=1,columns
       expected_array(i,j) = i*10+j
       actual_array(i,j) = i*10+j
     end do
   end do

   ! Write the two arrays to two separate files
   call write_array_to_file(expected_array, rows, columns, 'expected_sample')	
   call write_array_to_file(actual_array, rows, columns, 'actual_sample')	

   ! This function reads the arrays from both files back into memory
   ! and looks for a difference between the elements of greater than 0.01
   call compare_arrays_in_files('expected_sample', 'actual_sample', 0.01d0) 

   ! Now, perturb the actual array by a little bit, write it out, and
   ! compare it against the expected array.  The following call to
   ! compare_arrays_in_files should should write an error message to the screen.
   actual_array(4,4) = actual_array(4,4) + 0.5;
   call write_array_to_file(actual_array, rows, columns, 'actual_sample')	
   call compare_arrays_in_files('expected_sample', 'actual_sample', 0.01d0) 

   deallocate(expected_array)
   deallocate(actual_array)
end program

subroutine write_array_to_file(array, rows, columns, filename)
	implicit none
	integer,intent(in) :: rows, columns
    real(8),dimension(rows, columns),intent(in) :: array
    character(len=*),intent(in) :: filename
    integer i,j

    open(30,file=filename,status='unknown')
    write (30,*) rows, columns
    do i=1,rows
      write (30,*) (array(i,j), j=1,columns)
    end do
    close(30)
end subroutine write_array_to_file

subroutine compare_arrays_in_files(expected_filename, actual_filename, tolerance)
    implicit none
    character(len=*),intent(in) :: expected_filename, actual_filename
    real(8),intent(in) :: tolerance
    integer i, j, expected_rows, expected_columns, actual_rows, actual_columns
    integer arrays_differ
    double precision, allocatable :: expected_array(:,:), actual_array(:,:)

    open(30,file=expected_filename,status='unknown')
    ! By convention, the first line in the file has the number of rows and
    ! the number of columns.  Read in the number of rows and columns
    ! so we know how much memory to allocate.
    read(30,*) expected_rows, expected_columns
    allocate(expected_array(expected_rows, expected_columns))
    
    do i=1,expected_rows
      !do j=1,expected_columns
        !expected_array(i,j) = 0
      !end do
      read (30,*) (expected_array(i,j), j=1,expected_columns)
    end do
    close(30)

    open(30,file=actual_filename,status='unknown')
    read(30,*) actual_rows, actual_columns
    allocate(actual_array(actual_rows, actual_columns))
    do i=1,actual_rows
      !do j=1,actual_columns
        !actual_array(i,j) = 0
      !end do
      read (30,*) (actual_array(i,j), j=1,actual_columns)
    end do
    close(30)

    arrays_differ = 0
    do i=1,expected_rows
      do j=1,expected_columns
        if (abs((expected_array(i,j) - actual_array(i,j))/expected_array(i,j)) > abs(tolerance)) then
          print *, 'ERROR:  Arrays differ by more than the tolerance at array location', i, ',', j
          print *, 'Expected value: ', expected_array(i,j), ' Actual value: ', actual_array(i,j)
          arrays_differ = 1
        endif
      end do
    end do

    if(arrays_differ .eq. 0) then
      print *, 'SUCCESS:  the arrays in files ', expected_filename, ' and ', actual_filename, ' are within tolerance of ', tolerance
    endif

    deallocate(expected_array)
    deallocate(actual_array)
end subroutine compare_arrays_in_files
