

serial : serial.f90
	gfortran -o serial serial.f90

transpose : transpose.f90
	mpif90 -o transpose transpose.f90
