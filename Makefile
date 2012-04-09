

serial : serial.f90
	gfortran -o serial serial.f90

transpose : transpose.f90
	mpif90 -o transpose transpose.f90

simple_wktn : simple_wktn.f90
	mpif90 -o simple_wktn simple_wktn.f90

tridi_wktn : tridi_wktn.f90
	mpif90 -o tridi_wktn tridi_wktn.f90
