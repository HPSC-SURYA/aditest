

serial : serial.f90
	mpif90 -O3 -o serial serial.f90

compare	: compare.f90
	gfortran -o compare compare.f90

compare_output : compare_output.f90
	gfortran -o compare_output compare_output.f90

transpose : transpose.f90
	mpif90 -O3 -o transpose transpose.f90

simple_wktn : simple_wktn.f90
	mpif90 -o simple_wktn simple_wktn.f90

tridi_wktn : tridi_wktn.f90
	mpif90 -o tridi_wktn tridi_wktn.f90

tridi_wktn_block : tridi_wktn_block.f90
	mpif90 -o tridi_wktn_block tridi_wktn_block.f90

evolve_wktn : evolve_wktn.f90
	mpif90 -O3 -o evolve_wktn evolve_wktn.f90
