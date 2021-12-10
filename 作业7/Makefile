jsff7:jsff7.f90
	gfortran -o jsff7 -fbackslash jsff7.f90 jsff7_1.f90 subroutines.f90 && git commit -a -m 'make'
clean:
	rm jsff7
run:
	gfortran -o jsff7 -fbackslash jsff7.f90 jsff7_1.f90 subroutines.f90 && git commit -a -m 'make' && ./jsff7
gdb:
	gfortran -o jsff7 -g -fbackslash jsff7.f90 jsff7_1.f90 subroutines.f90 && git commit -a -m 'make' && gdb jsff7


