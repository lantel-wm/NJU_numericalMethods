jsff6: jsff6_1.f90
	gfortran -o jsff6 jsff6.f90 jsff6_1.f90 jsff6_2.f90 functions.f90
clean:
	rm jsff6
run:
	gfortran -o jsff6 jsff6.f90 jsff6_1.f90 jsff6_2.f90 functions.f90 && ./jsff6
gdb:
	gfortran -o jsff6 -g jsff6.f90 jsff6_1.f90 jsff6_2.f90 functions.f90 && gdb jsff6
