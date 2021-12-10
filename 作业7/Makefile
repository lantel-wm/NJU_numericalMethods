jsff7:jsff7.f90
	python ./xls2arr.py
	gfortran -o jsff7 -fbackslash jsff7.f90 jsff7_1.f90 subroutines.f90
	git commit -am 'make'
clean:
	rm jsff7
run:
	python ./xls2arr.py
	gfortran -o jsff7 -fbackslash jsff7.f90 jsff7_1.f90 subroutines.f90 
	git commit -am 'make run' 
	./jsff7
gdb:
	python ./xls2arr.py
	gfortran -o jsff7 -g -fbackslash jsff7.f90 jsff7_1.f90 subroutines.f90 
	git commit -am 'make gdb' 
	gdb jsff7


