jsff7:jsff7.f90
	python ./xls2arr.py
	gfortran -o jsff7 -fbackslash jsff7.f90 subroutines.f90
	git add .
	git commit -m 'make'
clean:
	rm jsff7 wind.png div_wind.png vor_wind.txt u.txt v.txt phi.txt up.txt vp.txt
run:
	python ./xls2arr.py
	gfortran -o jsff7 -fbackslash jsff7.f90 subroutines.f90 
	./jsff7
	python ./plot_wind.py
	git add .
	git commit --allow-empty -m 'make run' 
gdb:
	python ./xls2arr.py
	gfortran -o jsff7 -g -fbackslash jsff7.f90 subroutines.f90 
	git add .
	git commit -m 'make gdb' 
	gdb jsff7


