#HHL1DHyperspherical.x: modules_qd.o besselnew.o Bsplines.o  ../bspllib_22/bspline90_22.o
#	gfortran -O -ffixed-line-length-132 Bsplines.o modules_qd.o ../bspllib_22/bspline90_22.o besselnew.o -L/usr/local/lib/ -L/Users/mehtan/Code/ARPACK/ARPACK -larpack_OSX -framework accelerate -lm HHL1DHyperspherical.f -o HHL1DHyperspherical.x


HHL1DHyperspherical.x: HHL1DHyperspherical.f90 modules_qd.f90 besselnew.f90 Bsplines.f90 ./bspline/bspline90_22.f90
	ifort HHL1DHyperspherical.f90 modules_qd.f90 besselnew.f90 Bsplines.f90 ./bspline/bspline90_22.f90 -mkl -L/opt/ARPACK -larpack_Intel -no-wrap-margin 
