all : main.out

clean :
	rm main.out
main.out : main.f90 gtd2d_libinter_cfun.a
	gfortran -fno-underscoring -ffree-line-length-none -O3 -o main.out main.f90 gtd2d_libinter_cfun.a -lm

