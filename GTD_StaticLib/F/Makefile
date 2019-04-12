all : main.out

clean :
	rm main.out
main.out : main.f90 gtd_eig_cfun.a
	gfortran -fno-underscoring -ffree-line-length-none -O3 -o main.out main.f90 gtd_eig_cfun.a -lm

