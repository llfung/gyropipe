  PROGRAM MAIN
    IMPLICIT NONE
    DOUBLE PRECISION :: Gin(0:7), Diff(0:8), t_ini,t_fin
    INTEGER :: I
    Gin(0)=1d0
  	Gin(1)=-1d0
  	Gin(2)=-0.5d0
  	Gin(3)=4d0
  	Gin(4)=-3d0
  	Gin(5)=-1.5d0
  	Gin(6)=2d0
  	Gin(7)=1d0


	! call CPU_TIME(t_ini)
     call gtd_eig_cfun_initialize()
	! call CPU_TIME(t_fin)
	! print*, t_fin-t_ini

  call CPU_TIME(t_ini)
    do I=1,576
      call gtd_eig_cfun(Gin,Diff)
    end do
  call CPU_TIME(t_fin)

	print*, t_fin-t_ini
	! call CPU_TIME(t_ini)
     call gtd_eig_cfun_terminate()
	! call CPU_TIME(t_fin)
	! print*, t_fin-t_ini
	print*, Diff

  END PROGRAM MAIN
