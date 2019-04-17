  PROGRAM MAIN
    IMPLICIT NONE
    DOUBLE PRECISION :: Gin(0:2), Diff(0:5), t_ini,t_fin
    INTEGER :: I
    Gin(0)=1d0
  	Gin(1)=-1d0
  	Gin(2)=-0.5d0


	! call CPU_TIME(t_ini)
     call gtd2d_libinter_cfun_initialize()
	! call CPU_TIME(t_fin)
	! print*, t_fin-t_ini

  call CPU_TIME(t_ini)
    do I=1,576
      call gtd2d_libinter_cfun(Gin,Diff)
    end do
  call CPU_TIME(t_fin)

	print*, t_fin-t_ini
	! call CPU_TIME(t_ini)
     call gtd2d_libinter_cfun_terminate()
	! call CPU_TIME(t_fin)
	! print*, t_fin-t_ini
	print*, Diff

  END PROGRAM MAIN
