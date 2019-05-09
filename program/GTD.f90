!*************************************************************************
! GTD Table Lookup and interpolation
!*************************************************************************
#include "../parallel.h"
  module GTD
    use netcdf
    use velocity
    use variables
    implicit none
    save

    type (phys)    :: GTD_Drr, GTD_Drt, GTD_Drz
    type (phys)    :: GTD_Dtt, GTD_Dtz, GTD_Dzz

    type (phys_bc) :: GTD_Drr_bc, GTD_Drt_bc, GTD_Drz_bc
    type (phys_bc) :: GTD_er_bc

    type (phys)    :: GTD_er, GTD_et, GTD_ez

    type (coll)    :: GTD_er_col, GTD_et_col, GTD_ez_col, GTD_grade, GTD_lapH ! See nonlinear>non_temperature

    INTEGER        :: I,J,K
    LOGICAL        :: EXTRAPOLATE_FLAG

    ! DOUBLE PRECISION :: A(0:7), A_(0:7), Diff(0:8)
    DOUBLE PRECISION :: A(0:2), A_(0:2), Diff(0:4)
  contains
    subroutine GTD_compute_bc()
      if (mpi_rnk/=(_Nr-1)) return
      A(:)=0d0
      K=mes_D%pN
        do J=0,i_Th-1
          do I=0,i_pZ-1
            A_=A
            ! A(3)=vel_Grt%Re(I,J,K)/d_dr
            ! A(6)=vel_Grz%Re(I,J,K)/d_dr
            A(2)=vel_Grz%Re(I,J,K)/d_dr

            ! call gtd_eig_cfun(A,Diff)
            if (maxval(dabs(A_-A))/=0d0) call gtd2d_libinter_cfun(A,Diff)

            ! GTD_Drr_bc%Re(I,J)=Diff(0)
            ! GTD_Drt_bc%Re(I,J)=Diff(1)
            ! GTD_Drz_bc%Re(I,J)=Diff(2)
            !
            ! GTD_er_bc%Re(I,J)=Diff(6)

            GTD_Drr_bc%Re(I,J)=Diff(0)
            GTD_Drt_bc%Re(I,J)=0d0
            GTD_Drz_bc%Re(I,J)=Diff(1)

            GTD_er_bc%Re(I,J)=Diff(3)
            end do
        end do

    end subroutine GTD_compute_bc
    ! Main Algorithm
    subroutine GTD_compute()
      ! Calculate
      A_=0d0
      do K=1,mes_D%pN
        do I=0,i_pZ-1
          do J=0,i_Th-1
            A_=A
            ! A(0)=vel_Grr%Re(I,J,K)/d_dr
            ! A(1)=vel_Gtr%Re(I,J,K)/d_dr
            ! A(2)=vel_Gzr%Re(I,J,K)/d_dr
            ! A(3)=vel_Grt%Re(I,J,K)/d_dr
            ! A(4)=vel_Gtt%Re(I,J,K)/d_dr
            ! A(5)=vel_Gzt%Re(I,J,K)/d_dr
            ! A(6)=vel_Grz%Re(I,J,K)/d_dr
            ! A(7)=vel_Gtz%Re(I,J,K)/d_dr
            A(0)=vel_Grr%Re(I,J,K)/d_dr
            A(1)=vel_Gzr%Re(I,J,K)/d_dr
            A(2)=vel_Grz%Re(I,J,K)/d_dr
            if ((dabs(A(0))>2d0 .or. dabs(A(1))>2d0 .or. dabs(A(2))>5d0) .and. .NOT.(EXTRAPOLATE_FLAG)) then
		print*,' Extrapolating GTD!', A
		EXTRAPOLATE_FLAG=.TRUE.
            end if
	    if (maxval(dabs(A_-A))/=0d0) then
              ! call gtd_eig_cfun(A,Diff)
              call gtd2d_libinter_cfun(A,Diff)
            end if
            ! GTD_Drr%Re(I,J,K)=Diff(0)
            ! GTD_Drt%Re(I,J,K)=Diff(1)
            ! GTD_Drz%Re(I,J,K)=Diff(2)
            ! GTD_Dtt%Re(I,J,K)=Diff(3)
            ! GTD_Dtz%Re(I,J,K)=Diff(4)
            ! GTD_Dzz%Re(I,J,K)=Diff(5)
            ! GTD_er%Re(I,J,K)=Diff(6)
            ! GTD_et%Re(I,J,K)=Diff(7)
            ! GTD_ez%Re(I,J,K)=Diff(8)
            GTD_Drr%Re(I,J,K)=Diff(0)
            GTD_Drt%Re(I,J,K)=0d0
            GTD_Drz%Re(I,J,K)=Diff(1)
            GTD_Dtt%Re(I,J,K)=0d0
            GTD_Dtz%Re(I,J,K)=0d0
            GTD_Dzz%Re(I,J,K)=Diff(2)
            GTD_er%Re(I,J,K)=Diff(3)
            GTD_et%Re(I,J,K)=0d0
            GTD_ez%Re(I,J,K)=Diff(4)
          end do
        end do
      end do
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
    end subroutine GTD_compute

    function InvMat(MatIn)
      DOUBLE COMPLEX :: MatIn(3,3)
      DOUBLE COMPLEX :: InvMat(3,3)
      DOUBLE COMPLEX :: Matdet
      Matdet=MatIn(1,1)*MatIn(2,2)*MatIn(3,3) &
              +MatIn(1,2)*MatIn(2,3)*MatIn(3,1) &
              +MatIn(1,3)*MatIn(2,1)*MatIn(3,2) &
              -MatIn(1,1)*MatIn(2,3)*MatIn(3,2) &
              -MatIn(1,2)*MatIn(2,1)*MatIn(3,3) &
              -MatIn(1,3)*MatIn(2,2)*MatIn(3,1)
      if (abs(Matdet)<1d-8) print*, 'InvMat of W: det(W)=', Matdet
      InvMat(1,1)=MatIn(2,2)*MatIn(3,3)-MatIn(2,3)*MatIn(3,2)
      InvMat(1,2)=-MatIn(1,2)*MatIn(3,3)+MatIn(1,3)*MatIn(3,2)
      InvMat(1,3)=MatIn(1,2)*MatIn(2,3)-MatIn(1,3)*MatIn(2,2)
      InvMat(2,1)=-MatIn(2,1)*MatIn(3,3)+MatIn(2,3)*MatIn(3,1)
      InvMat(2,2)=MatIn(1,1)*MatIn(3,3)-MatIn(1,3)*MatIn(3,1)
      InvMat(2,3)=-MatIn(1,1)*MatIn(2,3)+MatIn(1,3)*MatIn(2,1)
      InvMat(3,1)=MatIn(2,1)*MatIn(3,2)-MatIn(2,2)*MatIn(3,1)
      InvMat(3,2)=-MatIn(1,1)*MatIn(3,2)+MatIn(1,2)*MatIn(3,1)
      InvMat(3,3)=MatIn(2,2)*MatIn(1,1)-MatIn(2,1)*MatIn(1,2)
      InvMat=InvMat/Matdet
      return
    end function
    subroutine GTD_precompute()
      ! call gtd_eig_cfun_initialize()
	    EXTRAPOLATE_FLAG=.FALSE.
      call gtd2d_libinter_cfun_initialize()
    end subroutine GTD_precompute
    subroutine GTD_closing()
      ! call gtd_eig_cfun_terminate()
      call gtd2d_libinter_cfun_terminate()
    end subroutine GTD_closing
  end module GTD
