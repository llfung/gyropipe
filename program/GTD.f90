!*************************************************************************
! GTD Table Lookup and interpolation
!*************************************************************************
#include "../parallel.h"
  module GTD
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

    INTEGER        :: II,JJ,KK
    LOGICAL        :: EXTRAPOLATE_FLAG

    DOUBLE PRECISION :: GTD_G(0:2), GTD_G_(0:2), Diff(0:4)
  contains
!------------------------------------------------------------------------
!  Algorithm for computing D_T from G=Grad(u) at BC for no-flux condition of H
!------------------------------------------------------------------------
    subroutine GTD_compute_bc()
      if (mpi_rnk/=(_Nr-1)) return
      GTD_G(:)=0d0
      KK=mes_D%pN
        do JJ=0,i_Th-1
          do II=0,i_pZ-1
            GTD_G_=GTD_G
            GTD_G(2)=vel_Grz%Re(II,JJ,KK)/d_dr

            ! call gtd_eig_cfun(GTD_G,Diff)
            if (maxval(dabs(GTD_G_-GTD_G))/=0d0) call gtd2d_libinter_cfun(GTD_G,Diff)

            GTD_Drr_bc%Re(II,JJ)=Diff(0)
            GTD_Drt_bc%Re(II,JJ)=0d0
            GTD_Drz_bc%Re(II,JJ)=Diff(1)
            GTD_er_bc%Re(II,JJ)=Diff(3)
            end do
        end do
    end subroutine GTD_compute_bc
!------------------------------------------------------------------------
!  Main Algorithm for computing D_T from G=Grad(u)
!------------------------------------------------------------------------
    subroutine GTD_compute()
      GTD_G_=0d0
      do KK=1,mes_D%pN
        do II=0,i_pZ-1
          do JJ=0,i_Th-1
            GTD_G_=GTD_G
            GTD_G(0)=vel_Grr%Re(II,JJ,KK)/d_dr
            GTD_G(1)=vel_Gzr%Re(II,JJ,KK)/d_dr
            GTD_G(2)=vel_Grz%Re(II,JJ,KK)/d_dr
            if ((dabs(GTD_G(0))>0.25d0 .or. dabs(GTD_G(1))>0.25d0 .or. (GTD_G(2))<-2.8d0 .or. (GTD_G(2))>1.4d0) .and. .NOT.(EXTRAPOLATE_FLAG)) then
              print*,' Extrapolating GTD!', GTD_G
              EXTRAPOLATE_FLAG=.TRUE.
            end if
            if (maxval(dabs(GTD_G_-GTD_G))/=0d0) then
              call gtd2d_libinter_cfun(GTD_G,Diff)
            end if
            GTD_Drr%Re(II,JJ,KK)=Diff(0)
            GTD_Drt%Re(II,JJ,KK)=0d0
            GTD_Drz%Re(II,JJ,KK)=Diff(1)
            GTD_Dtt%Re(II,JJ,KK)=0d0
            GTD_Dtz%Re(II,JJ,KK)=0d0
            GTD_Dzz%Re(II,JJ,KK)=Diff(2)
            GTD_er%Re(II,JJ,KK)=Diff(3)
            GTD_et%Re(II,JJ,KK)=0d0
            GTD_ez%Re(II,JJ,KK)=Diff(4)
          end do
        end do
      end do
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
    end subroutine GTD_compute
!------------------------------------------------------------------------
!  Initialise stuff for the module and the MATLAB Coder generated library
!------------------------------------------------------------------------
    subroutine GTD_precompute()
	    EXTRAPOLATE_FLAG=.FALSE.
      call gtd2d_libinter_cfun_initialize()
    end subroutine GTD_precompute
!------------------------------------------------------------------------
!  Decollocate memory stuff for the MATLAB Coder generated library
!------------------------------------------------------------------------
    subroutine GTD_closing()
      call gtd2d_libinter_cfun_terminate()
    end subroutine GTD_closing
  end module GTD
