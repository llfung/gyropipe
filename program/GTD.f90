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

    !DOUBLE PRECISION :: GTD_G(0:3), GTD_G_(0:3), Diff(0:4)
    INTEGER :: insize(1), outsize(1)
  contains
!------------------------------------------------------------------------
!  Algorithm for computing D_T from G=Grad(u) at BC for no-flux condition of H
!------------------------------------------------------------------------
    subroutine GTD_compute_bc()
      DOUBLE PRECISION :: loc_G(i_pZ*i_Th,0:3)
      DOUBLE PRECISION :: loc_Drr(i_pZ*i_Th,1),loc_Drz(i_pZ*i_Th,1),loc_Dzz(i_pZ*i_Th,1)
      DOUBLE PRECISION :: loc_er(i_pZ*i_Th,1),loc_ez(i_pZ*i_Th,1)
      if (mpi_rnk/=(_Nr-1)) return
      loc_G(:,:)=0d0

      insize(1)=i_pZ*i_Th

      call gtd2d_libinter_cfunvec(vel_Grr%Re(:,:,mes_D%pN)/d_dr,insize, &
      loc_G(:,1),insize, &
      vel_Grz%Re(:,:,mes_D%pN)/d_dr,insize, &
      loc_G(:,3),insize, &
      loc_Drr,outsize, &
      loc_Drz,outsize, &
      loc_Dzz,outsize, &
      loc_er,outsize, &
      loc_ez,outsize)

      GTD_Drr_bc%Re=RESHAPE(loc_Drr,(/i_pZ,i_Th/))
      GTD_Drt_bc%Re=0d0
      GTD_Drz_bc%Re=RESHAPE(loc_Drz,(/i_pZ,i_Th/))
      GTD_er_bc %Re=RESHAPE(loc_er,(/i_pZ,i_Th/))
    end subroutine GTD_compute_bc
!------------------------------------------------------------------------
!  Main Algorithm for computing D_T from G=Grad(u)
!------------------------------------------------------------------------
    subroutine GTD_compute()
      DOUBLE PRECISION :: loc_Drr(mes_D%pN*i_pZ*i_Th,1,1),loc_Drz(mes_D%pN*i_pZ*i_Th,1,1),loc_Dzz(mes_D%pN*i_pZ*i_Th,1,1)
      DOUBLE PRECISION :: loc_er(mes_D%pN*i_pZ*i_Th,1,1),loc_ez(mes_D%pN*i_pZ*i_Th,1,1)
      ! if ((MAXVAL(dabs(loc_G(:,0)))>0.25d0 .or. MAXVAL(dabs(loc_G(:,1)))>0.25d0 .or. MAXVAL(dabs(loc_G(:,3)))>0.25d0 .or. MINVAL(loc_G(:,2))<-2.8d0 .or. MAXVAL(loc_G(:,2))>1.4d0) .and. .NOT.(EXTRAPOLATE_FLAG)) then
      !   print*,' Extrapolating GTD!'
      !   EXTRAPOLATE_FLAG=.TRUE.
      ! end if

      insize(1)=mes_D%pN*i_pZ*i_Th

      call gtd2d_libinter_cfunvec(vel_Grr%Re/d_dr,insize, &
      vel_Gzr%Re/d_dr,insize, &
      vel_Grz%Re/d_dr,insize, &
      vel_Gzz%Re/d_dr,insize, &
      loc_Drr,outsize, &
      loc_Drz,outsize, &
      loc_Dzz,outsize, &
      loc_er,outsize, &
      loc_ez,outsize)

      GTD_Drr%Re=RESHAPE(loc_Drr,(/i_pZ,i_Th,mes_D%pN/))
      GTD_Drz%Re=RESHAPE(loc_Drz,(/i_pZ,i_Th,mes_D%pN/))
      GTD_Dzz%Re=RESHAPE(loc_Dzz,(/i_pZ,i_Th,mes_D%pN/))
      GTD_er %Re=RESHAPE( loc_er,(/i_pZ,i_Th,mes_D%pN/))
      GTD_ez %Re=RESHAPE( loc_ez,(/i_pZ,i_Th,mes_D%pN/))
      GTD_Drt%Re=0d0
      GTD_Dtt%Re=0d0
      GTD_Dtz%Re=0d0
      GTD_et %Re=0d0

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
    end subroutine GTD_compute
!------------------------------------------------------------------------
!  Initialise stuff for the module and the MATLAB Coder generated library
!------------------------------------------------------------------------
    subroutine GTD_precompute()
	    EXTRAPOLATE_FLAG=.FALSE.
      call gtd2d_libinter_cfunvec_initialize()
    end subroutine GTD_precompute
!------------------------------------------------------------------------
!  Decollocate memory stuff for the MATLAB Coder generated library
!------------------------------------------------------------------------
    subroutine GTD_closing()
      call gtd2d_libinter_cfunvec_terminate()
    end subroutine GTD_closing
  end module GTD
