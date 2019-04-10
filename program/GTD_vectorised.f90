!*************************************************************************
! GTD Table Lookup and interpolation
!*************************************************************************
#include "../parallel.h"
  module GTD
    use netcdf
    use velocity
    ! use linear_interpolation_module
    use bspline_module
    use bspline_kinds_module, only: wp
    use variables
    implicit none
    save


    type (phys) :: GTD_S
    type (phys) :: GTD_theta
    type (phys) :: GTD_sin,GTD_cos,GTD_coscos,GTD_sinsin,GTD_sincos
    type (phys) :: GTD_eig_imag
    type (phys) :: GTD_TSOURCE,GTD_FSOURCE,GTD_diag_bool

    type (phys) :: vel_G11, vel_G12, vel_G13
    type (phys) :: vel_G21, vel_G22, vel_G23
    type (phys) :: vel_G31, vel_G32, vel_G33

    type (phys) :: GTD_Drr, GTD_Drt, GTD_Drz
    type (phys) :: GTD_Dtt, GTD_Dtz, GTD_Dzz

    type (phys) :: GTD_er, GTD_et, GTD_ez
    type (phys) :: GTD_e1, GTD_e2, GTD_e3

    type (coll) :: GTD_er_col, GTD_et_col, GTD_ez_col, GTD_grade, GTD_lapH

    type (phys) :: GG1, GG2, GG3

    type (phys_cmp) :: W11, W12, W13
    type (phys_cmp) :: W21, W22, W23
    type (phys_cmp) :: W31, W32, W33
    type (phys_cmp) :: W1_norm, W2_norm, W3_norm

    INTEGER           :: I,J,K,LL,II,JJ,KK,RR,QQ,SS,PP
    INTEGER           :: INFO, ipvt(3)

    DOUBLE PRECISION :: BeReint(27),BeImint(27),BBReint(45),BBImint(45)
    DOUBLE PRECISION :: e123(3),ertz(3)
    DOUBLE COMPLEX :: Be(3,3),BB(3,3),W(3,3),Winv(3,3)
    DOUBLE COMPLEX :: BB_trans(3,3),Diff(3,3),Diff2(3,3),transform_matrix(3,3)
    type (phys) :: det
    type (phys_cmp) :: traceG2, halfdetB, phi
    type (phys_cmp) :: beta1, beta2, beta3
    DOUBLE COMPLEX :: beta_temp

    DOUBLE PRECISION :: Y,X(3)
    DOUBLE PRECISION, ALLOCATABLE :: fint_S(:),fint_eig(:),fint_theta(:)
    ! type(linear_interp_3d) :: fint_BeRe(27)
    ! type(linear_interp_3d) :: fint_BeIm(27)
    ! type(linear_interp_3d) :: fint_BBRe(45)
    ! type(linear_interp_3d) :: fint_BBIm(45)
    type(bspline_3d) :: fint_BeRe(27)
    type(bspline_3d) :: fint_BeIm(27)
    type(bspline_3d) :: fint_BBRe(45)
    type(bspline_3d) :: fint_BBIm(45)
    type(bspline_2d) :: fint_e(3)

    INTEGER :: iflag

    INTEGER           :: file_id
    INTEGER           :: S_len, theta_len, eig_len
    DOUBLE PRECISION              :: S_step, theta_step, eig_step
    DOUBLE PRECISION, ALLOCATABLE :: BeRe_mul(:,:,:,:), BeIm_mul(:,:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: BBRe_mul(:,:,:,:), BBIm_mul(:,:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: e1_mul(:,:), e2_mul(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: e3_mul(:,:)
  contains
    ! Main Algorithm
    subroutine GTD_compute()
      ! Calculate e_avg
      ! if (mpi_rnk==0) then
      ! print*, vel_curlr%Re(i_pZ-1,i_Th-1,:)
      ! print*, '--------------------------------------------'
      ! print*, vel_curlt_comb%Re(i_pZ-1,i_Th-1,:)
      ! print*, '--------------------------------------------'
      ! print*, vel_curlz%Re(i_pZ-1,i_Th-1,:)
      ! print*, '--------------------------------------------------------------------'
      ! endif
      GTD_S%Re(:,:,1:mes_D%pN)=dsqrt(vel_curlr%Re(:,:,1:mes_D%pN)**2d0+vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlz%Re(:,:,1:mes_D%pN)**2d0)
      GTD_theta%Re(:,:,1:mes_D%pN)=dacos(-vel_curlz%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))
      GTD_sin%Re(:,:,1:mes_D%pN)=vel_curlt_comb%Re(:,:,1:mes_D%pN)/dsqrt(vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      GTD_cos%Re(:,:,1:mes_D%pN)=vel_curlr%Re(:,:,1:mes_D%pN)/dsqrt(vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      call GTD_eavg()
      !print*,'ok 4', mpi_rnk
      ! if (mpi_rnk==0) print*, vel_Grr%Re(0,2,4), vel_Grt%Re(0,2,4), vel_Grz%Re(0,2,4)
      ! if (mpi_rnk==0) print*, vel_Gtr%Re(0,2,4), vel_Gtt%Re(0,2,4), vel_Gtz%Re(0,2,4)
      ! if (mpi_rnk==0) print*, vel_Gzr%Re(0,2,4), vel_Gzt%Re(0,2,4), vel_Gzz%Re(0,2,4)
      ! if (mpi_rnk==0) print*, vel_Grz%Re(0,2,:)
      ! vel_Gzr%Re=vel_Gzr%Re-vel_Grz%Re

      ! Check diagonalisability
      call GTD_prevar()
      ! if (mpi_rnk==0) then
      ! print*, traceG2%CMP(i_pZ-1,i_Th-1,:)
      ! print*, '--------------------------------------------'
      ! print*, halfdetB%CMP(i_pZ-1,i_Th-1,:)
      ! print*, '--------------------------------------------'
      ! print*, GTD_S%Re(i_pZ-1,i_Th-1,:)
      ! print*, '--------------------------------------------------------------------'
      ! endif
      GTD_diag_bool%Re=MERGE(GTD_TSOURCE%Re,GTD_FSOURCE%Re,ABS(traceG2%CMP)<1d-8 .or. ABS(abs(DBLE(halfdetB%CMP))-1d0)<1d-8 .or. GTD_S%Re<1d-8)
      ! print*, GTD_diag_bool%Re
      ! Compute eigenvalue
        ! Modify for non-diagonalisable
      vel_curlt_comb%Re=vel_curlt_comb%Re-2d0*GTD_diag_bool%Re
      vel_Gzr%Re=vel_Gzr%Re-GTD_diag_bool%Re
      vel_Grz%Re=vel_Grz%Re+GTD_diag_bool%Re

      GTD_S%Re(:,:,1:mes_D%pN)=dsqrt(vel_curlr%Re(:,:,1:mes_D%pN)**2d0+vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlz%Re(:,:,1:mes_D%pN)**2d0)
      GTD_theta%Re(:,:,1:mes_D%pN)=dacos(-vel_curlz%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))
      GTD_sin%Re(:,:,1:mes_D%pN)=vel_curlt_comb%Re(:,:,1:mes_D%pN)/dsqrt(vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      GTD_cos%Re(:,:,1:mes_D%pN)=vel_curlr%Re(:,:,1:mes_D%pN)/dsqrt(vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      GTD_sinsin%Re(:,:,1:mes_D%pN)=vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0/&
      (vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      GTD_coscos%Re(:,:,1:mes_D%pN)=vel_curlr%Re(:,:,1:mes_D%pN)**2d0/&
      (vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      GTD_sincos%Re(:,:,1:mes_D%pN)=vel_curlr%Re(:,:,1:mes_D%pN)*vel_curlt_comb%Re(:,:,1:mes_D%pN)/&
      (vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      !print*,'ok 5', mpi_rnk
      ! if (mpi_rnk==0) then
      !   print*, vel_Grr%Re(1,1,1), vel_Grt%Re(1,1,1), vel_Grz%Re(1,1,1)
      !   print*, vel_Gtr%Re(1,1,1), vel_Gtt%Re(1,1,1), vel_Gtz%Re(1,1,1)
      !   print*, vel_Gzr%Re(1,1,1), vel_Gzt%Re(1,1,1), vel_Gzz%Re(1,1,1)
      !   print*, GTD_S%Re(1,1,1)
      ! end if

      ! Compute eigenvalue
      call GTD_prevar()
      ! if (mpi_rnk==0) then
      ! print*, traceG2%CMP(i_pZ-1,i_Th-1,:)
      ! print*, '--------------------------------------------'
      ! print*, halfdetB%CMP(i_pZ-1,i_Th-1,:)
      ! print*, '--------------------------------------------'
      ! print*, GTD_S%Re(i_pZ-1,i_Th-1,:)
      ! endif
      call GTD_eigsys(0)
      ! if (mpi_rnk==0) then
      !   print*, GTD_Drr%Re(1,1,1), GTD_Drt%Re(1,1,1), GTD_Drz%Re(1,1,1)
      !   print*, GTD_Dtt%Re(1,1,1), GTD_Dtz%Re(1,1,1), GTD_Dzz%Re(1,1,1)
      ! end if

      !print*,'ok 6', mpi_rnk
      ! ReCompute eigenvalue for non-diagonalisable
      vel_curlt_comb%Re=vel_curlt_comb%Re+4d0*GTD_diag_bool%Re
      vel_Gzr%Re=vel_Gzr%Re+2d0*GTD_diag_bool%Re
      vel_Grz%Re=vel_Grz%Re-2d0*GTD_diag_bool%Re

      GTD_S%Re(:,:,1:mes_D%pN)=dsqrt(vel_curlr%Re(:,:,1:mes_D%pN)**2d0+vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlz%Re(:,:,1:mes_D%pN)**2d0)
      GTD_theta%Re(:,:,1:mes_D%pN)=dacos(-vel_curlz%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))
      GTD_sin%Re(:,:,1:mes_D%pN)=vel_curlt_comb%Re(:,:,1:mes_D%pN)/dsqrt(vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      GTD_cos%Re(:,:,1:mes_D%pN)=vel_curlr%Re(:,:,1:mes_D%pN)/dsqrt(vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      GTD_sinsin%Re(:,:,1:mes_D%pN)=vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0/&
      (vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      GTD_coscos%Re(:,:,1:mes_D%pN)=vel_curlr%Re(:,:,1:mes_D%pN)**2d0/&
      (vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)
      GTD_sincos%Re(:,:,1:mes_D%pN)=vel_curlr%Re(:,:,1:mes_D%pN)*vel_curlt_comb%Re(:,:,1:mes_D%pN)/&
      (vel_curlt_comb%Re(:,:,1:mes_D%pN)**2d0+vel_curlr%Re(:,:,1:mes_D%pN)**2d0)

      ! ReCompute eigenvalue for non-diagonalisable
      call GTD_prevar()
      call GTD_eigsys(1) ! With interpolation
      !print*,'ok 7', mpi_rnk
    end subroutine GTD_compute
    subroutine GTD_eavg()
      transform_matrix(:,:)=DCMPLX(0.0)
      transform_matrix(2,3)=DCMPLX(-1.0)
      do K=1,i_pN
        do J=0,i_Th-1
          do I=0,i_pZ-1
              X(2)=GTD_theta%Re(I,J,K)
              X(3)=GTD_S%Re(I,J,K)/d_dr
              do LL=1,3
                ! call fint_BBRe(LL)%evaluate(X(3),X(2),Y)
                call fint_e(LL)%evaluate(X(3),X(2),0,0,Y,iflag)
                e123(LL)=Y
              end do

              transform_matrix(1,1)=DCMPLX(GTD_sin%Re(I,J,K))
              transform_matrix(1,2)=DCMPLX(-GTD_cos%Re(I,J,K))
              transform_matrix(3,1)=DCMPLX(GTD_cos%Re(I,J,K))
              transform_matrix(3,2)=DCMPLX(GTD_sin%Re(I,J,K))

              ertz=MATMUL(e123,DBLE(transform_matrix))

              GTD_er%Re(I,J,K)=ertz(1)
              GTD_et%Re(I,J,K)=ertz(2)
              GTD_ez%Re(I,J,K)=ertz(3)

          end do
        end do
      end do
    end subroutine GTD_eavg
    subroutine GTD_prevar()
      traceG2%CMP=SQRT(DCMPLX( &
        (vel_Grr%Re**2d0+vel_Gtt%Re**2d0+vel_Gzz%Re**2d0 &
        +2d0*(vel_Grt%Re*vel_Gtr%Re &
            + vel_Grz%Re*vel_Gzr%Re &
            + vel_Gzt%Re*vel_Gtz%Re)) &
        /6d0))
      traceG2%CMP(:,:,1:mes_D%pN)=traceG2%CMP(:,:,1:mes_D%pN)/DCMPLX(GTD_S%Re(:,:,1:mes_D%pN))

      det%Re=(vel_Grr%Re*vel_Gtt%Re*vel_Gzz%Re &
        +vel_Grt%Re*vel_Gtz%Re*vel_Gzr%Re &
        +vel_Grz%Re*vel_Gtr%Re*vel_Gzt%Re &
        -vel_Grr%Re*vel_Gtz%Re*vel_Gzt%Re &
        -vel_Grt%Re*vel_Gtr%Re*vel_Gzz%Re &
        -vel_Grz%Re*vel_Gtt%Re*vel_Gzr%Re)
      det%Re(:,:,1:mes_D%pN)=det%Re(:,:,1:mes_D%pN)/(GTD_S%Re(:,:,1:mes_D%pN)**3d0)

      halfdetB%CMP(:,:,1:mes_D%pN)= DCMPLX(det%Re(:,:,1:mes_D%pN))/(traceG2%CMP(:,:,1:mes_D%pN)**3d0)/DCMPLX(2d0)
    end subroutine GTD_prevar
    subroutine GTD_eigsys(F)
      INTEGER, INTENT(IN) :: F
      phi%CMP(:,:,1:mes_D%pN)=ACOS(halfdetB%CMP(:,:,1:mes_D%pN))/DCMPLX(3d0)
      beta1%CMP=DCMPLX(2d0)*COS(phi%CMP)*traceG2%CMP
      beta2%CMP=DCMPLX(2d0)*COS(phi%CMP+DCMPLX(2d0*d_PI/3d0))*traceG2%CMP
      beta3%CMP=-beta2%CMP-beta1%CMP

      GTD_eig_imag%Re=max(ABS(imag(beta1%CMP)),ABS(imag(beta2%CMP)))
      ! Temporary solution to eigen swapping Problem
      ! Because the library assumed eigenvalue in the order of [+i,-i,0]

      do K=1,mes_D%pN
        do J=0,i_Th-1
          do I=0,i_pZ-1
            if (abs(imag(beta3%CMP(I,J,K)))>abs(imag(beta2%CMP(I,J,K))) &
            .or. abs(imag(beta3%CMP(I,J,K)))>abs(imag(beta1%CMP(I,J,K))) ) then
              if (abs(imag(beta1%CMP(I,J,K)))>=abs(imag(beta2%CMP(I,J,K)))) then
                beta_temp=beta2%CMP(I,J,K)
                beta2%CMP(I,J,K)=beta3%CMP(I,J,K)
                beta3%CMP(I,J,K)=beta_temp
              else
                beta_temp=beta1%CMP(I,J,K)
                beta1%CMP(I,J,K)=beta3%CMP(I,J,K)
                beta3%CMP(I,J,K)=beta_temp
              endif
              if ((imag(beta1%CMP(I,J,K)))<(imag(beta2%CMP(I,J,K)))) then
                beta_temp=beta1%CMP(I,J,K)
                beta1%CMP(I,J,K)=beta2%CMP(I,J,K)
                beta2%CMP(I,J,K)=beta_temp
              endif
            endif
          end do
        end do
      end do

      ! Rotate Matrices (Checked)
      vel_G11%Re=GTD_coscos%Re*vel_Gtt%Re-GTD_sincos%Re*(vel_Grt%Re+vel_Gtr%Re) &
      +GTD_sinsin%Re*vel_Grr%Re
      vel_G12%Re=GTD_cos%Re*vel_Gtz%Re-GTD_sin%Re*vel_Grz%Re
      vel_G13%Re=-GTD_coscos%Re*vel_Gtr%Re+GTD_sincos%Re*(vel_Grr%Re-vel_Gtt%Re) &
      +GTD_sinsin%Re*vel_Grt%Re
      vel_G21%Re=GTD_cos%Re*vel_Gzt%Re-GTD_sin%Re*vel_Gzr%Re
      vel_G22%Re=vel_Gzz%Re
      vel_G23%Re=-GTD_cos%Re*vel_Gzr%Re-GTD_sin%Re*vel_Gzt%Re
      vel_G31%Re=-GTD_coscos%Re*vel_Grt%Re+GTD_sincos%Re*(vel_Grr%Re-vel_Gtt%Re) &
      +GTD_sinsin%Re*vel_Gtr%Re
      vel_G32%Re=-GTD_sin%Re*vel_Gtz%Re-GTD_cos%Re*vel_Grz%Re
      vel_G33%Re=GTD_coscos%Re*vel_Grr%Re+GTD_sincos%Re*(vel_Grt%Re+vel_Gtr%Re) &
      +GTD_sinsin%Re*vel_Gtt%Re
      ! if (mpi_rnk==0 .and. F==0) then
      !   print*, vel_G11%Re(1,1,1), vel_G12%Re(1,1,1), vel_G13%Re(1,1,1)
      !   print*, vel_G21%Re(1,1,1), vel_G22%Re(1,1,1), vel_G23%Re(1,1,1)
      !   print*, vel_G31%Re(1,1,1), vel_G32%Re(1,1,1), vel_G33%Re(1,1,1)
      !   print*, GTD_S%Re(1,1,1)
      ! end if
      ! Compute eigenvector
      GG1%Re(:,:,1:mes_D%pN)=vel_G11%Re(:,:,1:mes_D%pN)*vel_G11%Re(:,:,1:mes_D%pN) &
      +vel_G12%Re(:,:,1:mes_D%pN)*vel_G21%Re(:,:,1:mes_D%pN) &
      +vel_G13%Re(:,:,1:mes_D%pN)*vel_G31%Re(:,:,1:mes_D%pN)
      GG2%Re(:,:,1:mes_D%pN)=vel_G21%Re(:,:,1:mes_D%pN)*vel_G11%Re(:,:,1:mes_D%pN) &
      +vel_G22%Re(:,:,1:mes_D%pN)*vel_G21%Re(:,:,1:mes_D%pN) &
      +vel_G23%Re(:,:,1:mes_D%pN)*vel_G31%Re(:,:,1:mes_D%pN)
      GG3%Re(:,:,1:mes_D%pN)=vel_G31%Re(:,:,1:mes_D%pN)*vel_G11%Re(:,:,1:mes_D%pN) &
      +vel_G32%Re(:,:,1:mes_D%pN)*vel_G21%Re(:,:,1:mes_D%pN) &
      +vel_G33%Re(:,:,1:mes_D%pN)*vel_G31%Re(:,:,1:mes_D%pN)

      W11%CMP(:,:,1:mes_D%pN)=DCMPLX(GG1%Re(:,:,1:mes_D%pN) &
      /(GTD_S%Re(:,:,1:mes_D%pN)**2d0))+beta1%CMP(:,:,1:mes_D%pN) &
      *DCMPLX(vel_G11%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))+beta2%CMP(:,:,1:mes_D%pN)*beta3%CMP(:,:,1:mes_D%pN)
      W21%CMP(:,:,1:mes_D%pN)=DCMPLX(GG2%Re(:,:,1:mes_D%pN) &
      /(GTD_S%Re(:,:,1:mes_D%pN)**2d0))+beta1%CMP(:,:,1:mes_D%pN) &
      *DCMPLX(vel_G21%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))
      W31%CMP(:,:,1:mes_D%pN)=DCMPLX(GG3%Re(:,:,1:mes_D%pN) &
      /(GTD_S%Re(:,:,1:mes_D%pN)**2d0))+beta1%CMP(:,:,1:mes_D%pN) &
      *DCMPLX(vel_G31%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))

      GG1%Re(:,:,1:mes_D%pN)=vel_G11%Re(:,:,1:mes_D%pN) &
      *vel_G12%Re(:,:,1:mes_D%pN)+vel_G12%Re(:,:,1:mes_D%pN) &
      *vel_G22%Re(:,:,1:mes_D%pN)+vel_G13%Re(:,:,1:mes_D%pN)*vel_G32%Re(:,:,1:mes_D%pN)
      GG2%Re(:,:,1:mes_D%pN)=vel_G21%Re(:,:,1:mes_D%pN) &
      *vel_G12%Re(:,:,1:mes_D%pN)+vel_G22%Re(:,:,1:mes_D%pN) &
      *vel_G22%Re(:,:,1:mes_D%pN)+vel_G23%Re(:,:,1:mes_D%pN)*vel_G32%Re(:,:,1:mes_D%pN)
      GG3%Re(:,:,1:mes_D%pN)=vel_G31%Re(:,:,1:mes_D%pN) &
      *vel_G12%Re(:,:,1:mes_D%pN)+vel_G32%Re(:,:,1:mes_D%pN) &
      *vel_G22%Re(:,:,1:mes_D%pN)+vel_G33%Re(:,:,1:mes_D%pN)*vel_G32%Re(:,:,1:mes_D%pN)

      W12%CMP(:,:,1:mes_D%pN)=DCMPLX(GG1%Re(:,:,1:mes_D%pN) &
      /(GTD_S%Re(:,:,1:mes_D%pN)**2d0))+beta2%CMP(:,:,1:mes_D%pN) &
      *DCMPLX(vel_G12%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))
      W22%CMP(:,:,1:mes_D%pN)=DCMPLX(GG2%Re(:,:,1:mes_D%pN) &
      /(GTD_S%Re(:,:,1:mes_D%pN)**2d0))+beta2%CMP(:,:,1:mes_D%pN) &
      *DCMPLX(vel_G22%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))+beta1%CMP(:,:,1:mes_D%pN)*beta3%CMP(:,:,1:mes_D%pN)
      W32%CMP(:,:,1:mes_D%pN)=DCMPLX(GG3%Re(:,:,1:mes_D%pN) &
      /(GTD_S%Re(:,:,1:mes_D%pN)**2d0))+beta2%CMP(:,:,1:mes_D%pN) &
      *DCMPLX(vel_G32%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))

      GG1%Re(:,:,1:mes_D%pN)=vel_G11%Re(:,:,1:mes_D%pN) &
      *vel_G13%Re(:,:,1:mes_D%pN)+vel_G12%Re(:,:,1:mes_D%pN)* &
      vel_G23%Re(:,:,1:mes_D%pN)+vel_G13%Re(:,:,1:mes_D%pN)*vel_G33%Re(:,:,1:mes_D%pN)
      GG2%Re(:,:,1:mes_D%pN)=vel_G21%Re(:,:,1:mes_D%pN) &
      *vel_G13%Re(:,:,1:mes_D%pN)+vel_G22%Re(:,:,1:mes_D%pN)* &
      vel_G23%Re(:,:,1:mes_D%pN)+vel_G23%Re(:,:,1:mes_D%pN)*vel_G33%Re(:,:,1:mes_D%pN)
      GG3%Re(:,:,1:mes_D%pN)=vel_G31%Re(:,:,1:mes_D%pN) &
      *vel_G13%Re(:,:,1:mes_D%pN)+vel_G32%Re(:,:,1:mes_D%pN)* &
      vel_G23%Re(:,:,1:mes_D%pN)+vel_G33%Re(:,:,1:mes_D%pN)*vel_G33%Re(:,:,1:mes_D%pN)

      W13%CMP(:,:,1:mes_D%pN)=DCMPLX(GG1%Re(:,:,1:mes_D%pN) &
      /(GTD_S%Re(:,:,1:mes_D%pN)**2d0))+beta3%CMP(:,:,1:mes_D%pN) &
      *DCMPLX(vel_G13%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))
      W23%CMP(:,:,1:mes_D%pN)=DCMPLX(GG2%Re(:,:,1:mes_D%pN) &
      /(GTD_S%Re(:,:,1:mes_D%pN)**2d0))+beta3%CMP(:,:,1:mes_D%pN) &
      *DCMPLX(vel_G23%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))
      W33%CMP(:,:,1:mes_D%pN)=DCMPLX(GG3%Re(:,:,1:mes_D%pN) &
      /(GTD_S%Re(:,:,1:mes_D%pN)**2d0))+beta3%CMP(:,:,1:mes_D%pN) &
      *DCMPLX(vel_G33%Re(:,:,1:mes_D%pN)/GTD_S%Re(:,:,1:mes_D%pN))+beta1%CMP(:,:,1:mes_D%pN)*beta2%CMP(:,:,1:mes_D%pN)

      do K=1,mes_D%pN
        do J=0,i_Th-1
          do I=0,i_pZ-1
            if (F==0 .or. GTD_diag_bool%Re(I,J,K)/=0d0) then
              X(1)=GTD_eig_imag%Re(I,J,K)
              X(2)=GTD_theta%Re(I,J,K)
              X(3)=GTD_S%Re(I,J,K)/d_dr

              ! If (I==1 .and. J==1 .and. K==1 .and. F==0 .and. mpi_rnk==0) print*, X
              W(1,1)=W11%CMP(I,J,K)
              W(2,1)=W21%CMP(I,J,K)
              W(3,1)=W31%CMP(I,J,K)
              W(1,2)=W12%CMP(I,J,K)
              W(2,2)=W22%CMP(I,J,K)
              W(3,2)=W32%CMP(I,J,K)
              W(1,3)=W13%CMP(I,J,K)
              W(2,3)=W23%CMP(I,J,K)
              W(3,3)=W33%CMP(I,J,K)

              ! NORMALISE EIGVEC
              W1_norm%CMP(I,J,K)=sqrt(W(1,1)*W(1,1)+W(2,1)*W(2,1)+W(3,1)*W(3,1))
              W2_norm%CMP(I,J,K)=sqrt(W(1,2)*W(1,2)+W(2,2)*W(2,2)+W(3,2)*W(3,2))
              W3_norm%CMP(I,J,K)=sqrt(W(1,3)*W(1,3)+W(2,3)*W(2,3)+W(3,3)*W(3,3))

              W(:,1)=W(:,1)/W1_norm%CMP(I,J,K)
              W(:,2)=W(:,2)/W2_norm%CMP(I,J,K)
              W(:,3)=W(:,3)/W3_norm%CMP(I,J,K)

              Winv=InvMat(W)
              ! If (I==1 .and. J==1 .and. K==1 .and. F==0 .and. mpi_rnk==0) then
              !   Print*,'W'
              !   print*, W(1,:)
              !   print*, W(2,:)
              !   print*, W(3,:)
              !   Print*,'Winv'
              !   print*, Winv(1,:)
              !   print*, Winv(2,:)
              !   print*, Winv(3,:)
              !   Print*,'------------------------------------------'
              ! end if
              ! print*,X
              do LL=1,27
                ! call fint_BeRe(LL)%evaluate(X(3),X(2),X(1),Y)
                call fint_BeRe(LL)%evaluate(X(3),X(2),X(1),0,0,0,Y,iflag)
                if (iflag/=0) print*, 'Interpolation error',iflag
                BeReint(LL)=Y
              end do
              do LL=1,27
                ! call fint_BeIm(LL)%evaluate(X(3),X(2),X(1),Y)
                call fint_BeIm(LL)%evaluate(X(3),X(2),X(1),0,0,0,Y,iflag)
                if (iflag/=0) print*, 'Interpolation error',iflag
                BeImint(LL)=Y
              end do
              do LL=1,45
                ! call fint_BBRe(LL)%evaluate(X(3),X(2),X(1),Y)
                call fint_BBRe(LL)%evaluate(X(3),X(2),X(1),0,0,0,Y,iflag)
                if (iflag/=0) print*, 'Interpolation error',iflag
                BBReint(LL)=Y
              end do
              do LL=1,45
                ! call fint_BBIm(LL)%evaluate(X(3),X(2),X(1),Y)
                call fint_BBIm(LL)%evaluate(X(3),X(2),X(1),0,0,0,Y,iflag)
                if (iflag/=0) print*, 'Interpolation error',iflag
                BBImint(LL)=Y
              end do

              Be(:,:)=(0d0,0d0)
              BB(:,:)=(0d0,0d0)
              ! Be (Checked)ó
              do RR=1,3
                do JJ=1,3
                  do qq=1,3
                    do PP=1,3
                      Be(pp,qq)=Be(pp,qq)+Winv(rr,pp)*W(jj,rr)* &
                      DCMPLX(BeReint(BeInd(RR,QQ,JJ)),BeImint(BeInd(RR,QQ,JJ)))
                    end do
                  end do
                end do
              end do
              ! BB (checked)
              do RR=1,3
                do SS=1,3
                  do JJ=1,3
                    do KK=1,3
                      BB(RR,SS)=BB(RR,SS)+W(KK,SS)*W(jj,rr)* &
                      DCMPLX(BBReint(BBInd(RR,SS,JJ,KK)),BBImint(BBInd(RR,SS,JJ,KK)))
                    end do
                  end do
                end do
              end do
              BB(:,1)=BB(:,1)*(beta1%CMP(I,J,K))
              BB(:,2)=BB(:,2)*(beta2%CMP(I,J,K))
              BB(:,3)=BB(:,3)*(beta3%CMP(I,J,K))
              BB_trans=MATMUL(TRANSPOSE(Winv),MATMUL(BB,Winv))
              ! If (I==1 .and. J==1 .and. K==1 .and. F==0 .and. mpi_rnk==0) then
              !   print*, Be(1,:)
              !   print*, Be(2,:)
              !   print*, Be(3,:)
              !   print*, BB_trans(1,:)
              !   print*, BB_trans(2,:)
              !   print*, BB_trans(3,:)
              ! end if
              !Checked
              Diff=Be+BB_trans*DCMPLX(X(3))

              transform_matrix(:,:)=DCMPLX(0.0)
              transform_matrix(1,1)=DCMPLX(GTD_sin%Re(I,J,K))
              transform_matrix(1,2)=DCMPLX(-GTD_cos%Re(I,J,K))
              transform_matrix(2,3)=DCMPLX(-1.0)
              transform_matrix(3,1)=DCMPLX(GTD_cos%Re(I,J,K))
              transform_matrix(3,2)=DCMPLX(GTD_sin%Re(I,J,K))

              Diff=(transpose(Diff)+Diff)/DCMPLX(2d0)
              Diff2=MATMUL(transpose(transform_matrix),Diff)
              Diff=MATMUL(Diff2,transform_matrix)

              !Checked
              if (F==0) then
                GTD_Drr%Re(I,J,K)=DBLE(Diff(1,1))
                GTD_Drt%Re(I,J,K)=DBLE(Diff(1,2))
                GTD_Drz%Re(I,J,K)=DBLE(Diff(1,3))
                GTD_Dtt%Re(I,J,K)=DBLE(Diff(2,2))
                GTD_Dtz%Re(I,J,K)=DBLE(Diff(2,3))
                GTD_Dzz%Re(I,J,K)=DBLE(Diff(3,3))
              else
                GTD_Drr%Re(I,J,K)=(GTD_Drr%Re(I,J,K)+DBLE(Diff(1,1)))/2d0
                GTD_Drt%Re(I,J,K)=(GTD_Drt%Re(I,J,K)+DBLE(Diff(1,2)))/2d0
                GTD_Drz%Re(I,J,K)=(GTD_Drz%Re(I,J,K)+DBLE(Diff(1,3)))/2d0
                GTD_Dtt%Re(I,J,K)=(GTD_Dtt%Re(I,J,K)+DBLE(Diff(2,2)))/2d0
                GTD_Dtz%Re(I,J,K)=(GTD_Dtz%Re(I,J,K)+DBLE(Diff(2,3)))/2d0
                GTD_Dzz%Re(I,J,K)=(GTD_Dzz%Re(I,J,K)+DBLE(Diff(3,3)))/2d0
              end if
            end if
          end do
        end do
      end do
    end subroutine GTD_eigsys
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
    function BeInd(RRR,QQQ,JJJ)
      INTEGER, INTENT(IN) :: RRR,QQQ,JJJ
      INTEGER :: BeInd
        BeInd=JJJ+3*(RRR-1)+9*(QQQ-1)
        RETURN
    end function BeInd
    function BBInd(RRR,SSS,JJJ,KKK)
      INTEGER, INTENT(IN) :: RRR,SSS,JJJ,KKK
      INTEGER :: BBInd
      INTEGER :: JR, KS, DIA
        JR=JJJ+3*RRR-3
        KS=KKK+3*SSS-3
        DIA=ABS(JR-KS)
        BBInd=dia*(19-dia)/2+min(jr,ks);
        RETURN
    end function BBInd
    subroutine GTD_precompute()
      GTD_TSOURCE%RE(:,:,:)=1d-1
      GTD_FSOURCE%RE(:,:,:)=0d0
      ! Open the CDF file
      call GTD_open()

      ! Read-in whole library
      call GTD_read('BeRe',BeRe_mul,27)
      call GTD_read('BeIm',BeIm_mul,27)
      call GTD_read('BBRe',BBRe_mul,45)
      call GTD_read('BBIm',BBIm_mul,45)
      call GTD_read_e('e1_array',e1_mul)
      call GTD_read_e('e2_array',e2_mul)
      call GTD_read_e('e3_array',e3_mul)
      ! Close the CDF file
      call GTD_close()

#ifdef _MPI
      call MPI_BCAST(BeRe_mul, 27*S_len*theta_len*eig_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(BeIm_mul, 27*S_len*theta_len*eig_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(BBRe_mul, 45*S_len*theta_len*eig_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(BBIm_mul, 45*S_len*theta_len*eig_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(e1_mul, S_len*theta_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(e2_mul, S_len*theta_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(e3_mul, S_len*theta_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
#endif
      if (mpi_rnk==0) print*, 'GTD variables Broadcase finished.'
      ! Initialise interpolation object
      ALLOCATE(fint_S(S_len))
      ALLOCATE(fint_theta(theta_len))
      ALLOCATE(fint_eig(eig_len))
      fint_S=(/ ((I-1)*S_step,I=1, S_len) /)
      fint_theta=(/( (I-1)*theta_step,I=1, theta_len )/)
      fint_eig=(/( (I-1)*eig_step,I=1, eig_len) /)

      Do I=1,27
        ! call fint_BeRe(I)%initialize(fint_S,fint_theta,fint_eig,BeRe_mul(:,:,:,I),iflag)
        call fint_BeRe(I)%initialize(fint_S,fint_theta,fint_eig,BeRe_mul(:,:,:,I),4,4,4,iflag)
        if (iflag/=0) print*, 'Initialising error',iflag
        ! call fint_BeIm(I)%initialize(fint_S,fint_theta,fint_eig,BeIm_mul(:,:,:,I),iflag)
        call fint_BeIm(I)%initialize(fint_S,fint_theta,fint_eig,BeIm_mul(:,:,:,I),4,4,4,iflag)
        if (iflag/=0) print*, 'Initialising error',iflag
      end do
      Do I=1,45
        ! call fint_BBRe(I)%initialize(fint_S,fint_theta,fint_eig,BBRe_mul(:,:,:,I),iflag)
        call fint_BBRe(I)%initialize(fint_S,fint_theta,fint_eig,BBRe_mul(:,:,:,I),4,4,4,iflag)
        if (iflag/=0) print*, 'Initialising error',iflag
        ! call fint_BBIm(I)%initialize(fint_S,fint_theta,fint_eig,BBIm_mul(:,:,:,I),iflag)
        call fint_BBIm(I)%initialize(fint_S,fint_theta,fint_eig,BBIm_mul(:,:,:,I),4,4,4,iflag)
        if (iflag/=0) print*, 'Initialising error',iflag
      end do

      call fint_e(1)%initialize(fint_S,fint_theta,e1_mul(:,:),4,4,iflag)
      call fint_e(2)%initialize(fint_S,fint_theta,e2_mul(:,:),4,4,iflag)
      call fint_e(3)%initialize(fint_S,fint_theta,e3_mul(:,:),4,4,iflag)

    end subroutine GTD_precompute
    subroutine GTD_open()
      INTEGER :: e , id, dimid
      if (mpi_rnk==0) then
        e=nf90_open('GTD_lib.cdf',nf90_nowrite, file_id) !f is the returned netCDF id
        if(e/=nf90_noerr) then
           stop 'io: file not found!'
        end if

        e=nf90_get_att(file_id,nf90_global,'S_step'    , S_step)     !function nf90_get_att(ncid, varid, name, values)
        e=nf90_get_att(file_id,nf90_global,'theta_step', theta_step) !function nf90_get_att(ncid, varid, name, values)
        e=nf90_get_att(file_id,nf90_global,'eig_step'  , eig_step)   !function nf90_get_att(ncid, varid, name, values)

        e=nf90_inq_dimid(file_id,'S', dimid)
        e=nf90_inquire_dimension(file_id, dimid, len=S_len)
        e=nf90_inq_dimid(file_id,'theta', dimid)
        e=nf90_inquire_dimension(file_id, dimid, len=theta_len)
        e=nf90_inq_dimid(file_id,'eig', dimid)
        e=nf90_inquire_dimension(file_id, dimid, len=eig_len)
      end if
#ifdef _MPI
      call MPI_BCAST(S_len, 1, MPI_INTEGER, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(theta_len, 1, MPI_INTEGER, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(eig_len, 1, MPI_INTEGER, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(S_step, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(theta_step, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(eig_step, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
#endif
    end subroutine GTD_open

    subroutine GTD_read(nm,res,res_len)
      INTEGER :: e, id, dimid
      character(*),     intent(in)  :: nm
      integer, intent(in)  :: res_len
      DOUBLE PRECISION, intent(out), allocatable :: res(:,:,:,:)

      allocate(res(S_len,theta_len,eig_len,res_len))
      if (mpi_rnk/=0) return
      e=nf90_inq_varid(file_id,nm, id)
      e=nf90_get_var(file_id,id,res)
      print*, nm, ' read.'
    end subroutine GTD_read
    subroutine GTD_read_e(nm,res)
      INTEGER :: e, id, dimid
      character(*),     intent(in)  :: nm
      DOUBLE PRECISION, intent(out), allocatable :: res(:,:)

      allocate(res(S_len,theta_len))
      if (mpi_rnk/=0) return
      e=nf90_inq_varid(file_id,nm, id)
      e=nf90_get_var(file_id,id,res)
      print*, nm, ' read.'
    end subroutine GTD_read_e

    subroutine GTD_close()
      INTEGER :: e
      if (mpi_rnk/=0) return
          e=nf90_close(file_id)
          print*, 'GTD_lib read and closed.'
    end subroutine GTD_close
  end module GTD
