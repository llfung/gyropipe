!*************************************************************************
! GTD Table Lookup and interpolation
!*************************************************************************
#include "../parallel.h"
  module GTD
    use netcdf
    use parameters
    implicit none
    save


    type (phys) :: GTD_S
    type (phys) :: GTD_theta
    type (phys) :: GTD_sin
    type (phys) :: GTD_cos
    type (phys) :: GTD_eig_imag


    INTEGER           :: file_id
    INTEGER           :: S_len, theta_len, eig_len
    REAL              :: S_step, theta_step, eig_step
    REAL, ALLOCATABLE :: BeRe_mul(:,:,:,:), BeIm_mul(:,:,:,:), BBRe_mul(:,:,:,:), BBIm_mul(:,:,:,:)
    REAL, ALLOCATABLE :: BeRe(:,:), BeIm(:,:), BBRe(:,:), BBIm(:,:)

    INTEGER, PARAMETER               :: NDIM = 3
    INTEGER, PARAMETER, DIMENSION(3) :: NDEG = (/ 1, 1, 1 /)
    INTEGER, PARAMETER               :: NEOPT= 2*(2*NDIM+3)+1
    INTEGER, PARAMETER, DIMENSION(3) :: LUP  = (/ 3, 3, 3 /)
    INTEGER                          :: NTAB(2*NDIM+1)
    REAL                             ::   XT(2*NDIM)
    INTEGER, PARAMETER, DIMENSION(3) :: IOPT = (/0,NEOPT,0/)
    REAL                             :: EOPT(NEOPT)
  contains
    subroutine GTD_precompute()
      ! Open the CDF file
      call GTD_open()
      ! Initialise interpolation variables for MATH77 lib
      XT(2)=eig_step
      XT(4)=theta_step
      XT(6)=S_step
      XT(1)=0.0
      XT(3)=0.0
      XT(5)=0.0

      NTAB(:)=0
      NTAB(1)=eig_len
      NTAB(2)=theta_len
      NTAB(3)=S_len
      ! Read-in whole library
      call GTD_read('BeRe',BeRe_mul,27)
      call GTD_read('BeIm',BeIm_mul,27)
      call GTD_read('BBRe',BBRe_mul,45)
      call GTD_read('BBIm',BBIm_mul,45)
      ! Reshape libraries
      ALLOCATE(BeRe(S_len*theta_len*eig_len,27))
      ALLOCATE(BeIm(S_len*theta_len*eig_len,27))
      ALLOCATE(BBRe(S_len*theta_len*eig_len,45))
      ALLOCATE(BBIm(S_len*theta_len*eig_len,45))
      BeRe=reshape(BeRe_mul(:,:,:,:),shape(BeRe))
      BeIm=reshape(BeIm_mul(:,:,:,:),shape(BeIm))
      BBRe=reshape(BBRe_mul(:,:,:,:),shape(BBRe))
      BBIm=reshape(BBIm_mul(:,:,:,:),shape(BBIm))
      ! Close the CDF file
      call GTD_close()

    end subroutine GTD_precompute

    subroutine GTD_trans()
      type (phys) :: S_horizontal, det
      DOUBLE COMPLEX :: p(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: phi(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: beta1(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: beta2(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE PRECISION :: GG1(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE PRECISION :: GG2(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE PRECISION :: GG3(0:i_pZ-1, 0:i_Th-1, i_pN)

      DOUBLE COMPLEX :: W11(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: W12(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: W13(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: W21(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: W22(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: W23(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: W31(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: W32(0:i_pZ-1, 0:i_Th-1, i_pN)
      DOUBLE COMPLEX :: W33(0:i_pZ-1, 0:i_Th-1, i_pN)

      GTD_S%Re=dsqrt(vel_curlr%Re**2d0+vel_curlt%Re**2d0+vel_curlz%Re**2d0)
      GTD_theta%Re=dacos(-vel_curlz%Re/GTD_S%Re)
      S_horizontal%Re=dsqrt(vel_curlr%Re**2d0+vel_curlt%Re**2d0)
      GTD_sin%Re=vel_curlt%Re/S_horizontal%Re
      GTD_cos%Re=vel_curlr%Re/S_horizontal%Re
      det%Re=( vel_Grr%Re*vel_Gtt%Re*vel_Gzz%Re &
              +vel_Grt%Re*vel_Gtz%Re*vel_Gzr%Re &
              +vel_Grz%Re*vel_Gtr%Re*vel_Gzt%Re &
              -vel_Grr%Re*vel_Gtz%Re*vel_Gzt%Re &
              -vel_Grt%Re*vel_Gtr%Re*vel_Gzz%Re &
              -vel_Grz%Re*vel_Gtt%Re*vel_Gzr%Re)/(GTD_S%Re**3d0)

      p=CDSQRT(DCMPLX( &
        (vel_Grr%Re**2d0+vel_Gtt%Re**2d0+vel_Gzz%Re**2d0 &
        +2d0*(vel_Grt%Re*vel_Gtr%Re &
            + vel_Grz%Re*vel_Gzr%Re &
            + vel_Gzt%Re*vel_Gtz%Re)) &
        /(GTD_S%Re**2d0)/6d0))

      phi=ACOS(DCMPLX(det%Re)/(p**3d0)/DCMPLX(2d0))/DCMPLX(3d0)
      beta1=DCMPLX(2d0)*COS(phi)*p
      beta2=DCMPLX(2d0)*COS(phi+DCMPLX(2d0*d_PI/3d0))*p
      beta3=-beta2-beta1
      GTD_eig_imag%Re=max(ABS(imag(beta1)),ABS(imag(beta2)))

      GG1=vel_Grr%Re*vel_Grr%Re+vel_Grt%Re*vel_Gtr%Re+vel_Grz%Re*vel_Gzr%Re
      GG2=vel_Gtr%Re*vel_Grr%Re+vel_Gtt%Re*vel_Gtr%Re+vel_Gtz%Re*vel_Gzr%Re
      GG3=vel_Gzr%Re*vel_Grr%Re+vel_Gzt%Re*vel_Gtr%Re+vel_Gzz%Re*vel_Gzr%Re

      W11=DCMPLX(GG1)+beta1*DCMPLX(vel_Grr%Re)+beta2*beta3
      W21=DCMPLX(GG2)+beta1*DCMPLX(vel_Gtr%Re)
      W31=DCMPLX(GG3)+beta1*DCMPLX(vel_Gzr%Re)
      W12=DCMPLX(GG1)+beta2*DCMPLX(vel_Grr%Re)+beta1*beta3
      W22=DCMPLX(GG2)+beta2*DCMPLX(vel_Gtr%Re)
      W32=DCMPLX(GG3)+beta2*DCMPLX(vel_Gzr%Re)
      W13=DCMPLX(GG1)+beta3*DCMPLX(vel_Grr%Re)+beta1*beta2
      W23=DCMPLX(GG2)+beta3*DCMPLX(vel_Gtr%Re)
      W33=DCMPLX(GG3)+beta3*DCMPLX(vel_Gzr%Re)


    end subroutine GTD_trans

    subroutine GTD_swim(G)
    end subroutine GTD_swim
    subroutine GTD_diff(G)
    end subroutine GTD_diff
    subroutine GTD_interpolate()

    end subroutine GTD_interpolate

    subroutine GTD_open()
      INTEGER :: e , id, dimid
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
    end subroutine GTD_open

    subroutine GTD_read(nm,res,res_len)
      INTEGER :: e, id, dimid
      character(*),     intent(in)  :: nm
      integer, intent(in)  :: res_len
      REAL, intent(out), allocatable :: res(:,:,:,:)
      e=nf90_inq_dimid(f,nm, dimid)
      allocate(res(S_len,theta_len,eig_len,res_len))

      e=nf90_inq_varid(f,nm, id)
      e=nf90_get_var(f,id,res)
    end subroutine GTD_read

    subroutine GTD_close()
          e=nf90_close(file_id)
    end subroutine GTD_close
  end module GTD
