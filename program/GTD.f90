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

    type (phys) :: vel_G11, vel_G12, vel_G13
    type (phys) :: vel_G21, vel_G22, vel_G23
    type (phys) :: vel_G31, vel_G32, vel_G33

    type (phys) :: GTD_D11, GTD_D12, GTD_D13
    type (phys) :: GTD_D22, GTD_D23, GTD_D33

    type (phys) :: GG1, GG2, GG3

    type (phys_cmp) :: W11, W12, W13
    type (phys_cmp) :: W21, W22, W23
    type (phys_cmp) :: W31, W32, W33
    type (phys_cmp) :: W1_norm, W2_norm, W3_norm



    DOUBLE PRECISION :: BeReint(27),BeImint(27),BBReint(45),BBImint(45)
    DOUBLE COMPLEX :: Be(3,3),BB(3,3),W(3,3),Winv(3,3)
    DOUBLE COMPLEX :: BB_trans(3,3),Diff(3,3),Diff2(3,3),transform_matrix(3,3)
    type (phys) :: det
    type (phys_cmp) :: p1, p2, phi
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
    INTEGER :: iflag

    INTEGER           :: file_id
    INTEGER           :: S_len, theta_len, eig_len
    DOUBLE PRECISION              :: S_step, theta_step, eig_step
    DOUBLE PRECISION, ALLOCATABLE :: BeRe_mul(:,:,:,:), BeIm_mul(:,:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: BBRe_mul(:,:,:,:), BBIm_mul(:,:,:,:)
  contains
    ! Main Algorithm
    subroutine GTD_compute()
      ! Compute eigenvalue
      det%Re=(vel_Grr%Re*vel_Gtt%Re*vel_Gzz%Re &
        +vel_Grt%Re*vel_Gtz%Re*vel_Gzr%Re &
        +vel_Grz%Re*vel_Gtr%Re*vel_Gzt%Re &
        -vel_Grr%Re*vel_Gtz%Re*vel_Gzt%Re &
        -vel_Grt%Re*vel_Gtr%Re*vel_Gzz%Re &
        -vel_Grz%Re*vel_Gtt%Re*vel_Gzr%Re)

      p1%CMP=SQRT(DCMPLX( &
        (vel_Grr%Re**2d0+vel_Gtt%Re**2d0+vel_Gzz%Re**2d0 &
        +2d0*(vel_Grt%Re*vel_Gtr%Re &
            + vel_Grz%Re*vel_Gzr%Re &
            + vel_Gzt%Re*vel_Gtz%Re)) &
        /6d0))
      p2%CMP= DCMPLX(det%Re)/(p1%CMP**3d0)/DCMPLX(2d0)

      phi%CMP=ACOS(p2%CMP)/DCMPLX(3d0)
      beta1%CMP=DCMPLX(2d0)*COS(phi%CMP)*p1%CMP
      beta2%CMP=DCMPLX(2d0)*COS(phi%CMP+DCMPLX(2d0*d_PI/3d0))*p1%CMP
      beta3%CMP=-beta2%CMP-beta1%CMP

      GTD_eig_imag%Re=max(ABS(imag(beta1%CMP)),ABS(imag(beta2%CMP)))
      ! Temporary solution to eigen swapping Problem
      ! Because the library assumed eigenvalue in the order of [+i,-i,0]
      do K=1,i_pN
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
      vel_G23%%ReRe=-%ReGTD_cos%Re*vel_Gzr%Re-GTD_sin%Re*vel_Gzt%Re
      vel_G31%Re=-GTD_coscos%Re*vel_Grt%Re+GTD_sincos%Re*(vel_Grr%Re-vel_Gtt%Re) &
      +GTD_sinsin%Re*vel_Gtr%Re
      vel_G32%Re=-GTD_sin%Re*vel_Gtz%Re-GTD_cos%Re*vel_Grz%Re
      vel_G33%Re=GTD_coscos%Re*vel_Grr%Re+GTD_sincos%Re*(vel_Grt%Re+vel_Gtr%Re) &
      +GTD_sinsin%Re*vel_Gtt%Re

      ! call save_vel_local('G123.cdf')
      ! Compute eigenvector
      GG1%Re=vel_G11%Re*vel_G11%Re+vel_G12%Re*vel_G21%Re+vel_G13%Re*vel_G31%Re
      GG2%Re=vel_G21%Re*vel_G11%Re+vel_G22%Re*vel_G21%Re+vel_G23%Re*vel_G31%Re
      GG3%Re=vel_G31%Re*vel_G11%Re+vel_G32%Re*vel_G21%Re+vel_G33%Re*vel_G31%Re

      W11%CMP=DCMPLX(GG1%Re)+beta1%CMP*DCMPLX(vel_G11%Re)+beta2%CMP*beta3%CMP
      W21%CMP=DCMPLX(GG2%Re)+beta1%CMP*DCMPLX(vel_G21%Re)
      W31%CMP=DCMPLX(GG3%Re)+beta1%CMP*DCMPLX(vel_G31%Re)

      W12%CMP=DCMPLX(GG1%Re)+beta2%CMP*DCMPLX(vel_G11%Re)+beta1%CMP*beta3%CMP
      W22%CMP=DCMPLX(GG2%Re)+beta2%CMP*DCMPLX(vel_G21%Re)
      W32%CMP=DCMPLX(GG3%Re)+beta2%CMP*DCMPLX(vel_G31%Re)

      W13%CMP=DCMPLX(GG1%Re)+beta3%CMP*DCMPLX(vel_G11%Re)+beta1%CMP*beta2%CMP
      W23%CMP=DCMPLX(GG2%Re)+beta3%CMP*DCMPLX(vel_G21%Re)
      W33%CMP=DCMPLX(GG3%Re)+beta3%CMP*DCMPLX(vel_G31%Re)


      do K=1,i_pN
        do J=0,i_Th-1
          do I=0,i_pZ-1
              X(1)=GTD_eig_imag%Re(I,J,K)
              X(2)=GTD_theta%Re(I,J,K)
              X(3)=GTD_S%Re(I,J,K)

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
              ! Winv=W
              ! call ZGESVD('N','N',3,3,Winv,3,SVD_S,SVD_U,1,SVD_VT,1,WORK,LWMAX,RWORK,INFO)
              ! rcond(I,J,K)=SVD_S(1)/SVD_S(3)
              Winv=W
              Call ZGETRF(3,3,Winv,3,IPVT,INFO)
              if (INFO/=0) print*, 'CGETRF error: ', Info, beta1(I,J,K)%CMP, beta2(I,J,K)%CMP,beta3(I,J,K)%CMP
              Call ZGETRI(3,Winv,3,IPVT,WORK,LWMAX,INFO)
              if (INFO/=0) print*, 'CGETRI error: ', Info, beta1(I,J,K)%CMP, beta2(I,J,K)%CMP,beta3(I,J,K)%CMP

              do L=1,27
                ! call fint_BeRe(L)%evaluate(X(3),X(2),X(1),Y)
                call fint_BeRe(L)%evaluate(X(3),X(2),X(1),0,0,0,Y,iflag)
                BeReint(L)=Y
              end do
              do L=1,27
                ! call fint_BeIm(L)%evaluate(X(3),X(2),X(1),Y)
                call fint_BeIm(L)%evaluate(X(3),X(2),X(1),0,0,0,Y,iflag)
                BeImint(L)=Y
              end do
              do L=1,45
                ! call fint_BBRe(L)%evaluate(X(3),X(2),X(1),Y)
                call fint_BBRe(L)%evaluate(X(3),X(2),X(1),0,0,0,Y,iflag)
                BBReint(L)=Y
              end do
              do L=1,45
                ! call fint_BBIm(L)%evaluate(X(3),X(2),X(1),Y)
                call fint_BBIm(L)%evaluate(X(3),X(2),X(1),0,0,0,Y,iflag)
                BBImint(L)=Y
              end do

              Be(:,:)=(0d0,0d0)
              BB(:,:)=(0d0,0d0)
              ! Be (Checked)
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

              !Checked
              Diff=Be+BB_trans*DCMPLX(X(3))

              transform_matrix(:,:)=DCMPLX(0.0)
              transform_matrix(1,1)=DCMPLX(GTD_sin(I,J,K))
              transform_matrix(1,2)=DCMPLX(-GTD_cos(I,J,K))
              transform_matrix(2,3)=DCMPLX(-1.0)
              transform_matrix(3,1)=DCMPLX(GTD_cos(I,J,K))
              transform_matrix(3,2)=DCMPLX(GTD_sin(I,J,K))

              Diff=(transpose(Diff)+Diff)/DCMPLX(2d0)
              Diff2=MATMUL(transpose(transform_matrix),Diff)
              Diff=MATMUL(Diff2,transform_matrix)

              !Checked
              GTD_D11%Re(I,J,K)=DBLE(Diff(1,1))
              GTD_D12%Re(I,J,K)=DBLE(Diff(1,2))
              GTD_D13%Re(I,J,K)=DBLE(Diff(1,3))
              GTD_D22%Re(I,J,K)=DBLE(Diff(2,2))
              GTD_D23%Re(I,J,K)=DBLE(Diff(2,3))
              GTD_D33%Re(I,J,K)=DBLE(Diff(3,3))

          end do
        end do
      end do

    end subroutine GTD_compute
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
      ! Open the CDF file
      call GTD_open()
      ! Read-in whole library
      call GTD_read('BeRe',BeRe_mul,27)
      call GTD_read('BeIm',BeIm_mul,27)
      call GTD_read('BBRe',BBRe_mul,45)
      call GTD_read('BBIm',BBIm_mul,45)

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
      ! Close the CDF file
      call GTD_close()
    end subroutine GTD_precompute
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
      e=nf90_inq_dimid(file_id,nm, dimid)
      allocate(res(S_len,theta_len,eig_len,res_len))

      e=nf90_inq_varid(file_id,nm, id)
      e=nf90_get_var(file_id,id,res)
    end subroutine GTD_read

    subroutine GTD_close()
      INTEGER :: e
          e=nf90_close(file_id)
    end subroutine GTD_close
  end module GTD
