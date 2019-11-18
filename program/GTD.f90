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

    INTEGER        :: II,JJ,KK, insize
    LOGICAL        :: EXTRAPOLATE_FLAG

    INTEGER            :: file_id
    INTEGER, PARAMETER :: G12_len = 97
    INTEGER, PARAMETER :: G21_len = 13
    INTEGER, PARAMETER :: G11_len = 13
    INTEGER, PARAMETER :: G22_len = 13
    INTEGER, PARAMETER :: omega_e3_len = 205

    DOUBLE PRECISION :: G12(G12_len),G21(G21_len),G11(G11_len),G22(G22_len)
    DOUBLE PRECISION :: G11_lim(2),G12_lim(2),G21_lim(2),G22_lim(2)
    DOUBLE PRECISION :: D11(G12_len,G21_len,G11_len,G22_len)
    DOUBLE PRECISION :: D12(G12_len,G21_len,G11_len,G22_len)
    DOUBLE PRECISION :: D22(G12_len,G21_len,G11_len,G22_len)
    DOUBLE PRECISION :: D11VV((G12_len+2)*(G21_len+2)*(G11_len+2)*(G22_len+2))
    DOUBLE PRECISION :: D12VV((G12_len+2)*(G21_len+2)*(G11_len+2)*(G22_len+2))
    DOUBLE PRECISION :: D22VV((G12_len+2)*(G21_len+2)*(G11_len+2)*(G22_len+2))

    DOUBLE PRECISION :: e1(omega_e3_len),e2(omega_e3_len),omg3(omega_e3_len)
  contains
!------------------------------------------------------------------------
!  Algorithm for computing D_T from G=Grad(u) at BC for no-flux condition of H
!------------------------------------------------------------------------
    subroutine GTD_compute_bc()
      DOUBLE PRECISION :: loc_G(i_pZ*i_Th)
      DOUBLE PRECISION :: loc_D11(i_pZ*i_Th,1,1),loc_D12(i_pZ*i_Th,1,1)
      DOUBLE PRECISION ::  loc_e1(i_pZ*i_Th,1,1)
      DOUBLE PRECISION :: loc_Drr(i_pZ*i_Th,1,1),loc_Drz(i_pZ*i_Th,1,1),loc_Dzz(i_pZ*i_Th,1,1)
      DOUBLE PRECISION ::  loc_er(i_pZ*i_Th,1,1),loc_ez(i_pZ*i_Th,1,1)

      if (mpi_rnk/=(_Nr-1)) return
      loc_G(:)=0d0

      insize=i_pZ*i_Th
      ! call interp4lin(G12,G21,G11,G22,D11, &
      ! -vel_Grz%Re(:,:,mes_D%pN)/d_dr,loc_G(:,1),vel_Grr%Re(:,:,mes_D%pN)/d_dr,loc_G(:,2), loc_D11, insize)
      ! call interp4lin(G12,G21,G11,G22,D12, &
      ! -vel_Grz%Re(:,:,mes_D%pN)/d_dr,loc_G(:,1),vel_Grr%Re(:,:,mes_D%pN)/d_dr,loc_G(:,2), loc_D12, insize)

      call interp4_interp(G12,G21,G11,G22,D11VV, &
      -vel_Grz%Re(:,:,mes_D%pN)/d_dr,loc_G,vel_Grr%Re(:,:,mes_D%pN)/d_dr,loc_G, loc_D11, insize)
      call interp4_interp(G12,G21,G11,G22,D12VV, &
      -vel_Grz%Re(:,:,mes_D%pN)/d_dr,loc_G,vel_Grr%Re(:,:,mes_D%pN)/d_dr,loc_G, loc_D12, insize)
      call interp1in(omg3,e1,-vel_Grz%Re(:,:,mes_D%pN)/d_dr,loc_e1, insize)


      GTD_Drr_bc%Re= RESHAPE(loc_D11,(/i_pZ,i_Th/))
      GTD_Drt_bc%Re=0d0
      GTD_Drz_bc%Re=-RESHAPE(loc_D12,(/i_pZ,i_Th/))
      GTD_er_bc %Re= RESHAPE(loc_e1,(/i_pZ,i_Th/))


    end subroutine GTD_compute_bc
!------------------------------------------------------------------------
!  Main Algorithm for computing D_T from G=Grad(u)
!------------------------------------------------------------------------
    subroutine GTD_compute()
      type(phys) :: loc_omg3
      DOUBLE PRECISION :: loc_D11(mes_D%pN*i_pZ*i_Th,1,1),loc_D12(mes_D%pN*i_pZ*i_Th,1,1),loc_D22(mes_D%pN*i_pZ*i_Th,1,1)
      DOUBLE PRECISION ::  loc_e1(mes_D%pN*i_pZ*i_Th,1,1), loc_e2(mes_D%pN*i_pZ*i_Th,1,1)


      if ((MINVAL( vel_Grr%Re)<G11_lim(1) .or. MAXVAL( vel_Grr%Re)>G11_lim(2) .or. &
           MINVAL( vel_Gzz%Re)<G22_lim(1) .or. MAXVAL( vel_Gzz%Re)>G22_lim(2) .or. &
           MINVAL(-vel_Grz%Re)<G12_lim(1) .or. MAXVAL(-vel_Grz%Re)>G12_lim(2) .or. &
           MINVAL(-vel_Gzr%Re)<G21_lim(1) .or. MAXVAL(-vel_Gzr%Re)>G21_lim(2) ) &
       .and. .NOT.(EXTRAPOLATE_FLAG)) then
        print*,' Extrapolating GTD!'

        if (MINVAL( vel_Grr%Re)<G11_lim(1)) then
          print*,' Grr: ', MINVAL( vel_Grr%Re),' smaller than ', G11_lim(1)
        end if
        if (MINVAL(-vel_Grz%Re)<G12_lim(1)) then
          print*,'-Grz: ', MINVAL(-vel_Grz%Re),' smaller than ', G12_lim(1)
        end if
        if (MINVAL(-vel_Gzr%Re)<G21_lim(1)) then
          print*,'-Gzr: ', MINVAL(-vel_Gzr%Re),' smaller than ', G21_lim(1)
        end if
        if (MINVAL( vel_Gzz%Re)<G22_lim(1)) then
          print*,' Gzz: ', MINVAL( vel_Gzz%Re),' smaller than ', G22_lim(1)
        end if
        if (MAXVAL( vel_Grr%Re)>G11_lim(2)) then
          print*,' Grr: ', MAXVAL( vel_Grr%Re),'  bigger than ', G11_lim(2)
        end if
        if (MAXVAL(-vel_Grz%Re)>G12_lim(2)) then
          print*,'-Grz: ', MAXVAL(-vel_Grz%Re),'  bigger than ', G12_lim(2)
        end if
        if (MAXVAL(-vel_Gzr%Re)>G21_lim(2)) then
          print*,'-Gzr: ', MAXVAL(-vel_Gzr%Re),'  bigger than ', G21_lim(2)
        end if
        if (MAXVAL( vel_Gzz%Re)>G22_lim(2)) then
          print*,' Gzz: ', MAXVAL( vel_Gzz%Re),'  bigger than ', G22_lim(2)
        end if
        EXTRAPOLATE_FLAG=.TRUE.
      end if

      insize= mes_D%pN*i_pZ*i_Th
      loc_omg3%Re=-vel_Grz%Re/d_dr+vel_Gzr%Re/d_dr

      ! call interp4lin(G12,G21,G11,G22,D11, &
      ! loc_G12,loc_G21,loc_G11,loc_G22, loc_D11, insize)
      ! call interp4lin(G12,G21,G11,G22,D12, &
      ! loc_G12,loc_G21,loc_G11,loc_G22, loc_D12, insize)
      ! call interp4lin(G12,G21,G11,G22,D22, &
      ! loc_G12,loc_G21,loc_G11,loc_G22, loc_D22, insize)

      call interp4_interp(G12,G21,G11,G22,D11VV, &
      -vel_Grz%Re/d_dr,-vel_Gzr%Re/d_dr,vel_Grr%Re/d_dr,vel_Gzz%Re/d_dr, loc_D11, insize)
      call interp4_interp(G12,G21,G11,G22,D12VV, &
      -vel_Grz%Re/d_dr,-vel_Gzr%Re/d_dr,vel_Grr%Re/d_dr,vel_Gzz%Re/d_dr, loc_D12, insize)
      call interp4_interp(G12,G21,G11,G22,D22VV, &
      -vel_Grz%Re/d_dr,-vel_Gzr%Re/d_dr,vel_Grr%Re/d_dr,vel_Gzz%Re/d_dr, loc_D22, insize)

      call interp1in(omg3,e1,loc_omg3%Re,loc_e1, insize)
      call interp1in(omg3,e2,loc_omg3%Re,loc_e2, insize)

      GTD_Drr%Re= RESHAPE(loc_D11,(/i_pZ,i_Th,mes_D%pN/))
      GTD_Drz%Re=-RESHAPE(loc_D12,(/i_pZ,i_Th,mes_D%pN/))
      GTD_Dzz%Re= RESHAPE(loc_D22,(/i_pZ,i_Th,mes_D%pN/))
      GTD_er %Re= RESHAPE( loc_e1,(/i_pZ,i_Th,mes_D%pN/))
      GTD_ez %Re=-RESHAPE( loc_e2,(/i_pZ,i_Th,mes_D%pN/))

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

      !Open the CDF file
      call GTD_open()

      ! Read-in whole library
      call GTD_read_e('G11_loop',G11,G11_len)
      call GTD_read_e('G12_loop',G12,G12_len)
      call GTD_read_e('G21_loop',G21,G21_len)
      call GTD_read_e('G22_loop',G22,G22_len)
      call GTD_read_e('omega_e3_col',omg3,omega_e3_len)

      call GTD_read('D11',D11)
      call GTD_read('D12',D12)
      call GTD_read('D22',D22)
      call GTD_read_e('e1_col',e1,omega_e3_len)
      call GTD_read_e('e2_col',e2,omega_e3_len)
      ! Close the CDF file
      call GTD_close()

#ifdef _MPI
      call MPI_BCAST(G11, G11_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(G12, G12_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(G21, G21_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(G22, G22_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)

      call MPI_BCAST(D11, G12_len*G21_len*G11_len*G22_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(D12, G12_len*G21_len*G11_len*G22_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(D22, G12_len*G21_len*G11_len*G22_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)

      call MPI_BCAST(omg3, omega_e3_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(  e1, omega_e3_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
      call MPI_BCAST(  e2, omega_e3_len, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
#endif
      if (mpi_rnk==0) print*, 'GTD variables Broadcase finished.'

      ! Calculate limits
      G11_lim(1)=G11(1)*d_dr
      G11_lim(2)=G11(G11_len)*d_dr
      G12_lim(1)=G12(1)*d_dr
      G12_lim(2)=G12(G12_len)*d_dr
      G21_lim(1)=G21(1)*d_dr
      G21_lim(2)=G21(G21_len)*d_dr
      G22_lim(1)=G22(1)*d_dr
      G22_lim(2)=G22(G22_len)*d_dr

      call interp4_libgen(G12,G21,G11,G22,D11,D11VV)
      call interp4_libgen(G12,G21,G11,G22,D12,D12VV)
      call interp4_libgen(G12,G21,G11,G22,D22,D22VV)
      ! call gtd2d_libinter_cfunvec_initialize()



    end subroutine GTD_precompute
!------------------------------------------------------------------------
!  Decollocate memory stuff for the MATLAB Coder generated library
!------------------------------------------------------------------------
    subroutine GTD_closing()
        ! call gtd2d_libinter_cfunvec_terminate()
    end subroutine GTD_closing


    subroutine GTD_open()
      INTEGER :: e , id, dimid

      if (mpi_rnk/=0) return

        e=nf90_open('GTD_lib.cdf',nf90_nowrite, file_id) !f is the returned netCDF id

        if(e/=nf90_noerr) then
           stop 'io: file not found!'
        end if
    end subroutine GTD_open

    subroutine GTD_read(nm,res)
      INTEGER :: e, id, dimid
      character(*),     intent(in)  :: nm
      DOUBLE PRECISION, intent(out) :: res(G12_len,G21_len,G11_len,G22_len)

      if (mpi_rnk/=0) return
      e=nf90_inq_varid(file_id,nm, id)
      e=nf90_get_var(file_id,id,res)
      print*, nm, ' read.'
    end subroutine GTD_read
    subroutine GTD_read_e(nm,res,inp_len)
      INTEGER :: e, id, dimid
      INTEGER, INTENT(IN) :: inp_len
      character(*),     intent(in)  :: nm
      DOUBLE PRECISION, intent(out) :: res(inp_len)

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
