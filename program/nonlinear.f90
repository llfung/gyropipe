!**************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
 module nonlinear
!*************************************************************************
   use velocity
   use temperature
   use transform
   implicit none
   save

   type (phys), private :: p
   type (spec), private :: s

 contains


!------------------------------------------------------------------------
!  Initialise nonlinear stuff
!------------------------------------------------------------------------
   subroutine non_precompute()
   end subroutine non_precompute


!-------------------------------------------------------------------------
!  nonlinear terms for velocity
!------------------------------------------------------------------------
!   subroutine non_velocity()
!      double precision :: d
!      _loop_km_vars
!         			! advection  u x curlu
!      p%Re = vel_y%Re*vel_curlz%Re - vel_z%Re*vel_curly%Re
!      call tra_phys2spec(p, s)
!      call var_spec2coll(s, vel_Nx)
!      p%Re = vel_z%Re*vel_curlx%Re - vel_x%Re*vel_curlz%Re
!      call tra_phys2spec(p, s)
!      call var_spec2coll(s, vel_Ny)
!      p%Re = vel_x%Re*vel_curly%Re - vel_y%Re*vel_curlx%Re
!      call tra_phys2spec(p, s)
!      call var_spec2coll(s, vel_Nz)
!				! force from C
!      d = d_Ra * d_Pr
!      _loop_km_begin
!         vel_Ny%Re(:,nh) = vel_Ny%Re(:,nh) + d*temp_T%Re(:,nh)
!         vel_Ny%Im(:,nh) = vel_Ny%Im(:,nh) + d*temp_T%Im(:,nh)
!      _loop_km_end

!      if(mpi_rnk/=0) return
!      				! from background C0
!      vel_Ny%Re(:,0) = vel_Ny%Re(:,0) + d*temp_T0
!				! zero mode real
!      vel_Nx%Im(:,0) = 0d0
!      vel_Ny%Im(:,0) = 0d0
!      vel_Nz%Im(:,0) = 0d0

!   end subroutine non_velocity


!-------------------------------------------------------------------------
!  nonlinear terms for velocity
!------------------------------------------------------------------------
   subroutine non_velocity()
      type (phys) :: p1,p2,p3
         			! advection  u x curlu
      p1%Re = vel_t%Re*vel_curlz%Re - vel_z%Re*vel_curlt%Re
      p2%Re = vel_z%Re*vel_curlr%Re - vel_r%Re*vel_curlz%Re
      p3%Re = vel_r%Re*vel_curlt%Re - vel_t%Re*vel_curlr%Re
      call tra_phys2coll(p1,vel_Nr, p2,vel_Nt, p3,vel_Nz)

      call non_addHPF()

   end subroutine non_velocity

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine non_addHPF()
      double precision :: a(i_N), b(i_N)
      double precision :: d
      _loop_km_vars
                  			! force from background HPF and from T
      b = -vel_Up
      d = d_Ri
!      if(mpi_rnk==0) print*, 'd=',d
      _loop_km_begin
         a = d_alpha*k * vel_U
         vel_Nr%Re(:,nh) = vel_Nr%Re(:,nh) + a*vel_ur%Im(:,nh)
         vel_Nr%Im(:,nh) = vel_Nr%Im(:,nh) - a*vel_ur%Re(:,nh)
         vel_Nt%Re(:,nh) = vel_Nt%Re(:,nh) + a*vel_ut%Im(:,nh)
         vel_Nt%Im(:,nh) = vel_Nt%Im(:,nh) - a*vel_ut%Re(:,nh)
         vel_Nz%Re(:,nh) = vel_Nz%Re(:,nh) + a*vel_uz%Im(:,nh)  &
                                           + b*vel_ur%Re(:,nh)  &
                                           + d*temp_tau%Re(:,nh)
         vel_Nz%Im(:,nh) = vel_Nz%Im(:,nh) - a*vel_uz%Re(:,nh)  &
                                           + b*vel_ur%Im(:,nh)  &
                                           + d*temp_tau%Im(:,nh)

      _loop_km_end
       					! additional pressure if fixed flx
      if(b_const_flux .and. mpi_rnk==0) then
         vel_Pr0 = dot_product(vel_uz%Re(i_N-i_KL:,0),mes_D%dr1(:,1))
         vel_Pr0 = vel_Pr0/(-2d0)
         vel_Nz%Re(:,0) = vel_Nz%Re(:,0) + 4d0*vel_Pr0/d_Re
      end if


      if(mpi_rnk/=0) return
      				! from background T0
      !vel_Nz%Re(:,0) = vel_Nz%Re(:,0) + d*temp_T0(:)
				! zero mode real
      vel_Nr%Im(:,0) = 0d0
      vel_Nt%Im(:,0) = 0d0
      vel_Nz%Im(:,0) = 0d0

   end subroutine non_addHPF


!-------------------------------------------------------------------------
!  nonlinear terms for the tempertaure
!-------------------------------------------------------------------------
   subroutine non_temperature()
      double precision :: a(i_N-1), c, beta
      INTEGER :: R_index
      _loop_km_vars

      beta=d_beta*d_Vs
				! advection temperature, -u.grad(tau)
      p%Re = -(vel_r%Re-beta*vel_curlt%Re)*temp_gradr%Re  &
            - (vel_t%Re+beta*vel_curlr%Re)*temp_gradt%Re  &
            - vel_z%Re*temp_gradz%Re  &
            - beta*vel_lapz%Re*temp_p%Re &
            - beta*vel_Up_phy%Re*temp_gradr%Re

      ! Boundary Condition at r=R
      if(mpi_rnk==_Np-1) then
        R_index=mes_D%pN+mes_D%pNi-1
        p%Re(:,:,R_index)=temp_gradr%Re(:,:,R_index)/d_Re/d_Pr/beta/ &
          (vel_Up_phy%Re(:,:,R_index)-vel_curlt%Re(:,:,R_index))*d_BC
      end if
      call tra_phys2spec(p, s)
      call var_spec2coll(s, temp_N)

      ! Due to B.C., only compute 1:(i_N-1)
      c = -beta*vel_Upp
      _loop_km_begin
         a = d_alpha*k * vel_U(1:i_N-1)

         temp_N%Re(1:(i_N-1),nh) = temp_N%Re(1:(i_N-1),nh) + a*temp_tau%Im(1:(i_N-1),nh)  &
                         + c*temp_tau%Re(1:(i_N-1),nh)

         temp_N%Im(1:(i_N-1),nh) = temp_N%Im(1:(i_N-1),nh) - a*temp_tau%Re(1:(i_N-1),nh)  &
                         + c*temp_tau%Im(1:(i_N-1),nh)

         ! a = d_alpha*k * vel_U(1:i_N)
         !
         ! temp_N%Re(1:(i_N),nh) = temp_N%Re(1:(i_N),nh) + a*temp_tau%Im(1:(i_N),nh)  &
         !                 + c*temp_tau%Re(1:(i_N),nh)
         !
         ! temp_N%Im(1:(i_N),nh) = temp_N%Im(1:(i_N),nh) - a*temp_tau%Re(1:(i_N),nh)  &
         !                 + c*temp_tau%Im(1:(i_N),nh)
      _loop_km_end
      				! zero mode real
      if(mpi_rnk/=0) return
      temp_N%Im(:,0) = 0d0
!      temp_N%Re(:,0) = temp_N%Re(:,0) + 4d0/(d_Re*d_Pr)

   end subroutine non_temperature


!*************************************************************************
 end module nonlinear
!*************************************************************************
