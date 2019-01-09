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

   type (phys), private :: p,p1,p2,p3
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
      type (phys) :: p11,p22,p33
      double precision :: d
       d= d_Ri /2d0/ d_nint
         			! advection  u x curlu
      p11%Re = vel_t%Re*vel_curlz%Re - vel_z%Re*vel_curlt%Re
      p22%Re = vel_z%Re*vel_curlr%Re - vel_r%Re*vel_curlz%Re
      p33%Re = vel_r%Re*vel_curlt%Re - vel_t%Re*vel_curlr%Re +d*dexp(temp_p%Re)
      call tra_phys2coll(p11,vel_Nr, p22,vel_Nt, p33,vel_Nz)

      call non_addHPF()

   end subroutine non_velocity

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine non_addHPF()
      double precision :: a(i_N), b(i_N)
      _loop_km_vars
                  			! force from background HPF and from T
      b = -vel_Up
      _loop_km_begin
         a = d_alpha*k * vel_U
         vel_Nr%Re(:,nh) = vel_Nr%Re(:,nh) + a*vel_ur%Im(:,nh)
         vel_Nr%Im(:,nh) = vel_Nr%Im(:,nh) - a*vel_ur%Re(:,nh)
         vel_Nt%Re(:,nh) = vel_Nt%Re(:,nh) + a*vel_ut%Im(:,nh)
         vel_Nt%Im(:,nh) = vel_Nt%Im(:,nh) - a*vel_ut%Re(:,nh)
         vel_Nz%Re(:,nh) = vel_Nz%Re(:,nh) + a*vel_uz%Im(:,nh)  &
                                           + b*vel_ur%Re(:,nh)
         vel_Nz%Im(:,nh) = vel_Nz%Im(:,nh) - a*vel_uz%Re(:,nh)  &
                                           + b*vel_ur%Im(:,nh)

      _loop_km_end
       					! additional pressure if fixed flx
      if(b_const_flux .and. mpi_rnk==0) then
         vel_Pr0 = dot_product(vel_uz%Re(i_N-i_KL:,0),mes_D%dr1(:,1))
         vel_Pr0 = vel_Pr0/(-2d0)
         vel_Nz%Re(:,0) = vel_Nz%Re(:,0) + 4d0*vel_Pr0/d_Re
      end if


      if(mpi_rnk/=0) return
				! zero mode real
      vel_Nr%Im(:,0) = 0d0
      vel_Nt%Im(:,0) = 0d0
      vel_Nz%Im(:,0) = 0d0

   end subroutine non_addHPF


!-------------------------------------------------------------------------
!  nonlinear terms for the tempertaure
!-------------------------------------------------------------------------
   subroutine non_temperature()
      double precision :: a(i_N), c, beta, Pe
      _loop_km_vars
      Pe=d_Pr*d_Re
      beta=d_beta*d_Vs
				! advection temperature, -u.grad(tau)
      ! p%Re = -(vel_r%Re-beta*vel_curlt%Re)*temp_gradr%Re  &
      !       - (vel_t%Re+beta*vel_curlr%Re)*temp_gradt%Re  &
      !       - vel_z%Re*temp_gradz%Re  &
      !       - beta*vel_lapz%Re &
      !       - beta*vel_Up_phy%Re*temp_gradr%Re &
      !       + (temp_gradr%Re*temp_gradr%Re &
      !       + temp_gradt%Re*temp_gradt%Re &
      !       +temp_gradz%Re*temp_gradz%Re)/Pe
      p1%Re=-(vel_r%Re+beta*(vel_Up_phy%Re-vel_curlt%Re)-temp_gradr%Re/Pe)
      p2%Re=-(vel_t%Re+beta*vel_curlr%Re-temp_gradt%Re/Pe)
      p3%Re=-(vel_z%Re-temp_gradz%Re/Pe)

      p%Re =   p1%Re*temp_gradr%Re  &
             + p2%Re*temp_gradt%Re  &
             + p3%Re*temp_gradz%Re  &
             - beta*vel_lapz%Re
      call tra_phys2spec(p, s)
      call var_spec2coll(s, temp_N)

      _loop_km_begin
         a = d_alpha*k * vel_U

         temp_N%Re(:,nh) = temp_N%Re(:,nh) + a*temp_tau%Im(:,nh)

         temp_N%Im(:,nh) = temp_N%Im(:,nh) - a*temp_tau%Re(:,nh)

      _loop_km_end
      				! zero mode real
      if(mpi_rnk/=0) return
      temp_N%Re(:,0) = temp_N%Re(:,0)-beta*vel_Upp
      temp_N%Im(:,0) = 0d0

   end subroutine non_temperature
   !------------------------------------------------------------------------
   !  get cfl max dt due to flow field
   !------------------------------------------------------------------------
      subroutine non_maxtstep()
         double precision :: d,mx, dt(6),dt_(6), r(i_N)
         integer :: n, n__

         r = mes_D%r(:,1)
         dt = 1d99

         do n = 1, mes_D%pN
            n__ = n+mes_D%pNi-1

            if(n__==1) then
               d = r(2) - r(1)
            else if(n__==i_N) then
               d = r(i_N) - r(i_N-1)
            else
               d = min( r(n__)-r(n__-1), r(n__+1)-r(n__) )
            end if
            mx = maxval( dabs(vel_r%Re(:,:,n)) )
            if(mx/=0d0) dt(1) = min( dt(1), d/mx )
            mx = maxval( dabs(p1%Re(:,:,n)) )
            if(mx/=0d0) dt(4) = min( dt(4), d/mx )

            d = 2d0*d_PI/dble(i_Th*i_Mp) 		!---  *r_(n)? ---
            mx = maxval( dabs(vel_t%Re(:,:,n)) )
            if(mx/=0d0) dt(2) = min( dt(2), d/mx )
            mx = maxval( dabs(p2%Re(:,:,n)) )
            if(mx/=0d0) dt(5) = min( dt(5), d/mx )

            d = 2d0*d_PI/(d_alpha*i_Z)
            mx = maxval( dabs(vel_z%Re(:,:,n) + vel_U(n__)) )
            if(mx/=0d0) dt(3) = min( dt(3), d/mx )
            mx = maxval( dabs(p3%Re(:,:,n) + vel_U(n__)) )
            if(mx/=0d0) dt(6) = min( dt(6), d/mx )
         end do

#ifdef _MPI
         call mpi_allreduce(dt, dt_, 6, mpi_double_precision,  &
            mpi_min, mpi_comm_world, mpi_er)
         dt = dt_
#endif
         tim_cfl_dt = minval(dt)
         if(tim_cfl_dt==dt(1)) tim_cfl_dir=1
         if(tim_cfl_dt==dt(2)) tim_cfl_dir=2
         if(tim_cfl_dt==dt(3)) tim_cfl_dir=3
         if(tim_cfl_dt==dt(4)) tim_cfl_dir=4
         if(tim_cfl_dt==dt(5)) tim_cfl_dir=5
         if(tim_cfl_dt==dt(6)) tim_cfl_dir=6

      end subroutine non_maxtstep
!*************************************************************************
 end module nonlinear
!*************************************************************************
