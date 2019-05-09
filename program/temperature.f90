!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
 module temperature
!*************************************************************************
   use transform
   use timestep
   use velocity
   use GTD
   implicit none
   save

   type (phys) :: temp_gradr
   type (phys) :: temp_gradt
   type (phys) :: temp_gradz
   type (phys) :: temp_p   !temperature perturbation (physical)
   type (coll) :: temp_tau !temperature perturbation
   type (coll) :: temp_N !nonlinear terms temperature eq temp_n = -u.grad(tau)


   ! For Boundary Condition
   ! type (phys) :: temp_er_Drr ! e_r/Drr, for b.c.
   ! type (coll) :: temp_er_Drr_col
   type (phys_bc) :: temp_bc
   ! type (phys_bc) :: temp_bc_Drr
   ! type (phys_bc) :: temp_bc_er
   type (coll_bc) :: temp_bc_col

   double precision :: temp_T0(i_N) !temperature basic state T0 = 1 -r^2
   double precision :: temp_T0p(i_N) !temperature gradient basic state dT/dr = -2r
   double precision ::  d_nint

   type (lumesh), private :: LD(0:i_pH1)!lhs matrix
   type (mesh),   private :: Lt(0:i_pH1) !rhs matrix for timestepping
   type (coll),   private :: N_, T_ !predictor

   type (coll), private :: c1,c2,c3
   type (phys), private :: p
   type (spec), private :: s

 contains

!------------------------------------------------------------------------
!  initialise field
!------------------------------------------------------------------------
   subroutine temp_precompute()
      call var_coll_init(temp_tau)
      temp_bc%Re=0d0
      temp_bc_col%Re=0d0
      temp_bc_col%Im=0d0
      temp_T0  = 0d0
      temp_T0p = 0d0

      if (mpi_rnk==0) then
        temp_tau%Re(:,0)=(- mes_D%r(:,2))*d_dr/d_Vs
        d_nint = dot_product(dexp(temp_tau%Re(:,0)),mes_D%intrdr)
      end if
#ifdef _MPI
      call MPI_BCAST(d_nint, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_er)
#endif
    !  temp_T0  =  1d0 - mes_D%r(:,2)	! 1 - r^2
    !  temp_T0p = - 2d0 * mes_D%r(:,1)	! dT/dr
   end subroutine temp_precompute


!------------------------------------------------------------------------
!  precomputation of matrices for timestepping
!------------------------------------------------------------------------
   subroutine temp_matrices()
      double precision :: d1, d2
      integer :: nl,nr,n
      ! _loop_km_vars

			! lhs matrices
      d1 =  1d0/tim_dt
      d2 = -d_implicit/d_Pe_dm
      call tim_lumesh_init_mod(0,1,d1,d2, LD)

      			! timestepping matrices for rhs
      d1 =  1d0/tim_dt
      d2 =  (1d0-d_implicit)/d_Pe_dm
      call tim_mesh_init(0,d1,d2, Lt)

      call temp_adjustMean(0)
      ! For B.C. at r=R
      ! _loop_km_begin
      !   do j = 1-i_KL, i_N
      !      nl = max(1,j-i_KL)
      !      nr = min(j+i_KL,i_N)
      !      do n = nl, nr
      !         if (n==i_N) then
      !           Lt(nh)%M(i_KL+1+n-j,j)=0d0
      !         end if
      !      end do
      !   end do
      ! _loop_km_end
   end subroutine temp_matrices

!------------------------------------------------------------------------
!  adjust mean Temp; fix to d_nint
!------------------------------------------------------------------------
   subroutine temp_adjustMean(F)
      integer, intent(in) :: F
      double precision, save :: Ti(i_N), d1,d2,d3
      integer :: info

      if(mpi_rnk/=0) return

      if(F==0) then
         Ti(:)   = 1d0
         Ti(i_N) = 0d0
         call dgbtrs('N', i_N, i_KL, i_KL, 1, LD(0)%M, 3*i_KL+1,  &
                     LD(0)%ipiv, Ti, i_N, info)
         if(info/=0) stop 'temp_adjustMean: err 1'
         d1 = 2d0*dot_product(Ti, mes_D%intrdr)
         if(d1==0d0) stop 'temp_adjustMean: err 2'

      else if(F==1) then
         d2 = 2d0*(dot_product(dexp(temp_tau%Re(:,0)), mes_D%intrdr)-d_nint)
         d3 = -d2/d1 ! Need to work out how to compensate for exp(H) lost.
         temp_tau%Re(:,0) = dlog(dexp(temp_tau%Re(:,0)) + d3*Ti)
      end if

   end subroutine temp_adjustMean


!------------------------------------------------------------------------
!  Evaluate in physical space  grad(tau)
!------------------------------------------------------------------------
   subroutine temp_transform()

      call var_coll_grad(temp_tau, c1,c2,c3)
      call tra_coll2phys(c1,temp_gradr, c2,temp_gradt, c3,temp_gradz)

      call tra_coll2phys(temp_tau,temp_p)

!#ifdef _MPI
!
!        if(mpi_rnk==0) d_nint = dot_product(dexp(temp_tau%Re(:,0)),mes_D%intrdr)
!        call mpi_bcast(d_nint,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
!
!#else
!      d_nint = dot_product(dexp(temp_tau%Re(:,0)),mes_D%intrdr)
!#endif
   end subroutine temp_transform


!------------------------------------------------------------------------
!  Advance equations one timestep.
!------------------------------------------------------------------------
   subroutine temp_step()
				    ! get rhs = A u_ + N
      call tim_meshmult(0,Lt,T_,temp_N, temp_tau) !tim_meshmult(S,A,b,c, d)
						!  multiply  d = A b + c
						!  S=0, b even for m even; S=1, b odd for m even

      call temp_tempbc(temp_tau)
        		! invert
            ! Modify LD to BC need (LD=LD_original*Drr+Drt*im+Drz*ialpha)
            ! if (mpi_rnk==0) print*, 'before', temp_tau%Re(:,0)
      call tim_lumesh_invert(0,LD, temp_tau)
            ! if (mpi_rnk==0) print*, 'after', temp_tau%Re(:,0)

      call temp_adjustMean(1)

      if(mpi_rnk==0)  &
         temp_tau%Im(:,0) = 0d0

   end subroutine temp_step
!------------------------------------------------------------------------
!  predictor with euler nonlinear terms
!------------------------------------------------------------------------
   subroutine temp_predictor()

      call var_coll_copy(temp_tau, T_)
      call var_coll_copy(temp_N, N_)
      call temp_step()

   end subroutine temp_predictor


!------------------------------------------------------------------------
!  corrector iteration with Crank-Nicolson non-lin term
!------------------------------------------------------------------------
   subroutine temp_corrector()

      call tim_nlincorr(N_, temp_N)
      call temp_step()

   end subroutine temp_corrector

!------------------------------------------------------------------------
!  B.C. no-flux at r=1
!------------------------------------------------------------------------
  subroutine temp_tempbc(ain)
     type (coll), intent(inout) :: ain
     integer :: n

     do n = 0, var_H%pH1
        ain%Re(i_N, n ) = temp_bc_col%Re(1,n)
        ain%Im(i_N, n ) = temp_bc_col%Im(1,n)
     end do

  end subroutine temp_tempbc

!*************************************************************************
 end module temperature
!*************************************************************************
