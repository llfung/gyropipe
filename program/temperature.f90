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
   implicit none
   save

   type (phys) :: temp_gradr
   type (phys) :: temp_gradt
   type (phys) :: temp_gradz
   type (phys) :: temp_p   !temperature perturbation (physical)
   type (coll) :: temp_tau !temperature perturbation
   type (coll) :: temp_N !nonlinear terms temperature eq temp_n = -u.grad(tau)

  ! double precision :: temp_T0(i_N) !temperature basic state T0 = 1 -r^2
  ! double precision :: temp_T0p(i_N) !temperature gradient basic state dT/dr = -2r

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
      if (mpi_rnk/=0) return
      temp_tau%Re(:,0)=(1d0 - mes_D%r(:,2))*d_Re*d_Pr*d_beta*d_Vs
    !  temp_T0  =  1d0 - mes_D%r(:,2)	! 1 - r^2
    !  temp_T0p = - 2d0 * mes_D%r(:,1)	! dT/dr
   end subroutine temp_precompute


!------------------------------------------------------------------------
!  precomputation of matrices for timestepping
!------------------------------------------------------------------------
   subroutine temp_matrices()
      double precision :: d1, d2
      integer :: j,nl,nr,n
      _loop_km_vars

			! lhs matrices
      d1 =  1d0/tim_dt
      d2 = -d_implicit/(d_Re*d_Pr)
      call tim_lumesh_init_mod(0,1,d1,d2, LD)

      			! timestepping matrices for rhs
      d1 =  1d0/tim_dt
      d2 =  (1d0-d_implicit)/(d_Re*d_Pr)
      call tim_mesh_init(0,d1,d2, Lt)

      ! For B.C. at r=R
      _loop_km_begin
        do j = 1-i_KL, i_N
           nl = max(1,j-i_KL)
           nr = min(j+i_KL,i_N)
           do n = nl, nr
              if (n==i_N) then
                Lt(nh)%M(i_KL+1+n-j,j)=0d0
              end if
           end do
        end do
      _loop_km_end
   end subroutine temp_matrices


!------------------------------------------------------------------------
!  Evaluate in physical space  grad(tau)
!------------------------------------------------------------------------
   subroutine temp_transform()

      call var_coll_grad(temp_tau, c1,c2,c3)
      call tra_coll2phys(c1,temp_gradr, c2,temp_gradt, c3,temp_gradz)

      call tra_coll2phys(temp_tau,temp_p)

   end subroutine temp_transform


!------------------------------------------------------------------------
!  Advance equations one timestep.
!------------------------------------------------------------------------
   subroutine temp_step()
				     	! get rhs = A u_ + N
       ! print*, 'T Non-linear Re: ', maxval(temp_N%Re), mpi_rnk
       ! print*, 'T Non-linear Im: ', maxval(temp_N%Im), mpi_rnk
      call tim_meshmult(0,Lt,T_,temp_N, temp_tau) !tim_meshmult(S,A,b,c, d)
						!  multiply  d = A b + c
						!  S=0, b even for m even; S=1, b odd for m even

       ! print*, 'Before B.C. : ', maxval(temp_tau%Re), mpi_rnk
      call temp_tempbc(temp_tau)
        				! invert
       ! print*, 'Before: ', maxval(temp_tau%Re), mpi_rnk
      call tim_lumesh_invert(0,LD, temp_tau)
       ! print*, 'After: ', maxval(temp_tau%Re), mpi_rnk
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
  subroutine temp_tempbc(a)
     type (coll), intent(inout) :: a
     integer :: n
     double precision :: fac, gradrRe,gradrIm
     fac=d_Re*d_Pr*d_beta*d_Vs
     do n = 0, var_H%pH1
        gradrRe=dot_product(mes_D%dr1(1:1+i_KL,1),vel_uz%Re(i_N-i_KL:i_N,n))
        gradrIm=dot_product(mes_D%dr1(1:1+i_KL,1),vel_uz%Im(i_N-i_KL:i_N,n))
        a%Re(i_N, n ) = fac*(vel_Up_col%Re(i_N,n)+gradrRe)
        a%Im(i_N, n ) = fac*(vel_Up_col%Im(i_N,n)+gradrIm)
     end do
  end subroutine temp_tempbc

!*************************************************************************
 end module temperature
!*************************************************************************
