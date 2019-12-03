!***************************************************************************
! parameters
!***************************************************************************
!
! N		n in [1,N] radial r_n
! L		l in [0,L) theta
! M		m in [0,M) phi
! Mp		index m=0,1,2,.. corresponds to m=0,Mp,2Mp,...
!
!***************************************************************************
#include "../parallel.h"
 module parameters
!***************************************************************************
   implicit none
   save

   integer,          parameter :: i_N           = 80
   integer,          parameter :: i_K           = 16
   integer,          parameter :: i_M           = 1
   integer,          parameter :: i_Mp          = 1

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   double precision            :: d_Re          = 0.401070456591576d0
   double precision            :: d_Ri          = 4.322797117392170d0
   double precision, parameter :: d_Vs          = 0.314159265358979d0        ! Changes with flow rate, as always normalised by bulk flow
   double precision, parameter :: d_dr          = 0.668211770763543d0
   double precision            :: d_Pe          = 1d0/d_Vs/d_Vs*d_dr
   double precision            :: d_Pe_dm       = 3d0/d_Vs/d_Vs*d_dr
   double precision            :: d_alpha       = 0.1d0
   logical,          parameter :: b_const_flux  = .true.
   logical,          parameter :: b_mirrorsym   = .false.
   logical,          parameter :: b_shiftrefl   = .false.
   logical,          parameter :: b_shiftrott   = .false.
   logical,          parameter :: b_loadstate   = .true.
   logical,          parameter :: b_const_mean  = .false.
   double precision, parameter :: d_minE3d      = 1d-5

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   integer,          parameter :: i_save_rate1  = 5000!100
   integer,          parameter :: i_save_rate2  = 500
   integer,          parameter :: i_maxtstep    = 49000000
   double precision, parameter :: d_maxt        = -1d0
   double precision, parameter :: d_cpuhours    = 71.7d0
   double precision, parameter :: d_time        = -1d0
   double precision, parameter :: d_timestep    = -1d0 !1d-3 !0.05d0
   double precision, parameter :: d_maxdt       = 1d-4 !1d99
   double precision, parameter :: d_dterr       = 1d-8 !1d99
   double precision, parameter :: d_dterr_scalar= 2d-9 !1d99
   double precision, parameter :: d_dterr_bc    = 1d99 !1d99
   double precision, parameter :: d_courant     = 0.45d0
   double precision, parameter :: d_implicit    = 0.51d0
   integer,          parameter :: d_dterr_transient_iter = 6

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   integer,          parameter :: i_KL  = 4
   integer,          parameter :: i_K1  = i_K-1
   integer,          parameter :: i_M1  = i_M-1
   integer,          parameter :: i_Z   = 3*i_K !3/2 dealiasing rule, Can't be changed
   integer,          parameter :: i_Th  = 3*i_M !3/2 dealiasing rule, Can't be changed
   integer,          parameter :: i_H1  = (2*i_K1+1)*i_M-i_K1-1
   integer,          parameter :: i_pN  = (_Nr+i_N-1)/_Nr
   integer,          parameter :: i_pH1 = (_Nr+_Hs1)/_Nr-1
   integer,          parameter :: i_pZ  = i_Z/_Ns
   double precision, parameter :: d_PI  = 3.1415926535897931d0

 contains

!---------------------------------------------------------------------------
!  check params are ok
!---------------------------------------------------------------------------
   subroutine par_precompute()
      if(modulo(i_Z,_Ns)/=0) stop '_Ns must divide i_Z'
      if(modulo(i_M,_Ns)/=0) stop '_Ns must divide i_M'
      if(i_pH1==0) stop 'Too many cores='
   end subroutine par_precompute


!***************************************************************************
 end module parameters
!***************************************************************************
