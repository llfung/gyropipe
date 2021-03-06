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

   integer,	     parameter :: i_N           = 60
   integer,          parameter :: i_K           = 18
   integer,          parameter :: i_M           = 32
   integer,          parameter :: i_Mp          = 1

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   double precision            :: d_Re          = 0.4
   double precision            :: d_Pr          = 2.519526329d1
   double precision            :: d_Ri          = 20.58d0
   double precision            :: d_Vs          = 0.63d0
   double precision            :: d_beta        = 1.0251d-1
   double precision            :: d_alpha       = 1.25d0
   logical,          parameter :: b_const_flux  = .true.
   logical,          parameter :: b_mirrorsym   = .false.
   logical,          parameter :: b_shiftrefl   = .false.
   logical,          parameter :: b_shiftrott   = .false.
   double precision, parameter :: d_minE3d      = 1d-5

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   integer,          parameter :: i_save_rate1  = 100!100
   integer,          parameter :: i_save_rate2  = 100
   integer,          parameter :: i_maxtstep    = 100000
   double precision, parameter :: d_maxt        = -1d0
   double precision, parameter :: d_cpuhours    = 1d99 !90d0
   double precision, parameter :: d_time        = 0d0 !-1d0
   double precision, parameter :: d_timestep    = 0.01d0 !-1d0
   double precision, parameter :: d_maxdt       = 0.02d0 !1d99
   double precision, parameter :: d_dterr       = 1d-5 !1d99
   double precision, parameter :: d_courant     = 0.5d0
   double precision, parameter :: d_implicit    = 0.51d0
   double precision, parameter :: d_BC          = 1d5

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   integer,          parameter :: i_KL  = 4
   integer,          parameter :: i_K1  = i_K-1
   integer,          parameter :: i_M1  = i_M-1
   integer,          parameter :: i_Z   = 3*i_K
   integer,          parameter :: i_Th  = 3*i_M
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
   end subroutine par_precompute


!***************************************************************************
 end module parameters
!***************************************************************************
