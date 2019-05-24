!*************************************************************************
! main.f90 (executable)
!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
 PROGRAM MAIN
!*************************************************************************
   use nonlinear
   use io
   use GTD
   implicit none
   real :: d_start, d_stop

   call initialise()

!    call vel_imposesym()
!    call vel_transform()
!
! #ifdef _MPI
!    call mpi_barrier(mpi_comm_world, mpi_er)
! #endif
!
!    call GTD_compute()
!

      						! main loop
   do while(.not.terminate())
                  ! Predictor Step
      call vel_imposesym()
      call vel_transform()
      call temp_transform()
      call non_velocity()
      call non_temperature()
      call var_null(2)          ! where you put util()
      if(d_timestep<0d0) then
         call non_maxtstep()    ! Automated timestep size calculation->tim_cfl_dt
         call tim_new_tstep()   ! Other conditions that limits the timestep size
         if(tim_new_dt)  &      ! If we have used a new timestpe compared to the last one
            call vel_matrices() !    precomputation of matrices for timestepping
            call temp_matrices()
      end if
      call var_null(1)          ! where you put util()
      call io_write2files()     ! I/O files
      call vel_predictor()      ! vel_step() called inside (the main stepping algorithm)
         call vel_transform_bc()
         call non_temperature_bc(0)
      call temp_predictor()
      tim_it = 1

      call vel_transform()
      call temp_transform()
                   ! Corrector Step iteration
      do while(tim_it/=0)
        ! call temp_transform_gradr()
        ! if (mpi_rnk==0) then
        !   print*, tim_it, temp_gradr_coll%Re(:,0)
        !   print*, tim_it, temp_tau_test%Re(:,0)
        !   print*, tim_it, temp_gradr_coll%Re(:,0)-temp_tau_test%Re(:,0)
        ! end if

         call non_velocity()
         call non_temperature()
         call var_null(2)       ! where you put util()
         call vel_corrector()   ! vel_step() called inside (the main stepping algorithm)
           call vel_transform_bc()
           call non_temperature_bc(0)
         call temp_corrector()

         call vel_transform()
         call temp_transform()
         call non_temperature_bc(1)
         call temp_check_bc()
         call tim_check_cgce()  ! check if we can exit corrector iter. If true, tim_it=0.
      end do

      tim_t    = tim_t    + tim_dt
      tim_step = tim_step + 1

   end do
      						! end main loop

   call cleanup()
   stop

 contains


!-------------------------------------------------------------------------
!  Termaination conditions
!-------------------------------------------------------------------------
   logical function terminate()
      logical :: file_exist

      if(mpi_rnk==0) then
         terminate = .false.

         if(tim_step==i_maxtstep) then
            terminate = .true.
            print*, 'maxtstep reached!'
         end if

         if(d_maxt>0d0 .and. tim_t>=d_maxt) then
            terminate = .true.
            print*, 'maxt reached!'
         end if

         call clk_time(d_stop)
         if(dble(d_stop-d_start)/36d2 >= d_cpuhours) then
            terminate = .true.
            print*, 'cpuhours reached!'
         end if

         if( modulo(tim_step,i_save_rate2)==0) then
            inquire(file='RUNNING', exist=file_exist)
            if(.not. file_exist) then
               terminate = .true.
               print*, 'RUNNING deleted !'
            end if
         end if
      end if

#ifdef _MPI
      call mpi_bcast(terminate,1,mpi_logical, 0,mpi_comm_world,mpi_er)
#endif

   end function terminate


!-------------------------------------------------------------------------
!  Initialisation
!-------------------------------------------------------------------------
   subroutine initialise()
      logical :: file_exist

      call mpi_precompute()
      if(mpi_rnk==0) call system('touch PRECOMPUTING')
      call system('echo $HOSTNAME > HOST')

      if(mpi_rnk==0)  print*, 'precomputing function requisites...'
      call par_precompute()
      call mes_precompute()
      call var_precompute()
      call tra_precompute()
      call tim_precompute()
      call vel_precompute()
      call temp_precompute()
      call non_precompute()
      call  io_precompute()
      call GTD_precompute()

      if (b_loadstate) then
          if(mpi_rnk==0)  print*, 'loading state...'
          tim_dt = 1d99
          call io_load_state()
      else
          if (mpi_rnk==0)  print*, 'Using default initial condition...'
          if (d_timestep>0) then
            tim_dt=d_timestep
          else
            tim_dt = d_maxdt/1d4
          end if
      end if
      call vel_matrices()
      call temp_matrices()

      if(mpi_rnk==0)  print*, 'initialising output files...'
      call io_openfiles()

      if(mpi_rnk==0) then
         open (99, file='PRECOMPUTING')
         close(99, status='delete')
         open (99, file='RUNNING')
         write(99,*) 'This file indicates a job running in this directory.'
         write(99,*) 'Delete this file to cleanly terminate the process.'
         close(99)
         print*, 'timestepping.....'
      end if

      call clk_time(d_start)
   end subroutine initialise

!-------------------------------------------------------------------------
!  Program closure
!-------------------------------------------------------------------------
   subroutine cleanup()
      logical :: file_exist

      if(mpi_rnk==0) then
         print*, 'cleanup...'
         call clk_time(d_stop)
         print*, ' sec/step  = ', (d_stop-d_start)/real(tim_step)
         print*, ' CPU time raw =', (d_stop-d_start)
#ifdef _CPUTIME
         print*, ' CPU time  = ', int((d_stop-d_start)/6d1), ' mins.'
#else
         print*, ' WALL time = ', int((d_stop-d_start)/6d1), ' mins.'
#endif
      end if

      call io_save_state()
      call io_closefiles()

      if(mpi_rnk==0) then
         inquire(file='RUNNING', exist=file_exist)
         if(file_exist) open(99, file='RUNNING')
         if(file_exist) close(99, status='delete')
      end if

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
      call mpi_finalize(mpi_er)
#endif
      call GTD_closing()
      if(mpi_rnk==0) print*, '...done!'

   end subroutine cleanup


!-------------------------------------------------------------------------
   subroutine clk_time(t)
      real, intent(out) :: t
#ifdef _CPUTIME
      call cpu_time(t)
#else
      integer, save :: ct,ctrt,ctmx,ct_=0,wrap=0
      call system_clock(ct,ctrt,ctmx)
      if(ct<ct_) wrap = wrap + 1
      ct_ = ct
      t = (real(ct)+real(ctmx)*wrap)/real(ctrt)
#endif
   end subroutine clk_time

!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************
