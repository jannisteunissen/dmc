! Simple diffusion Monte Carlo method
!
! Originally inspired by the demo codes from Alfonso Annarelli
!
! Author: Jannis Teunissen
program dmc
  use iso_fortran_env, only: error_unit, real32, real64, int32, int64
  use m_settings
  use m_config
  use m_mser

  implicit none

  integer, parameter :: dp = real64

  type walkers_t
     integer               :: n
     real(fp), allocatable :: w(:)
     real(fp), allocatable :: x(:, :, :)
     real(fp), allocatable :: phi(:)
     integer, allocatable  :: region(:)
     ! For updating population
     real(fp), allocatable :: x_new(:, :, :)
     real(fp), allocatable :: phi_new(:)
     integer, allocatable  :: region_new(:)
     real(dp), allocatable :: cdf(:)
     ! Random seed state
     integer(int64), allocatable :: s(:, :)
  end type walkers_t

  type(cfg_t)           :: cfg
  type(walkers_t)       :: walkers
  integer               :: i, n_steps, max_steps, steps_log
  real(dp)              :: dt, dt_step, time, end_time
  real(dp)              :: kappa
  real(dp)              :: sum_w, n_eff, n_eff_frac, trial_energy
  logical               :: mser_done
  real(dp)              :: sum_w_prev, growth_energy
  real(dp), allocatable :: growth_energy_array(:)
  integer               :: batch_size, trunc_obs, status
  real(dp)              :: min_mser
  logical, parameter    :: correct_autocorr = .true.
  integer(int64)        :: s_sequential(2)

  integer :: update_interval

  call cfg_add(cfg, "initial_distribution", "uniform", &
       "Initial walker distribution")
  call cfg_add(cfg, "dt", 1e-2_dp, &
       "Time step (atomic units)")
  call cfg_add(cfg, "end_time", 40.0_dp, &
       "End time (atomic units)")
  call cfg_add(cfg, "n_walkers", 1000, &
       "Number of walkers")
  call cfg_add(cfg, "initial_energy", 0.0_dp, &
       "Initial energy estimate (atomic units)")
  call cfg_add(cfg, "kappa", 2.0_dp, &
       "Coefficient for updating trial energy")
  call cfg_add(cfg, "n_eff_frac", 0.8_dp, &
       "Coefficient between 0 and 1 that determines when to resample")
  call cfg_add(cfg, "batch_size", 8, &
       "Batch size for estimate of local energy")
  call cfg_add(cfg, "update_interval", 10, &
       "How often to update trial energy (# steps)")

  call cfg_update_from_arguments(cfg)
  call cfg_check(cfg)

  call cfg_get(cfg, "end_time", end_time)
  call cfg_get(cfg, "dt", dt)
  call cfg_get(cfg, "initial_energy", trial_energy)
  trial_energy_fp = real(trial_energy, fp)
  call cfg_get(cfg, "kappa", kappa)
  call cfg_get(cfg, "n_eff_frac", n_eff_frac)
  call cfg_get(cfg, "batch_size", batch_size)
  call cfg_get(cfg, "update_interval", update_interval)

  sqrt_dt = real(sqrt(dt), fp)
  dt_step = dt * update_interval
  max_steps = ceiling(end_time/dt_step)
  steps_log = ceiling(max_steps*1e-2_dp)

  call walkers_initialize(cfg, walkers)

  allocate(growth_energy_array(max_steps))

  mser_done  = .false.
  time       = 0.0_dp
  sum_w_prev = sum(walkers%w)

  do n_steps = 1, max_steps

     do i = 1, update_interval
        call walkers_update(walkers, sqrt_dt, trial_energy_fp)
     end do

     call get_stats(walkers, sum_w, n_eff)

     growth_energy_array(n_steps) = trial_energy - log(sum_w/sum_w_prev)/dt_step
     sum_w_prev = sum_w

     if (mser_done) then
        growth_energy = growth_energy + &
             (growth_energy_array(n_steps) - growth_energy) / (n_steps - trunc_obs)
     end if

     if (n_eff < n_eff_frac * walkers%n) then
        call systematic_resample(walkers)
     end if

     if (mod(n_steps, steps_log) == 0) then
        if (.not. mser_done) then
           call mser_analysis(n_steps, growth_energy_array, batch_size, &
                correct_autocorr, trunc_obs, growth_energy, min_mser, status)
           mser_done = (trunc_obs > 0 .and. trunc_obs < n_steps/10)
        end if

        write(*, "(I4,'%',I8,' ',L,I8,F12.6,F12.6)") &
             ceiling((n_steps*1e2_dp)/max_steps), n_steps, &
             mser_done, trunc_obs, time, growth_energy
     end if

     time = time + dt_step
     trial_energy = growth_energy - log(sum_w/real(walkers%n, dp)) / (kappa*dt_step)
     trial_energy_fp = real(trial_energy, fp)
  end do

  call mser_analysis(n_steps-1, growth_energy_array, batch_size, &
       correct_autocorr, trunc_obs, growth_energy, min_mser, status)
  write(*, "(A,2F12.5)") " Final result:  ", growth_energy, sqrt(min_mser)
  write(*, "(A,E12.4)") " Total updates: ", n_steps * real(walkers%n, dp)

contains

  subroutine walkers_initialize(cfg, walkers)
    type(cfg_t), intent(inout)     :: cfg
    type(walkers_t), intent(inout) :: walkers
    character(len=20)              :: initial_distribution
    integer                        :: i, n, idim
    integer(int64)                 :: initial_seed(2), rand_int64
    real(dp)                       :: r(2)

    call cfg_get(cfg, "n_walkers", walkers%n)

    allocate(walkers%w(walkers%n))
    allocate(walkers%x(n_dim, n_particles, walkers%n))
    allocate(walkers%phi(walkers%n))
    allocate(walkers%x_new(n_dim, n_particles, walkers%n))
    allocate(walkers%phi_new(walkers%n))
    allocate(walkers%cdf(walkers%n))
    allocate(walkers%s(2, walkers%n))
    allocate(walkers%region(walkers%n))
    allocate(walkers%region_new(walkers%n))

    ! Set initial seeds
    call random_seed()

    call random_number(r)
    initial_seed = transfer(r, initial_seed)

    call random_number(r)
    s_sequential = transfer(r, initial_seed)

    walkers%s(:, 1) = initial_seed

    do i = 2, walkers%n
       walkers%s(:, i) = walkers%s(:, i-1)
       call jump(walkers%s(1, i), walkers%s(2, i))
    end do

    call cfg_get(cfg, "initial_distribution", initial_distribution)

    select case (initial_distribution)
    case ("uniform")
       do i = 1, walkers%n
          do n = 1, n_particles
             do idim = 1, n_dim
                call next(walkers%s(1, i), walkers%s(2, i), rand_int64)
                walkers%x(idim, n, i) = real(1.0_dp * (uni01_64(rand_int64) - 0.5_dp), fp)
             end do
          end do
       end do
    case default
       error stop "Invalid initial_distribution. Options: uniform"
    end select

    ! Compute initial potential and set weight
    do i = 1, walkers%n
       call compute_potential(walkers%x(:, :, i), walkers%phi(i))
       walkers%w(i) = 1.0_dp
       walkers%region(i) = spin_region(walkers%x(:, :, i))
    end do

    ! OpenACC: store on device
    !$acc enter data copyin(walkers)
    !$acc enter data copyin(walkers%w)
    !$acc enter data copyin(walkers%x, walkers%x_new)
    !$acc enter data copyin(walkers%phi, walkers%phi_new)
    !$acc enter data copyin(walkers%region, walkers%region_new)
    !$acc enter data copyin(walkers%cdf)
    !$acc enter data copyin(walkers%s)

  end subroutine walkers_initialize

  subroutine walkers_update(walkers, sqrt_dt, trial_energy)
    type(walkers_t), intent(inout) :: walkers
    real(fp), intent(in)           :: sqrt_dt, trial_energy
    integer                        :: i, n, idim, ix, old_region
    real(fp)                       :: phi_avg, exp_arg
    real(fp)                       :: r1, r2

    !$omp parallel do private(phi_avg, n, r1, r2, idim, ix, exp_arg, old_region)
    !$acc parallel loop private(phi_avg, n, r1, r2, idim, ix, exp_arg, old_region)
    do i = 1, walkers%n
       phi_avg = walkers%phi(i)
       old_region = walkers%region(i)

       do n = 1, n_particles
          do idim = 1, n_dim, 2
             call box_muller_32(walkers%s(1, i), walkers%s(2, i), r1, r2)
             walkers%x(idim, n, i) = walkers%x(idim, n, i) + sqrt_dt * r1
             if (idim + 1 <= n_dim) then
                walkers%x(idim+1, n, i) = walkers%x(idim+1, n, i) + sqrt_dt * r2
             end if
          end do
       end do

       walkers%region(i) = spin_region(walkers%x(:, :, i))

       if (walkers%region(i) /= old_region) then
          ! Kill walker
          walkers%w(i) = 0
       else
          ! Compute new potential
          call compute_potential(walkers%x(:, :, i), walkers%phi(i))
       end if

       phi_avg = 0.5_fp * (phi_avg + walkers%phi(i))

       exp_arg = real(-dt * (phi_avg - trial_energy), fp)
       exp_arg = min(exp_arg, 3.0_fp)
       walkers%w(i) = walkers%w(i) * exp(exp_arg)
    end do
    !$acc end parallel loop
  end subroutine walkers_update

  subroutine get_stats(walkers, sum_w, n_eff)
    type(walkers_t), intent(inout) :: walkers
    real(dp), intent(out)          :: sum_w
    real(dp), intent(out)          :: n_eff
    real(dp)                       :: sum_w2

    sum_w  = 0.0_dp
    sum_w2 = 0.0_dp

    !$omp parallel do reduction(+: sum_w, sum_w2)
    !$acc parallel loop reduction(+: sum_w, sum_w2) default(present)
    do i = 1, walkers%n
       sum_w  = sum_w  + walkers%w(i)
       sum_w2 = sum_w2 + walkers%w(i)**2
    end do
    !$acc end parallel loop

    n_eff = sum_w**2 / sum_w2
  end subroutine get_stats

  subroutine systematic_resample(wlk)
    type(walkers_t), intent(inout) :: wlk
    integer                        :: i, j, n
    integer(int64)                 :: rand_int64
    real(dp)                       :: u0, step, mean_w, pos, prefix_sum

    n = walkers%n

    !$acc update self(wlk%w)

    prefix_sum = 0.0_dp
    !$omp parallel do reduction(inscan, +:prefix_sum)
    do i = 1, n
       prefix_sum = prefix_sum + wlk%w(i)
       !$omp scan inclusive(prefix_sum)
       wlk%cdf(i) = prefix_sum
    end do

    !$acc update device(wlk%cdf)
    step = wlk%cdf(n) / n

    call next(s_sequential(1), s_sequential(2), rand_int64)
    u0 = uni01_64(rand_int64) * step

    ! Resample: pick ancestors via search on the CDF
    !$omp parallel do private(pos, j)
    !$acc parallel loop private(pos, j) default(present)
    do i = 1, n
       pos = u0 + (i-1) * step
       j = upper_bound(n, wlk%cdf, pos) ! first index with cdf(j) >= pos

       wlk%x_new(:, :, i) = wlk%x(:, :, j)
       wlk%phi_new(i)     = wlk%phi(j)
       wlk%region_new(i)  = wlk%region(j)
    end do
    !$acc end parallel loop

    mean_w = step

    ! Copy back
    !$omp parallel do
    !$acc parallel loop default(present)
    do i = 1, n
       wlk%x(:, :, i)  = wlk%x_new(:, :, i)
       wlk%phi(i)      = wlk%phi_new(i)
       wlk%region(i)   = wlk%region_new(i)
       wlk%w(i)        = real(mean_w, fp)
    end do
    !$acc end parallel loop
  end subroutine systematic_resample

  pure integer function upper_bound(n, cdf, val) result(lo)
    !$acc routine seq
    integer, intent(in)  :: n
    real(dp), intent(in) :: cdf(n), val
    integer              :: hi, mid

    lo = 1
    hi = n
    do while (lo < hi)
       mid = (lo + hi) / 2
       if (cdf(mid) < val) then
          lo = mid + 1
       else
          hi = mid
       end if
    end do
  end function upper_bound

  pure subroutine compute_potential(x, phi)
    !$acc routine seq
    real(fp), intent(in)  :: x(n_dim, n_particles)
    real(fp), intent(out) :: phi
    integer               :: n, m

    phi = 0.0_fp

#if defined(POT_ATOM)
    ! Electron-nucleus terms
    do n = 1, n_particles
       phi = phi - atom_z / norm2(x(:, n))
    end do

    ! Electron-electron terms
    do n = 1, n_particles
       do m = n+1, n_particles
          phi = phi + 1.0_fp / norm2(x(:, n) - x(:, m))
       end do
    end do
#elif defined(POT_HARMONIC)
    do n = 1, n_particles
       phi = phi + 0.5_fp * sum(x(:, n)**2)
    end do
#else
#error "No potential selected (define POT_ATOM or POT_HARMONIC)"
#endif
  end subroutine compute_potential

  pure integer function spin_region(x) result(r)
    !$acc routine seq
    real(fp), intent(in) :: x(n_dim, n_particles)
    real(fp)             :: s

#if defined(REGION_HO2D)
    s = norm2(x(:, 3)) - 2.0_fp/3.0_fp
#else if defined(REGION_NONE)
    s = 1.0_fp
#endif
    r = merge(1, -1, s >= 0.0_fp)
  end function spin_region

  include 'rng.f90'

end program dmc
