! Simple diffusion Monte Carlo method, inspired by the demo codes from Alfonso
! Annarelli, alfonso.annarelli@gmail.com
!
! Author: Jannis Teunissen
program qmc
  use iso_fortran_env, only: error_unit, real64
  use m_config
  use m_random
  use m_mser

  implicit none

  integer, parameter  :: dp = real64

  type walkers_t
     integer               :: n
     integer               :: ndim
     integer               :: branch_limit
     real(dp)              :: n_eff_frac = 0.8_dp
     real(dp), allocatable :: w(:)
     real(dp), allocatable :: x(:, :)
     real(dp), allocatable :: phi(:)
     ! For updating population
     real(dp), allocatable :: x_new(:, :), phi_new(:), cdf(:)
  end type walkers_t

  type(cfg_t)           :: cfg
  type(rng_t)           :: rng
  type(walkers_t)       :: walkers
  integer               :: n_steps, max_steps
  real(dp)              :: dt, time, end_time
  real(dp)              :: trial_energy
  real(dp)              :: mean_local_energy
  real(dp), allocatable :: mean_local_energy_array(:)
  real(dp)              :: ratio, kappa
  real(dp)              :: sum_w, n_eff

  integer  :: batch_size, trunc_obs, status
  real(dp) :: min_mser
  logical, parameter  :: correct_autocorr = .true.

  call cfg_add(cfg, "potential_type", "harmonic", &
       "Type of potential")
  call cfg_add(cfg, "initial_distribution", "uniform", &
       "Initial walker distribution")
  call cfg_add(cfg, "ndim", 3, &
       "Number of dimensions for wave function (D * n_particles)")
  call cfg_add(cfg, "dt", 5e-3_dp, &
       "Time step (atomic units)")
  call cfg_add(cfg, "end_time", 20.0_dp, &
       "End time (atomic units)")
  call cfg_add(cfg, "num_walkers", 100000, &
       "Target number of walkers")
  call cfg_add(cfg, "walkers_branch_limit", 3, &
       "Maximum branching ratio (number of copies) per walker")
  call cfg_add(cfg, "initial_energy", 0.0_dp, &
       "Initial energy estimate (atomic units)")
  call cfg_add(cfg, "kappa", 2.0_dp, &
       "Coefficient for updating trial energy")
  call cfg_add(cfg, "verbose", 1, &
       "Verbosity (0: silent, 1: info)")
  call cfg_add(cfg, "batch_size", 8, &
       "Batch size for estimate of local energy")

  call cfg_update_from_arguments(cfg)

  call cfg_get(cfg, "end_time", end_time)
  call cfg_get(cfg, "dt", dt)
  call cfg_get(cfg, "initial_energy", trial_energy)
  call cfg_get(cfg, "kappa", kappa)
  call cfg_get(cfg, "batch_size", batch_size)

  max_steps = ceiling(end_time/dt)

  call rng%set_random_seed()

  call walkers_initialize(cfg, walkers)

  allocate(mean_local_energy_array(max_steps))

  time = 0.0_dp
  n_steps = 0
  mean_local_energy = 0.0_dp

  do n_steps = 1, max_steps
     call walkers_update(walkers, dt, trial_energy, &
          mean_local_energy_array(n_steps), sum_w, n_eff)

     if (n_eff < walkers%n_eff_frac * walkers%n) then
        call systematic_resample(walkers, rng, sum_w)
     end if

     if (mod(n_steps, max_steps/100) == 0) then
        call mser_analysis(n_steps, mean_local_energy_array, batch_size, &
             correct_autocorr, trunc_obs, mean_local_energy, min_mser, status)
        print *, "AVG", n_steps, time, trunc_obs, mean_local_energy, sqrt(min_mser)
     end if

     time = time + dt
     ratio = sum_w / walkers%n
     trial_energy = mean_local_energy - 1/(kappa*dt) * log(ratio)
  end do

  write(*, "(A,E12.4)") "Total updates: ", n_steps * real(walkers%n, dp)

contains

  subroutine walkers_initialize(cfg, walkers)
    type(cfg_t), intent(inout)     :: cfg
    type(walkers_t), intent(inout) :: walkers
    character(len=20)              :: initial_distribution
    integer                        :: i, idim

    call cfg_get(cfg, "num_walkers", walkers%n)
    call cfg_get(cfg, "ndim", walkers%ndim)
    call cfg_get(cfg, "walkers_branch_limit", walkers%branch_limit)

    allocate(walkers%w(walkers%n))
    allocate(walkers%x(walkers%ndim, walkers%n))
    allocate(walkers%phi(walkers%n))

    allocate(walkers%x_new(walkers%ndim, walkers%n))
    allocate(walkers%phi_new(walkers%n))
    allocate(walkers%cdf(walkers%n))

    call cfg_get(cfg, "initial_distribution", initial_distribution)

    select case (initial_distribution)
    case ("uniform")
       do i = 1, walkers%n
          do idim = 1, walkers%ndim
             walkers%x(idim, i) = rng%unif_01() - 0.5_dp
          end do
       end do
    case default
       error stop "Invalid initial_distribution. Options: uniform"
    end select

    ! Compute initial potential and set weight
    do i = 1, walkers%n
       call compute_potential(walkers, i)
       walkers%w(i) = 1.0_dp
    end do

  end subroutine walkers_initialize

  subroutine walkers_update(walkers, dt, trial_energy, local_energy, &
       sum_w, n_eff)
    type(walkers_t), intent(inout) :: walkers
    real(dp), intent(in)           :: dt, trial_energy
    real(dp), intent(out)          :: local_energy
    real(dp), intent(out)          :: sum_w
    real(dp), intent(out)          :: n_eff
    integer                        :: i, idim
    real(dp)                       :: phi_old, phi_avg, sqrt_dt
    real(dp)                       :: sum_w2, sum_we

    sqrt_dt = sqrt(dt)
    sum_w  = 0.0_dp
    sum_we = 0.0_dp
    sum_w2 = 0.0_dp

    do i = 1, walkers%n
       phi_old = walkers%phi(i)
       do idim = 1, walkers%ndim
          walkers%x(idim, i) = walkers%x(idim, i) + sqrt_dt * rng%normal()
       end do
       call compute_potential(walkers, i)
       phi_avg = 0.5_dp * (walkers%phi(i) + phi_old)

       walkers%w(i) = walkers%w(i) * &
            exp(-dt * (phi_avg - trial_energy))

       sum_w  = sum_w  + walkers%w(i)
       sum_we = sum_we + walkers%w(i) * phi_avg
       sum_w2 = sum_w2 + walkers%w(i)**2
    end do

    local_energy = sum_we / sum_w
    n_eff = sum_w**2 / sum_w2
  end subroutine walkers_update

  subroutine systematic_resample(walkers, rng, sum_w)
    type(walkers_t), intent(inout) :: walkers
    type(rng_t), intent(inout)     :: rng
    real(dp), intent(in)           :: sum_w
    integer                        :: i, j, n
    real(dp)                       :: u0, step, mean_w, pos

    n = walkers%n

    associate (cdf => walkers%cdf, x_new => walkers%x_new, phi_new => walkers%phi_new)
      cdf(1) = walkers%w(1)
      do i = 2, n
         cdf(i) = cdf(i-1) + walkers%w(i)
      end do

      step = sum_w / n
      u0   = rng%unif_01() * step        ! single random offset

      ! --- gather step: fully parallel, each thread independent ---
      !$acc data copyin(cdf, walkers%x, walkers%phi) &
      !$acc      copyout(x_new, phi_new)
      !$acc parallel loop private(pos, j)
      do i = 1, n
         pos = u0 + (i-1) * step
         j   = upper_bound(cdf, pos)      ! first index with cdf(j) >= pos
         x_new(:, i) = walkers%x(:, j)
         phi_new(i)  = walkers%phi(j)     ! copy phi, do NOT recompute
      end do

      mean_w = sum_w / n
      do i = 1, n
         walkers%x(:, i) = x_new(:, i)
         walkers%phi(i)  = phi_new(i)
         walkers%w(i)    = mean_w
      end do
    end associate
  end subroutine systematic_resample

  pure integer function upper_bound(cdf, val) result(lo)
    !$acc routine seq
    real(dp), intent(in) :: cdf(:), val
    integer :: hi, mid
    lo = 1; hi = size(cdf)
    do while (lo < hi)
       mid = (lo + hi) / 2
       if (cdf(mid) < val) then
          lo = mid + 1
       else
          hi = mid
       end if
    end do
  end function upper_bound

  subroutine compute_potential(walkers, i)
    type(walkers_t), intent(inout) :: walkers
    integer, intent(in)            :: i

    ! walkers%phi(i) = 0.5_dp * sum(walkers%x(:, i)**2)
    walkers%phi(i) = -1.0_dp / norm2(walkers%x(:, i))
  end subroutine compute_potential

end program qmc
