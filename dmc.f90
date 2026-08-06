! Simple diffusion Monte Carlo method, inspired by the demo codes from Alfonso
! Annarelli, alfonso.annarelli@gmail.com
!
! Author: Jannis Teunissen
program qmc
  use iso_fortran_env, only: error_unit, real32, real64, int32, int64
  use m_config
  use m_mser

  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: fp = real32

  integer, parameter :: potential_atom     = 1
  integer, parameter :: potential_harmonic = 2
  integer            :: potential_type

  type walkers_t
     integer               :: n
     integer               :: n_dim
     integer               :: n_particles
     integer               :: n_spin_up
     integer               :: n_spin_down
     real(dp)              :: n_eff_frac = 0.8_dp
     real(fp), allocatable :: w(:)
     real(fp), allocatable :: x(:, :, :)
     real(fp), allocatable :: phi(:)
     integer, allocatable  :: regions(:, :)
     ! For updating population
     real(fp), allocatable :: x_new(:, :, :)
     real(fp), allocatable :: phi_new(:)
     integer, allocatable  :: regions_new(:, :)
     real(dp), allocatable :: cdf(:)
     ! Random seed state
     integer(int64), allocatable :: s(:, :)
  end type walkers_t

  type(cfg_t)       :: cfg
  type(walkers_t)   :: walkers
  integer           :: n_steps, max_steps
  real(dp)          :: dt, time, end_time
  real(dp)          :: trial_energy
  real(dp)          :: kappa
  real(dp)          :: sum_w, n_eff, tmp_float64
  character(len=40) :: potential_name
  logical           :: mser_done

  real(dp)              :: sum_w_prev, growth_energy
  real(dp), allocatable :: growth_energy_array(:)

  real(fp) :: atom_z
  integer  :: batch_size, trunc_obs, status
  real(dp) :: min_mser
  logical, parameter :: correct_autocorr = .true.

  call cfg_add(cfg, "potential", "atom", &
       "Type of potential")
  call cfg_add(cfg, "atom_z", 1.0_dp, &
       "Atomic number")
  call cfg_add(cfg, "initial_distribution", "uniform", &
       "Initial walker distribution")
  call cfg_add(cfg, "n_dim", 3, &
       "Number of dimensions for wave function")
  call cfg_add(cfg, "n_particles", 1, &
       "Number of particles")
  call cfg_add(cfg, "n_spin_up", 0, &
       "Number of particles with spin up")
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
  call cfg_add(cfg, "verbose", 1, &
       "Verbosity (0: silent, 1: info)")
  call cfg_add(cfg, "batch_size", 8, &
       "Batch size for estimate of local energy")

  call cfg_update_from_arguments(cfg)
  call cfg_check(cfg)

  call cfg_get(cfg, "end_time", end_time)
  call cfg_get(cfg, "dt", dt)
  call cfg_get(cfg, "initial_energy", trial_energy)
  call cfg_get(cfg, "kappa", kappa)
  call cfg_get(cfg, "batch_size", batch_size)

  call cfg_get(cfg, "potential", potential_name)
  call cfg_get(cfg, "atom_z", tmp_float64)
  atom_z = real(tmp_float64, fp)

  select case (potential_name)
  case ("atom")
     potential_type = potential_atom
  case ("harmonic")
     potential_type = potential_harmonic
  case default
     error stop "Unknown potential type"
  end select

  max_steps = ceiling(end_time/dt)

  call walkers_initialize(cfg, walkers)

  allocate(growth_energy_array(max_steps))

  mser_done = .false.
  time = 0.0_dp
  sum_w_prev = sum(walkers%w)

  do n_steps = 1, max_steps
     call walkers_update(walkers, dt, trial_energy, sum_w, n_eff)

     growth_energy_array(n_steps) = trial_energy - log(sum_w/sum_w_prev)/dt
     sum_w_prev = sum_w

     if (mser_done) then
        growth_energy = growth_energy + &
             (growth_energy_array(n_steps) - growth_energy) / (n_steps - trunc_obs)
     end if

     if (n_eff < walkers%n_eff_frac * walkers%n) then
        call systematic_resample(walkers, sum_w)
     end if

     if (mod(n_steps, max_steps/100) == 0) then
        if (.not. mser_done) then
           call mser_analysis(n_steps, growth_energy_array, batch_size, &
                correct_autocorr, trunc_obs, growth_energy, min_mser, status)
           mser_done = (trunc_obs > 0 .and. trunc_obs < n_steps/10)
        end if

        print *, mser_done, n_steps, trunc_obs, time, growth_energy
     end if

     time = time + dt
     trial_energy = growth_energy - log(sum_w/real(walkers%n, dp)) / (kappa*dt)
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

    call cfg_get(cfg, "n_walkers", walkers%n)
    call cfg_get(cfg, "n_dim", walkers%n_dim)
    call cfg_get(cfg, "n_particles", walkers%n_particles)
    call cfg_get(cfg, "n_spin_up", walkers%n_spin_up)
    walkers%n_spin_down = walkers%n_particles - walkers%n_spin_up

    allocate(walkers%w(walkers%n))
    allocate(walkers%x(walkers%n_dim, walkers%n_particles, walkers%n))
    allocate(walkers%phi(walkers%n))
    allocate(walkers%x_new(walkers%n_dim, walkers%n_particles, walkers%n))
    allocate(walkers%phi_new(walkers%n))
    allocate(walkers%cdf(walkers%n))
    allocate(walkers%s(2, walkers%n))
    allocate(walkers%regions(2, walkers%n))
    allocate(walkers%regions_new(2, walkers%n))

    ! Set initial seeds
    initial_seed = [1234_int64, -89234_int64]

    walkers%s(:, 1) = initial_seed

    do i = 2, walkers%n
       walkers%s(:, i) = walkers%s(:, i-1)
       call jump(walkers%s(1, i), walkers%s(2, i))
    end do

    call cfg_get(cfg, "initial_distribution", initial_distribution)

    select case (initial_distribution)
    case ("uniform")
       do i = 1, walkers%n
          do n = 1, walkers%n_particles
             do idim = 1, walkers%n_dim
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
       call compute_potential(walkers, i)
       walkers%w(i) = 1.0_dp
       call set_regions(walkers, i)
    end do

  end subroutine walkers_initialize

  subroutine walkers_update(walkers, dt, trial_energy, sum_w, n_eff)
    type(walkers_t), intent(inout) :: walkers
    real(dp), intent(in)           :: dt, trial_energy
    real(dp), intent(out)          :: sum_w
    real(dp), intent(out)          :: n_eff
    integer                        :: i, n, idim, ix, old_regions(2)
    real(fp)                       :: sqrt_dt, phi_avg, exp_arg
    real(dp)                       :: sum_w2
    real(fp)                       :: rr(walkers%n_dim, walkers%n_particles)

    sqrt_dt = sqrt(real(dt, fp))
    sum_w  = 0.0_dp
    sum_w2 = 0.0_dp

    !$omp parallel do private(phi_avg, n, rr, idim, ix, exp_arg, old_regions) &
    !$omp &reduction(+: sum_w, sum_w2)
    do i = 1, walkers%n
       phi_avg = walkers%phi(i)
       old_regions = walkers%regions(:, i)

       call get_normal_numbers(walkers%n_dim, walkers%n_particles, &
            walkers%s(1, i), walkers%s(2, i), rr)

       walkers%x(:, :, i) = walkers%x(:, :, i) + sqrt_dt * rr

       call set_regions(walkers, i)

       if (any(walkers%regions(:, i) /= old_regions)) then
          ! Kill walker
          walkers%w(i) = 0
       else
          ! Compute new potential
          call compute_potential(walkers, i)
       end if

       phi_avg = 0.5_fp * (phi_avg + walkers%phi(i))

       exp_arg = real(-dt * (phi_avg - trial_energy), fp)
       exp_arg = min(exp_arg, 5.0_fp)
       walkers%w(i) = walkers%w(i) * exp(exp_arg)

       sum_w  = sum_w  + walkers%w(i)
       sum_w2 = sum_w2 + walkers%w(i)**2
    end do

    n_eff = sum_w**2 / sum_w2
  end subroutine walkers_update

  subroutine set_regions(walkers, i)
    type(walkers_t), intent(inout) :: walkers
    integer, intent(in)            :: i

    associate (n_up => walkers%n_spin_up, n_down => walkers%n_spin_down, &
         n_dim => walkers%n_dim, regions => walkers%regions(:, i))
      regions(:) = 0
      if (n_up > 0) then
         regions(1) = spin_up_region(n_up, n_dim, walkers%x(:, 1:n_up, i))
      end if
      if (n_down > 0) then
         regions(2) = spin_down_region(n_down, n_dim, walkers%x(:, n_up+1:, i))
      end if
    end associate
  end subroutine set_regions

  integer function spin_up_region(n, ndim, x)
    integer, intent(in) :: n, ndim
    real(fp)            :: x(ndim, n)
    spin_up_region = spin_region_ho_2d(n, ndim, x)
  end function spin_up_region

  integer function spin_down_region(n, ndim, x)
    integer, intent(in) :: n, ndim
    real(fp)            :: x(ndim, n)
    spin_down_region = spin_region_ho_2d(n, ndim, x)
  end function spin_down_region

  integer function spin_region_ho_2d(n, ndim, x)
    integer, intent(in) :: n, ndim
    real(fp), intent(in):: x(ndim, n)
    real(fp)            :: s

    s = 1
    ! s = x(1, 1)
    ! s = norm2(x(:, 1))
    ! s = (s - 0.25_fp) * (s - 0.5_fp)

    if (s < 0.0_fp) then
       spin_region_ho_2d = -1
    else
       spin_region_ho_2d = 1
    end if
  end function spin_region_ho_2d

  subroutine systematic_resample(walkers, sum_w)
    type(walkers_t), intent(inout) :: walkers
    real(dp), intent(in)           :: sum_w
    integer                        :: i, j, n
    integer(int64)                 :: rand_int64
    real(dp)                       :: u0, step, mean_w, pos, prefix_sum

    n = walkers%n

    associate (cdf => walkers%cdf, x_new => walkers%x_new, &
         phi_new => walkers%phi_new, regions_new => walkers%regions_new)

      prefix_sum = 0.0
      !$omp parallel do reduction(inscan, +:prefix_sum)
      do i = 1, n
         prefix_sum = prefix_sum + walkers%w(i)
         !$omp scan inclusive(prefix_sum)
         cdf(i) = prefix_sum
      end do

      step = sum_w / n

      call next(walkers%s(1, 1), walkers%s(2, 1), rand_int64)
      u0 = uni01_64(rand_int64) * step        ! single random offset

      do i = 1, n
         pos = u0 + (i-1) * step
         j   = upper_bound(cdf, pos)      ! first index with cdf(j) >= pos
         x_new(:, :, i) = walkers%x(:, :, j)
         phi_new(i)  = walkers%phi(j)
         regions_new(:, i) = walkers%regions(:, j)
      end do

      mean_w = sum_w / n
      do i = 1, n
         walkers%x(:, :, i) = x_new(:, :, i)
         walkers%phi(i)  = phi_new(i)
         walkers%regions(:, i)  = regions_new(:, i)
         walkers%w(i)    = real(mean_w, fp)
      end do
    end associate
  end subroutine systematic_resample

  pure integer function upper_bound(cdf, val) result(lo)
    !$acc routine seq
    real(dp), intent(in) :: cdf(:), val
    integer              :: hi, mid
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

  pure subroutine compute_potential(walkers, i)
    type(walkers_t), intent(inout) :: walkers
    integer, intent(in)            :: i

    select case (potential_type)
    case (potential_atom)
       call compute_potential_atom(walkers, i)
    case (potential_harmonic)
       call compute_potential_harmonic(walkers, i)
    end select
  end subroutine compute_potential

  pure subroutine compute_potential_atom(walkers, i)
    type(walkers_t), intent(inout) :: walkers
    integer, intent(in)            :: i
    integer                        :: n, m

    walkers%phi(i) = 0.0_dp

    ! Electron-nucleus terms
    do n = 1, walkers%n_particles
       walkers%phi(i) = walkers%phi(i) - atom_z/norm2(walkers%x(:, n, i))
    end do

    ! Electron-electron terms
    do n = 1, walkers%n_particles
       do m = n+1, walkers%n_particles
          walkers%phi(i) = walkers%phi(i) + &
               1/norm2(walkers%x(:, n, i) - walkers%x(:, m, i))
       end do
    end do
  end subroutine compute_potential_atom

  pure subroutine compute_potential_harmonic(walkers, i)
    type(walkers_t), intent(inout) :: walkers
    integer, intent(in)            :: i
    integer                        :: n

    walkers%phi(i) = 0.0_dp

    do n = 1, walkers%n_particles
       walkers%phi(i) = walkers%phi(i) + 0.5_fp * sum(walkers%x(:, n, i)**2)
    end do
  end subroutine compute_potential_harmonic

  subroutine get_normal_numbers(nx, ny, s1, s2, rr)
    integer, intent(in)           :: nx, ny
    integer(int64), intent(inout) :: s1, s2
    real(fp), intent(out)         :: rr(nx, ny)
    real(fp)                      :: tmp(nx*ny+1)
    integer                       :: i

    do i = 1, size(tmp)/2
       call box_muller_32(s1, s2, tmp(2*i-1), tmp(2*i))
    end do
    rr = reshape(tmp(1:nx*ny), [nx, ny])
  end subroutine get_normal_numbers

  include 'rng.f90'

end program qmc
