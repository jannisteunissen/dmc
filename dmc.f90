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
     integer               :: n_max
     integer               :: n_target
     integer               :: ndim
     integer               :: branch_limit
     integer, allocatable  :: w(:)
     real(dp), allocatable :: x(:, :)
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
  real(dp), allocatable :: potential(:, :)

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
  call cfg_add(cfg, "num_walkers_target", 100000, &
       "Target number of walkers")
  call cfg_add(cfg, "max_walkers_ratio", 4.0_dp, &
       "Maximum ratio of num_walkers/num_walkers_target")
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

  allocate(potential(walkers%n_max, 2))
  allocate(mean_local_energy_array(max_steps))

  time = 0.0_dp
  n_steps = 0
  mean_local_energy = 0.0_dp

  do n_steps = 1, max_steps
     ! Compute potential at current position
     call compute_potential(walkers, potential(:, 1))

     call walkers_diffuse(walkers, dt)

     ! Compute potential at new position
     call compute_potential(walkers, potential(:, 2))

     call walkers_branch(walkers, dt, potential, trial_energy, &
          mean_local_energy_array(n_steps))

     if (mod(n_steps, max_steps/100) == 0) then
        call mser_analysis(n_steps, mean_local_energy_array, batch_size, &
             correct_autocorr, trunc_obs, mean_local_energy, min_mser, status)
        print *, "AVG", n_steps, time, trunc_obs, mean_local_energy, sqrt(min_mser)
     end if

     time = time + dt

     ratio = walkers%n / real(walkers%n_target, dp)
     trial_energy = mean_local_energy - 1/(kappa*dt) * log(ratio)
  end do

  write(*, "(A,E12.4)") "Total updates: ", n_steps * real(walkers%n_target, dp)

contains

  subroutine walkers_initialize(cfg, walkers)
    type(cfg_t), intent(inout)     :: cfg
    type(walkers_t), intent(inout) :: walkers
    character(len=20)              :: initial_distribution
    integer                        :: i, idim
    real(dp)                       :: max_walkers_ratio

    call cfg_get(cfg, "num_walkers_target", walkers%n_target)
    call cfg_get(cfg, "max_walkers_ratio", max_walkers_ratio)
    call cfg_get(cfg, "ndim", walkers%ndim)
    call cfg_get(cfg, "walkers_branch_limit", walkers%branch_limit)

    walkers%n = walkers%n_target
    walkers%n_max = ceiling(max_walkers_ratio * walkers%n_target)
    allocate(walkers%w(walkers%n_max))
    allocate(walkers%x(walkers%ndim, walkers%n_max))

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

  end subroutine walkers_initialize

  subroutine walkers_diffuse(walkers, dt)
    type(walkers_t), intent(inout) :: walkers
    real(dp), intent(in)           :: dt
    integer                        :: i, idim
    real(dp)                       :: sqrt_dt

    sqrt_dt = sqrt(dt)
    do i = 1, walkers%n
       do idim = 1, walkers%ndim
          ! Note: rng%normal() is not thread-safe
          walkers%x(idim, i) = walkers%x(idim, i) + sqrt_dt * rng%normal()
       end do
    end do
  end subroutine walkers_diffuse

  subroutine walkers_branch(walkers, dt, potential, trial_energy, local_energy)
    type(walkers_t), intent(inout) :: walkers
    real(dp), intent(in)           :: dt
    real(dp), intent(in)           :: potential(walkers%n_max, 2)
    real(dp), intent(in)           :: trial_energy
    real(dp), intent(out)          :: local_energy
    integer                        :: i, j, n_old, ix_tail
    real(dp)                       :: phi, w, exp_arg

    ! First pass: compute weights
    local_energy = 0.0_dp
    n_old = walkers%n
    do i = 1, n_old
       phi = 0.5_dp * (potential(i, 1) + potential(i, 2))
       exp_arg = min(-dt * (phi - trial_energy), 20.0_dp)
       w = exp(exp_arg) + rng%unif_01()
       walkers%w(i) = min(int(w), walkers%branch_limit)
       local_energy = local_energy + phi
    end do

    local_energy = local_energy / n_old

    ! Place extra copies (w > 1) at the tail
    ix_tail = n_old
    do i = 1, n_old
       do j = 2, walkers%w(i)  ! only extra copies
          ix_tail = ix_tail + 1
          walkers%x(:, ix_tail) = walkers%x(:, i)
          walkers%w(ix_tail) = 1  ! mark copies as valid
       end do
    end do

    ! Compact: remove walkers with w=0 by swapping from tail
    i = 1
    do while (i <= ix_tail)
       if (walkers%w(i) == 0) then
          ! Skip tail walkers that are also dead
          do while (ix_tail > i .and. walkers%w(ix_tail) == 0)
             ix_tail = ix_tail - 1
          end do
          if (i >= ix_tail) then
             ix_tail = i - 1
             exit
          end if
          walkers%x(:, i) = walkers%x(:, ix_tail)
          walkers%w(i) = 1
          ix_tail = ix_tail - 1
       end if
       i = i + 1
    end do

    walkers%n = ix_tail
  end subroutine walkers_branch

  ! subroutine walkers_branch(walkers, dt, potential, trial_energy, local_energy)
  !   type(walkers_t), intent(inout) :: walkers
  !   real(dp), intent(in)           :: dt
  !   real(dp), intent(in)           :: potential(walkers%n_max, 2)
  !   real(dp), intent(in)           :: trial_energy
  !   real(dp), intent(out)          :: local_energy
  !   integer                        :: i, j, n_old, n_new, offset
  !   real(dp)                       :: phi, w, exp_arg
  !   integer, allocatable           :: counts(:)   ! offspring count per walker
  !   integer, allocatable           :: offsets(:)  ! exclusive prefix sum
  !   real(dp), allocatable          :: x_new(:, :) ! destination buffer

  !   n_old = walkers%n
  !   allocate(counts(n_old), offsets(n_old))

  !   ! -------- Pass 1: parallel MAP (weights + local energy) --------
  !   ! Each iteration independent -> GPU: one thread per walker.
  !   ! local_energy accumulation -> GPU: parallel REDUCTION.
  !   local_energy = 0.0_dp
  !   do i = 1, n_old
  !      phi     = 0.5_dp * (potential(i, 1) + potential(i, 2))
  !      exp_arg = min(-dt * (phi - trial_energy), 20.0_dp)
  !      w       = exp(exp_arg) + rng%unif_01()   ! use per-walker RNG stream on GPU
  !      counts(i) = min(int(w), walkers%branch_limit)
  !      local_energy = local_energy + phi
  !   end do
  !   local_energy = local_energy / n_old

  !   ! -------- Pass 2: exclusive SCAN (prefix sum of counts) --------
  !   ! GPU: CUB/Thrust exclusive_scan. Serial reference below.
  !   offset = 0
  !   do i = 1, n_old
  !      offsets(i) = offset          ! 0-based write position for walker i
  !      offset     = offset + counts(i)
  !   end do
  !   n_new = offset                  ! total survivors = last offset + last count

  !   if (n_new > walkers%n_max) then
  !      error stop "walker overflow: increase max_walkers_ratio"
  !   end if

  !   ! -------- Pass 3: SCATTER (write offspring to contiguous slots) --------
  !   ! Each source walker writes counts(i) copies into disjoint positions,
  !   ! so writes never collide -> fully parallel.
  !   ! Use a separate destination buffer (double-buffering) to avoid overlap.
  !   allocate(x_new(walkers%ndim, n_new))
  !   do i = 1, n_old
  !      do j = 1, counts(i)
  !         x_new(:, offsets(i) + j) = walkers%x(:, i)
  !      end do
  !   end do

  !   ! -------- Commit new generation --------
  !   walkers%x(:, 1:n_new) = x_new(:, 1:n_new)
  !   walkers%w(1:n_new)    = 1
  !   walkers%n             = n_new

  !   deallocate(counts, offsets, x_new)
  ! end subroutine walkers_branch

  subroutine compute_potential(walkers, phi)
    type(walkers_t), intent(in) :: walkers
    real(dp), intent(out)       :: phi(walkers%n_max)
    integer                     :: i

    do i = 1, walkers%n
       ! phi(i) = 0.5_dp * sum(walkers%x(:, i)**2)
       phi(i) = -1.0_dp / norm2(walkers%x(:, i))
    end do
  end subroutine compute_potential

end program qmc
