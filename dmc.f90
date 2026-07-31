! Simple diffusion Monte Carlo method, inspired by the demo codes from Alfonso
! Annarelli, alfonso.annarelli@gmail.com
!
! Author: Jannis Teunissen
program qmc
  use iso_fortran_env, only: error_unit, real32, real64, int32, int64
  use m_config
  use m_random
  use m_mser

  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32

  type walkers_t
     integer               :: n
     integer               :: ndim
     real(dp)              :: n_eff_frac = 0.8_dp
     real(sp), allocatable :: w(:)
     real(sp), allocatable :: x(:, :)
     real(sp), allocatable :: phi(:)
     ! For updating population
     real(sp), allocatable :: x_new(:, :)
     real(sp), allocatable :: phi_new(:)
     real(dp), allocatable :: cdf(:)
     ! Random seed state
     integer(int64), allocatable :: s(:, :)
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
    integer(int64)                 :: initial_seed(2)

    call cfg_get(cfg, "num_walkers", walkers%n)
    call cfg_get(cfg, "ndim", walkers%ndim)

    allocate(walkers%w(walkers%n))
    allocate(walkers%x(walkers%ndim, walkers%n))
    allocate(walkers%phi(walkers%n))

    allocate(walkers%x_new(walkers%ndim, walkers%n))
    allocate(walkers%phi_new(walkers%n))
    allocate(walkers%cdf(walkers%n))

    allocate(walkers%s(2, walkers%n))

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

    ! Set initial seeds
    initial_seed = [1234_int64, -89234_int64]

    walkers%s(:, 1) = initial_seed

    do i = 2, walkers%n
       walkers%s(:, i) = walkers%s(:, i-1)
       call jump(walkers%s(1, i), walkers%s(2, i))
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
    real(dp)                       :: phi_avg, sqrt_dt, exp_arg
    real(dp)                       :: sum_w2, sum_we
    real(real32)                   :: rr(2*((walkers%ndim+1)/2))

    sqrt_dt = sqrt(dt)
    sum_w  = 0.0_dp
    sum_we = 0.0_dp
    sum_w2 = 0.0_dp

    !$omp parallel do private(phi_avg, idim, exp_arg), reduction(+: sum_w, sum_we, sum_w2)
    do i = 1, walkers%n
       phi_avg = walkers%phi(i)

       do idim = 1, (walkers%ndim+1)/2
          call box_muller_32(walkers%s(1, i), walkers%s(2, i), &
               rr(2*idim-1), rr(2*idim))
       end do

       do idim = 1, walkers%ndim
          walkers%x(idim, i) = walkers%x(idim, i) + sqrt_dt * rr(idim)
       end do

       call compute_potential(walkers, i)
       phi_avg = 0.5_dp * (phi_avg + walkers%phi(i))

       exp_arg = -dt * (phi_avg - trial_energy)
       exp_arg = min(exp_arg, 2.0_dp)
       walkers%w(i) = walkers%w(i) * exp(exp_arg)

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
    walkers%phi(i) = -1 / norm2(walkers%x(:, i))
  end subroutine compute_potential

  pure subroutine next(s1, s2, res)
    !$acc routine seq
    integer(int64), intent(inout) :: s1, s2
    integer(int64), intent(out)   :: res
    integer(int64)                :: t1, t2

    t1  = s1
    t2  = s2
    res = t1 + t2
    t2  = ieor(t1, t2)
    s1  = ieor(ieor(ishftc(t1, 24), t2), ishft(t2, 16))
    s2  = ishftc(t2, 37)
  end subroutine next

  ! This is the jump function for the generator. It is equivalent
  ! to 2^64 calls to next(); it can be used to generate 2^64
  ! non-overlapping subsequences for parallel computations.
  pure subroutine jump(s1, s2)
    !$acc routine seq
    integer(int64), intent(inout) :: s1, s2
    integer                       :: i, b
    integer(int64)                :: t1, t2, dummy
    integer(int64), parameter     :: jmp_c(2) = &
         [-2337365368286915419_int64, 1659688472399708668_int64]

    t1 = 0
    t2 = 0
    do i = 1, 2
       do b = 0, 63
          if (iand(jmp_c(i), ishft(1_int64, b)) /= 0) then
             t1 = ieor(t1, s1)
             t2 = ieor(t2, s2)
          end if
          call next(s1, s2, dummy)
       end do
    end do

    s1 = t1
    s2 = t2
  end subroutine jump

  ! A [0, 1) random number
  pure function uni01_64(x) result(u)
    !$acc routine seq
    integer(int64), intent(in) :: x
    integer(int64)             :: y
    real(real64)               :: u
    y = ior(ishft(1023_int64, 52), ishft(x, -12))
    u = transfer(y, u)
  end function uni01_64

  ! A (0, 1] random number
  pure function uni01o_64(x) result(u)
    !$acc routine seq
    integer(int64), intent(in) :: x
    integer(int64)             :: y
    real(real64)               :: u
    y = ior(ishft(1023_int64, 52), ishft(x, -12))
    u = 2.0_real64 - transfer(y, u)
  end function uni01o_64

  ! A [0, 1) single-precision random number from 32 random bits
  pure function uni01_32(x) result(u)
    !$acc routine seq
    integer(int32), intent(in) :: x
    integer(int32)             :: y
    real(real32)               :: u
    ! IEEE-754 single: exponent bias 127 -> 127<<23 = 0x3F800000
    ! use top 23 bits of x as mantissa (shift right by 9)
    y = ior(ishft(127_int32, 23), ishft(x, -9))
    u = transfer(y, u) - 1
  end function uni01_32

  ! A (0, 1] single-precision random number from 32 random bits
  pure function uni01o_32(x) result(u)
    !$acc routine seq
    integer(int32), intent(in) :: x
    integer(int32)             :: y
    real(real32)               :: u
    y = ior(ishft(127_int32, 23), ishft(x, -9))
    u = 2.0_real32 - transfer(y, u)
  end function uni01o_32

  pure subroutine box_muller_64(s1, s2, z1, z2)
    !$acc routine seq
    integer(int64), intent(inout) :: s1, s2
    real(real64),   intent(out)   :: z1, z2
    integer(int64)                :: x
    real(real64)                  :: u1, u2, r, theta
    real(real64), parameter       :: two_pi = 8 * atan(1.0_real64)

    call next(s1, s2, x)
    u1 = uni01o_64(x) ! (0, 1]  -> safe for log

    call next(s1, s2, x)
    u2 = uni01o_64(x)

    r     = sqrt(-2 * log(u1))
    theta = two_pi * u2
    z1    = r * cos(theta)
    z2    = r * sin(theta)
  end subroutine box_muller_64

  pure subroutine box_muller_32(s1, s2, z1, z2)
    !$acc routine seq
    integer(int64), intent(inout) :: s1, s2
    real(real32),   intent(out)   :: z1, z2
    integer(int64)                :: x
    integer(int32)                :: xhi, xlo
    real(real32)                  :: u1, u2, r, theta
    real(real32), parameter       :: two_pi = 8 * atan(1.0_real32)

    call next(s1, s2, x)
    xhi = int(ishft(x, -32), int32)   ! top 32 bits
    xlo = int(x, int32)               ! bottom 32 bits

    u1 = uni01o_32(xhi) ! (0, 1]  -> safe for log
    u2 = uni01o_32(xlo)

    r     = sqrt(-2.0_real32 * log(u1))
    theta = two_pi * u2
    z1    = r * cos(theta)
    z2    = r * sin(theta)
  end subroutine box_muller_32

end program qmc
