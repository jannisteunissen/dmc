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

  ! Guesses for orbital parameters
  ! First index: orbital (1 = 1s, 2 = 2s)
  ! Second index: atomic number Z
  real(fp), parameter :: zeta(2, 10) = reshape([ &
       1.0000_fp, 0.5000_fp, & ! H
       1.6875_fp, 0.0575_fp, & ! He
       2.6906_fp, 1.2792_fp, & ! Li
       3.6848_fp, 1.9120_fp, & ! Be
       4.6795_fp, 2.5762_fp, & ! B
       5.6727_fp, 3.2166_fp, & ! C
       6.6651_fp, 3.8340_fp, & ! N
       7.6579_fp, 4.4531_fp, & ! O
       8.6501_fp, 5.0743_fp, & ! F
       9.6421_fp, 5.6982_fp  & ! Ne
       ], [2, 10])

  type(cfg_t)           :: cfg
  type(walkers_t)       :: walkers
  integer               :: i, z_int, i_step, max_steps, steps_log
  real(dp)              :: dt, dt_step, time, end_time
  real(fp)              :: dt_fp, sqrt_dt, trial_energy_fp
  real(dp)              :: kappa
  real(dp)              :: sum_w, n_eff, n_eff_frac
  real(dp)              :: trial_energy, local_energy
  logical               :: mser_done
  real(dp)              :: sum_w_prev, growth_energy, clock_count_rate
  real(dp), allocatable :: growth_energy_array(:)
  real(dp), allocatable :: local_energy_array(:)
  integer               :: batch_size, trunc_obs, status
  real(dp)              :: min_mser, tmp
  logical, parameter    :: correct_autocorr = .true.
  integer(int64)        :: s_sequential(2), t0, t1
  integer               :: dummy_intvec(0)

  real(dp) :: w_target, w_avg
  real(dp) :: tau_avg, gamma
  real(dp) :: corr, corr_max

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

  call cfg_add(cfg, "Z", 1, "Atom Z number")
  call cfg_add(cfg, "orb_up", dummy_intvec, &
       "Spin up orbitals (1=1s, 2=2s, 3=2px, 4=2py, 5=2pz)", &
       dynamic_size=.true.)
  call cfg_add(cfg, "orb_dn", dummy_intvec, &
       "Spin down orbitals (1=1s, 2=2s, 3=2px, 4=2py, 5=2pz)", &
       dynamic_size=.true.)
  call cfg_add(cfg, "z1", -1.0_dp, &
       "Parameter for 1s orbital (< 0 means auto value)")
  call cfg_add(cfg, "z2", -1.0_dp, &
       "Parameter for 2s orbital (< 0 means auto value)")

  call cfg_update_from_arguments(cfg)
  call cfg_check(cfg)

  call cfg_get_size(cfg, "orb_up", n_up)
  call cfg_get(cfg, "orb_up", orb_up(1:n_up))
  call cfg_get_size(cfg, "orb_dn", n_dn)
  call cfg_get(cfg, "orb_dn", orb_dn(1:n_dn))

  call cfg_get(cfg, "Z", z_int)

  if (z_int < 0) then
     ! Set z-value automatically from number of orbitals
     z_int = n_up + n_dn
  else if (n_up == 0 .and. n_dn == 0) then
     ! Set orbitals automatically from Z
     n_up = (z_int+1)/2
     n_dn = z_int/2
     do i = 1, n_up
        orb_up(i) = i
     end do
     do i = 1, n_dn
        orb_dn(i) = i
     end do
  end if
  atom_z = z_int
  !$acc update device(atom_z)

  call cfg_get(cfg, "z1", tmp)
  z1 = real(tmp, fp)
  if (z1 < 0) z1 = zeta(1, z_int)
  call cfg_get(cfg, "z2", tmp)
  z2 = real(tmp, fp)
  if (z2 < 0) z2 = zeta(2, z_int)
  !$acc update device(z1, z2)

  if (n_up + n_dn /= n_particles) then
     print *, "n_up = ", n_up, "n_dn = ", n_dn, "n_particles = ", n_particles
     error stop "n_up + n_dn /= n_particles"
  end if
  !$acc update device(n_up, n_dn, orb_up, orb_dn)

  call cfg_get(cfg, "Z", z_int)
  atom_z = z_int

  call cfg_get(cfg, "end_time", end_time)
  call cfg_get(cfg, "dt", dt)
  call cfg_get(cfg, "initial_energy", trial_energy)
  trial_energy_fp = real(trial_energy, fp)
  call cfg_get(cfg, "kappa", kappa)
  call cfg_get(cfg, "n_eff_frac", n_eff_frac)
  call cfg_get(cfg, "batch_size", batch_size)
  call cfg_get(cfg, "update_interval", update_interval)

  dt_fp     = real(dt, fp)
  sqrt_dt   = real(sqrt(dt), fp)
  dt_step   = dt * update_interval
  max_steps = ceiling(end_time/dt_step)
  steps_log = max(1, ceiling(max_steps*1e-2_dp))
  tau_avg   = 0.1_dp
  gamma     = dt_step / tau_avg ! tau_avg ~ several autocorrelation times
  corr_max  = 10_dp / dt_step        ! max |E_T| shift per unit time

  call walkers_initialize(cfg, walkers)

  ! Print settings
  write(*, '(A)') repeat("=", 50)
  write(*, '(A)') " Diffusion Monte Carlo - Settings"
  write(*, '(A)') repeat("=", 50)
  write(*, '(A, F12.6)')   " Atom Z number       : ", atom_z
  write(*, '(A, I0)')      " Spin up orbitals    : ", n_up
  if (n_up > 0) write(*, '(A, *(I0, 1X))') "   orb_up            : ", orb_up(1:n_up)
  write(*, '(A, I0)')      " Spin down orbitals  : ", n_dn
  if (n_dn > 0) write(*, '(A, *(I0, 1X))') "   orb_dn            : ", orb_dn(1:n_dn)
  write(*, '(A, I0)')      " Total particles     : ", n_particles
  write(*, '(A, F12.6)')   " z1 (1s param)       : ", z1
  write(*, '(A, F12.6)')   " z2 (2s param)       : ", z2
  write(*, '(A, ES12.4)')  " Time step (dt)      : ", dt
  write(*, '(A, ES12.4)')  " End time            : ", end_time
  write(*, '(A, I0)')      " N walkers           : ", walkers%n
  write(*, '(A, F12.6)')   " Initial energy      : ", trial_energy
  write(*, '(A, F12.6)')   " kappa               : ", kappa
  write(*, '(A, F12.6)')   " n_eff_frac          : ", n_eff_frac
  write(*, '(A, I0)')      " batch_size          : ", batch_size
  write(*, '(A, I0)')      " update_interval     : ", update_interval
  write(*, '(A)') repeat("=", 50)

  allocate(growth_energy_array(max_steps))
  allocate(local_energy_array(max_steps))

  growth_energy = trial_energy
  local_energy  = trial_energy
  mser_done     = .false.
  time          = 0.0_dp
  w_target      = sum(walkers%w) ! fixed reference
  sum_w_prev    = w_target
  w_avg         = w_target       ! running average of total weight

  call system_clock(t0, clock_count_rate)

  do i_step = 1, max_steps

     call walkers_update(walkers%n, walkers%x, walkers%phi, walkers%w, &
          walkers%region, walkers%s, dt_fp, sqrt_dt, trial_energy_fp, &
          update_interval)

     call get_stats(walkers%n, walkers%w, walkers%phi, &
          sum_w, n_eff, local_energy_array(i_step))

     growth_energy_array(i_step) = trial_energy - log(sum_w/sum_w_prev)/dt_step
     sum_w_prev = sum_w

     if (mser_done) then
        growth_energy = growth_energy + &
             (growth_energy_array(i_step) - growth_energy) / (i_step - trunc_obs)
        local_energy = local_energy + &
             (local_energy_array(i_step) - local_energy) / (i_step - trunc_obs)
     end if

     if (n_eff < n_eff_frac * walkers%n) then
        call systematic_resample(walkers)
     end if

     if (mod(i_step, steps_log) == 0) then
        if (.not. mser_done) then
           call mser_analysis(i_step, growth_energy_array, batch_size, &
                correct_autocorr, trunc_obs, growth_energy, min_mser, status)
           mser_done = (trunc_obs > 0 .and. trunc_obs < i_step/10)
        end if
        write(*, "(I4,'%',I8,' ',L,I8,F12.6,F12.6)") &
             ceiling((i_step*1e2_dp)/max_steps), i_step, &
             mser_done, trunc_obs, time, growth_energy
     end if

     time = time + dt_step

     ! Exponential moving average of total weight (integral/slow term)
     w_avg = (1.0_dp - gamma)*w_avg + gamma*sum_w

     corr = -log(w_avg/w_target) / (kappa * dt_step)
     corr = sign(min(abs(corr), corr_max), corr)
     trial_energy = growth_energy + corr
     trial_energy_fp = real(trial_energy, fp)
  end do

  call system_clock(t1)

  call mser_analysis(i_step-1, growth_energy_array, batch_size, &
       correct_autocorr, trunc_obs, growth_energy, min_mser, status)
  write(*, "(A,2F12.5)") " Final result:  ", growth_energy, sqrt(min_mser)
  tmp = max_steps * real(walkers%n, dp) * update_interval
  write(*, "(A,E12.4)") " Total updates: ", tmp
  write(*, "(A,E12.4)") " Updates/s:    ", tmp * clock_count_rate / (t1 - t0)

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
    initial_seed = int(r * real(huge(1_int64), dp), int64)

    call random_number(r)
    s_sequential = int(r * real(huge(1_int64), dp), int64)

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
       walkers%region(i) = nodal_region(walkers%x(:, :, i))
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

  subroutine walkers_update(n, x, phi, w, region, s, &
       dt, sqrt_dt, trial_energy, n_steps)
    integer, intent(in)           :: n
    real(fp), intent(inout)       :: x(n_dim, n_particles, n)
    real(fp), intent(inout)       :: phi(n)
    real(fp), intent(inout)       :: w(n)
    integer, intent(inout)        :: region(n)
    integer(int64), intent(inout) :: s(2, n)
    integer, intent(in)           :: n_steps
    real(fp), intent(in)          :: dt, sqrt_dt, trial_energy
    integer                       :: i, it, idim, p, reg, reg_old
    integer(int64)                :: s1, s2
    real(fp)                      :: phil, wl, phi_old, e, r1, r2

    !$acc parallel loop default(present) &
    !$acc private(s1, s2, phil, wl, phi_old, reg, reg_old, r1, r2, e)
    do i = 1, n
       phil = phi(i)
       wl   = w(i)
       reg  = region(i)
       s1   = s(1, i)
       s2   = s(2, i)

       do it = 1, n_steps

          phi_old = phil
          reg_old = reg

          do p = 1, n_particles
             do idim = 1, n_dim - 1, 2 ! full pairs
#ifdef FLOAT_32
                call box_muller_32(s1, s2, r1, r2)
#else
                call box_muller_64(s1, s2, r1, r2)
#endif
                x(idim,   p, i) = x(idim,   p, i) + sqrt_dt*r1
                x(idim+1, p, i) = x(idim+1, p, i) + sqrt_dt*r2
             end do

             if (iand(n_dim, 1) == 1) then ! odd n_dim, do last dimension
#ifdef FLOAT_32
                call box_muller_32(s1, s2, r1, r2)
#else
                call box_muller_64(s1, s2, r1, r2)
#endif
                x(n_dim, p, i) = x(n_dim, p, i) + sqrt_dt*r1
             end if
          end do

          reg = nodal_region(x(:, :, i))
          if (reg /= reg_old) wl = 0

          call compute_potential(x(:, :, i), phil)
          e  = min(-dt*(0.5_fp*(phi_old + phil) - trial_energy), 3.0_fp)
          wl = wl * exp(e)
       end do

       phi(i)     = phil
       w(i)       = wl
       region(i)  = reg
       s(1, i)    = s1
       s(2, i)    = s2
    end do
  end subroutine walkers_update

  subroutine get_stats(n, w, phi, sum_w, n_eff, E_L)
    integer, intent(in)   :: n
    real(fp), intent(in)  :: w(n)
    real(fp), intent(in)  :: phi(n)
    real(dp), intent(out) :: sum_w
    real(dp), intent(out) :: n_eff
    real(dp), intent(out) :: E_L
    real(dp)              :: sum_w2
    integer               :: i

    sum_w  = 0.0_dp
    sum_w2 = 0.0_dp
    E_L = 0.0_dp

    !$omp parallel do reduction(+: sum_w, sum_w2, E_L)
    !$acc parallel loop reduction(+: sum_w, sum_w2, E_L) default(present)
    do i = 1, n
       sum_w  = sum_w  + w(i)
       sum_w2 = sum_w2 + w(i)**2
       E_L = E_L + w(i) * phi(i)
    end do
    !$acc end parallel loop

    n_eff = sum_w**2 / sum_w2
    E_L = E_L/sum_w
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
    mean_w = step

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
    integer               :: n, m, k
    real(fp)              :: d2, dxk
    real(fp), parameter   :: eps2 = 1.0e-12_fp   ! softening

    phi = 0.0_fp

#if defined(POT_ATOM)
    ! Electron-nucleus terms
    do n = 1, n_particles
       d2 = 0.0_fp
       do k = 1, n_dim
          d2 = d2 + x(k, n)**2
       end do
       phi = phi - atom_z / sqrt(d2 + eps2)
    end do

    ! Electron-electron terms
    do n = 1, n_particles
       do m = n+1, n_particles
          d2 = 0.0_fp
          do k = 1, n_dim
             dxk = x(k,n) - x(k,m)
             d2  = d2 + dxk**2
          end do
          phi = phi + 1.0_fp / sqrt(d2 + eps2)
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

  ! sign of det of the Slater matrix A(i,j) = phi_{orb(j)}(elec i)
  pure integer function nodal_region(x) result(r)
    !$acc routine seq
    real(fp), intent(in) :: x(n_dim, n_particles)
    integer              :: s_up, s_dn
    s_up = det_sign(n_up, x(:, 1:n_up),           orb_up(1:n_up), z1, z2)
    s_dn = det_sign(n_dn, x(:, n_up+1:n_up+n_dn), orb_dn(1:n_dn), z1, z2)
    r = s_up * s_dn
  end function nodal_region

  ! sign of det of Slater matrix A(i,j) = phi_{orb(j)}(elec i)
  ! Uses closed-form determinants (n<=5). Only the sign is needed;
  ! a common positive per-row factor exp(-z2*r/2) is divided out.
  pure integer function det_sign(n, xe, orb, z1, z2) result(sg)
    !$acc routine seq
    integer, intent(in)  :: n
    real(fp), intent(in) :: xe(n_dim, n)
    real(fp), intent(in) :: z1 ! 1s exponent
    real(fp), intent(in) :: z2 ! 2s/2p exponent
    integer,  intent(in) :: orb(:)

    real(fp)           :: a(5, 5)
    real(fp)           :: r, x(n_dim), det
    integer            :: i, j

    if (n == 0) then
       sg = 1
       return
    end if

    ! build Slater matrix (per-row positive factor removed)
    do i = 1, n
       x = xe(:, i)
       r = norm2(x)

       do j = 1, n
          select case (orb(j))
          case (1); a(i,j) = exp(-(z1 - 0.5_fp*z2) * r) ! 1s
          case (2); a(i,j) = 1.0_fp - 0.5_fp*z2*r ! 2s
          case (3); a(i,j) = x(1)              ! 2px
          case (4); a(i,j) = x(2)              ! 2py
          case (5); a(i,j) = x(3)              ! 2pz
          end select
       end do
    end do

    ! closed-form determinant by size
    select case (n)
    case (1)
       det = a(1,1)
    case (2)
       det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    case (3)
       det = det3(a(1:3,1:3))
    case (4)
       det = det4(a(1:4,1:4))
    case (5)
       det = det5(a(1:5,1:5))
    case default
       det = 1.0_fp
    end select

    if (det > 0.0_fp) then
       sg = 1
    else
       sg = -1
    end if
  end function det_sign

  include 'rng.f90'

end program

