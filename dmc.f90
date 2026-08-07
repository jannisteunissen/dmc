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
  integer               :: i_step, max_steps, steps_log
  real(dp)              :: dt, dt_step, time, end_time
  real(fp)              :: dt_fp, sqrt_dt, trial_energy_fp
  real(dp)              :: kappa
  real(dp)              :: sum_w, n_eff, n_eff_frac, trial_energy
  logical               :: mser_done
  real(dp)              :: sum_w_prev, growth_energy, clock_count_rate
  real(dp), allocatable :: growth_energy_array(:)
  integer               :: batch_size, trunc_obs, status
  real(dp)              :: min_mser, tmp
  logical, parameter    :: correct_autocorr = .true.
  integer(int64)        :: s_sequential(2), t0, t1

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

  dt_fp     = real(dt, fp)
  sqrt_dt   = real(sqrt(dt), fp)
  dt_step   = dt * update_interval
  max_steps = ceiling(end_time/dt_step)
  steps_log = ceiling(max_steps*1e-2_dp)

  call walkers_initialize(cfg, walkers)

  allocate(growth_energy_array(max_steps))

  mser_done  = .false.
  time       = 0.0_dp
  sum_w_prev = sum(walkers%w)

  call system_clock(t0, clock_count_rate)

  do i_step = 1, max_steps

     call walkers_update(walkers%n, walkers%x, walkers%phi, walkers%w, &
          walkers%region, walkers%s, dt_fp, sqrt_dt, trial_energy_fp, update_interval)

     call get_stats(walkers%n, walkers%w, sum_w, n_eff)

     growth_energy_array(i_step) = trial_energy - log(sum_w/sum_w_prev)/dt_step
     sum_w_prev = sum_w

     if (mser_done) then
        growth_energy = growth_energy + &
             (growth_energy_array(i_step) - growth_energy) / (i_step - trunc_obs)
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
     trial_energy = growth_energy - log(sum_w/real(walkers%n, dp)) / (kappa*dt_step)
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
    real(fp)                      :: phil, wl, phi_old, r1, r2, e

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
                call box_muller_32(s1, s2, r1, r2)
                x(idim,   p, i) = x(idim,   p, i) + sqrt_dt*r1
                x(idim+1, p, i) = x(idim+1, p, i) + sqrt_dt*r2
             end do

             if (iand(n_dim, 1) == 1) then ! odd n_dim, do last dimension
                call box_muller_32(s1, s2, r1, r2)
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

  subroutine get_stats(n, w, sum_w, n_eff)
    integer, intent(in)   :: n
    real(fp), intent(in)  :: w(n)
    real(dp), intent(out) :: sum_w
    real(dp), intent(out) :: n_eff
    real(dp)              :: sum_w2
    integer               :: i

    sum_w  = 0.0_dp
    sum_w2 = 0.0_dp

    !$omp parallel do reduction(+: sum_w, sum_w2)
    !$acc parallel loop reduction(+: sum_w, sum_w2) default(present)
    do i = 1, n
       sum_w  = sum_w  + w(i)
       sum_w2 = sum_w2 + w(i)**2
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

    phi = 0.0_fp

#if defined(POT_ATOM)
    ! Electron-nucleus terms
    do n = 1, n_particles
       d2 = 0.0_fp
       do k = 1, n_dim
          d2 = d2 + x(k, n)**2
       end do
       phi = phi - atom_z / sqrt(d2)
    end do

    ! Electron-electron terms
    do n = 1, n_particles
       do m = n+1, n_particles
          d2 = 0.0_fp
          do k = 1, n_dim
             dxk = x(k,n) - x(k,m)
             d2  = d2 + dxk**2
          end do
          phi = phi + 1.0_fp / sqrt(d2)
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
    integer :: s_up, s_dn

    ! orbital codes: 1=1s, 2=2s, 3=2px, 4=2py, 5=2pz
#if NPART == 3
    integer, parameter :: orb_up(*) = [1,2],       orb_dn(*) = [1]
#elif NPART == 4
    integer, parameter :: orb_up(*) = [1,2],       orb_dn(*) = [1,2]
#elif NPART == 5
    integer, parameter :: orb_up(*) = [1,2,3],     orb_dn(*) = [1,2]
#elif NPART == 6
    integer, parameter :: orb_up(*) = [1,2,3,4],   orb_dn(*) = [1,2]
#elif NPART == 7
    integer, parameter :: orb_up(*) = [1,2,3,4,5], orb_dn(*) = [1,2]
#elif NPART == 8
    integer, parameter :: orb_up(*) = [1,2,3,4,5], orb_dn(*) = [1,2,3]
#elif NPART == 9
    integer, parameter :: orb_up(*) = [1,2,3,4,5], orb_dn(*) = [1,2,3,4]
#elif NPART == 10
    integer, parameter :: orb_up(*) = [1,2,3,4,5], orb_dn(*) = [1,2,3,4,5]
#endif

#if NPART < 3

#elif NPART <= 10
    integer, parameter :: n_up = size(orb_up), n_dn = size(orb_dn)

    ! Zeta values indexed by atomic number (Z = 3..10)
    real(fp), parameter :: z1_tab(3:10) = &
      [2.6906_fp, 3.6848_fp, 4.6795_fp, 5.6727_fp, 6.6651_fp, 7.6579_fp, 8.6501_fp, 9.6421_fp]
    real(fp), parameter :: z2_tab(3:10) = &
      [1.2792_fp, 1.9120_fp, 2.5762_fp, 3.2166_fp, 3.8340_fp, 4.4531_fp, 5.0743_fp, 5.6982_fp]
    real(fp), parameter :: z1 = z1_tab(ATOMZ)
    real(fp), parameter :: z2 = z2_tab(ATOMZ)

    s_up = det_sign(n_up, x(:, 1:n_up),           orb_up, z1, z2)
    s_dn = det_sign(n_dn, x(:, n_up+1:n_up+n_dn), orb_dn, z1, z2)

    r = s_up * s_dn
#else
    error stop "Not implemented"
#endif
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

    real(fp)           :: a(n, n)
    real(fp)           :: rk, ek, det
    integer            :: i, j

    ! build Slater matrix (per-row positive factor removed)
    do i = 1, n
       rk = norm2(xe(:, i))
       ek = exp(-(z1 - 0.5_fp*z2) * rk)          ! 1s, relative factor
       do j = 1, n
          select case (orb(j))
          case (1); a(i,j) = ek                    ! 1s
          case (2); a(i,j) = 1.0_fp - 0.5_fp*z2*rk ! 2s
          case (3); a(i,j) = xe(1, i)              ! 2px
          case (4); a(i,j) = xe(2, i)              ! 2py
          case (5); a(i,j) = xe(3, i)              ! 2pz
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
       det = 0.0_fp
    end select

    if (det > 0.0_fp) then
       sg = 1
    else
       sg = -1
    end if
  end function det_sign

  !----------------------------------------------------------------
  ! Fixed-size determinant helpers (cofactor expansion)
  !----------------------------------------------------------------
  pure real(fp) function det3(m) result(d)
    !$acc routine seq
    real(fp), intent(in) :: m(3,3)
    d = m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2)) &
      - m(1,2)*(m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
      + m(1,3)*(m(2,1)*m(3,2) - m(2,2)*m(3,1))
  end function det3

  pure real(fp) function minor3(m, rskip, cskip) result(d)
    !$acc routine seq
    real(fp), intent(in) :: m(4, 4)
    integer,  intent(in) :: rskip, cskip
    real(fp)             :: s(3,3)
    integer              :: ri, ci, i, j
    ri = 0
    do i = 1, size(m,1)
       if (i == rskip) cycle
       ri = ri + 1
       ci = 0
       do j = 1, size(m,2)
          if (j == cskip) cycle
          ci = ci + 1
          s(ri,ci) = m(i,j)
       end do
    end do
    d = det3(s)
  end function minor3

  pure real(fp) function det4(m) result(d)
    !$acc routine seq
    real(fp), intent(in) :: m(4,4)
    ! expand along first row
    d =  m(1,1)*minor3(m,1,1) - m(1,2)*minor3(m,1,2) &
       + m(1,3)*minor3(m,1,3) - m(1,4)*minor3(m,1,4)
  end function det4

  pure real(fp) function minor4(m, cskip) result(d)
    !$acc routine seq
    real(fp), intent(in) :: m(5,5)
    integer,  intent(in) :: cskip
    real(fp) :: s(4,4)
    integer  :: ci, i, j
    do i = 1, 4                 ! drop row 1
       ci = 0
       do j = 1, 5
          if (j == cskip) cycle
          ci = ci + 1
          s(i,ci) = m(i+1,j)
       end do
    end do
    d = det4(s)
  end function minor4

  pure real(fp) function det5(m) result(d)
    !$acc routine seq
    real(fp), intent(in) :: m(5,5)
    ! expand along first row
    d =  m(1,1)*minor4(m,1) - m(1,2)*minor4(m,2) &
       + m(1,3)*minor4(m,3) - m(1,4)*minor4(m,4) &
       + m(1,5)*minor4(m,5)
  end function det5

  include 'rng.f90'

end program dmc
