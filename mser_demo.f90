program mser_demo
  use m_mser
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)

  integer, parameter  :: n_total = 10000
  integer  :: batch_size
  real(dp) :: data(n_total)
  integer  :: trunc_batches, trunc_obs
  real(dp) :: steady_state_mean, min_mser_val
  integer  :: i, status
  real(dp) :: r1, r2, u1, u2, noise
  logical  :: correct_autocorr

  ! --- Generate synthetic data with initial transient ---
  ! Model: Y(t) = 10*exp(-t/500) + 5.0 + noise
  ! True steady-state mean = 5.0
  ! The exponential decay creates a warm-up bias

  correct_autocorr = .false.

  ! Simple Box-Muller for Gaussian noise
  do i = 1, n_total
    call random_number(r1)
    call random_number(r2)
    ! Ensure r1 > 0 for log
    r1 = max(r1, 1.0e-15_dp)
    u1 = sqrt(-2.0_dp * log(r1)) * cos(2.0_dp * 3.14159265358979_dp * r2)
    noise = u1 * 1.0_dp  ! std dev = 1.0

    data(i) = 10.0_dp * exp(-real(i, dp) / 500.0_dp) + 5.0_dp + noise
  end do

  ! --- Run MSER analysis with different batch sizes ---
  print *, "=============================================="
  print *, " MSER Analysis with Controllable Batch Means"
  print *, "=============================================="
  print *, ""
  print *, "Synthetic data: Y(t) = 10*exp(-t/500) + 5.0 + N(0,1)"
  print *, "True steady-state mean = 5.0"
  print *, "Total observations:", n_total
  print *, ""

  ! Test with different batch sizes
  do batch_size = 1, 100, 10
    if (batch_size == 1) then
      ! batch_size = 1 is equivalent to MSER without batching
    end if

    call mser_analysis(n_total, data, batch_size, correct_autocorr, &
         trunc_batches, trunc_obs, steady_state_mean, min_mser_val, status)

    print '(A, I5, A, I5, A, I6, A, F10.5, A, ES12.5)', &
      " Batch size=", batch_size, &
      "  Trunc(batches)=", trunc_batches, &
      "  Trunc(obs)=", trunc_obs, &
      "  SS Mean=", steady_state_mean, &
      "  MSER=", min_mser_val
  end do

  print *, ""
  print *, "----------------------------------------------"
  print *, " Detailed run with batch_size = 10"
  print *, "----------------------------------------------"

  batch_size = 10
  call mser_analysis(n_total, data, batch_size, correct_autocorr, &
       trunc_batches, trunc_obs, steady_state_mean, min_mser_val, status)

  print '(A, I6)',    "  Optimal truncation (observations): ", trunc_obs
  print '(A, I6)',    "  Optimal truncation (batches):      ", trunc_batches
  print '(A, F12.6)', "  Estimated steady-state mean:       ", steady_state_mean
  print '(A, ES14.6)', "  Minimum MSER statistic:            ", min_mser_val
  print '(A, I6)',    "  Observations used for estimation:  ", n_total - trunc_obs
  print *, ""

  ! Compute naive mean (no truncation) for comparison
  print '(A, F12.6)', "  Naive mean (no truncation):        ", &
    sum(data) / real(n_total, dp)
  print *, ""

    ! --- Test mser_auto_batch with autocorrelated data ---
  print *, "=============================================="
  print *, " Auto-Batch MSER with Autocorrelated Data"
  print *, "=============================================="
  print *, ""

  block
    real(dp) :: ac_data(n_total)
    real(dp) :: ac_ssm, ac_stderr
    integer  :: ac_bs, ac_trunc, ac_status
    real(dp) :: ar_coeff, transient, ar_state

    correct_autocorr = .true.

    ! Generate AR(1) process with initial transient:
    !   X(t) = phi * X(t-1) + eps(t)  +  transient(t)
    !
    ! where phi = 0.8 (strong positive autocorrelation),
    !       eps ~ N(0, 1),
    !       transient(t) = 20*exp(-t/300)
    !
    ! True steady-state mean = 0.0
    ! Marginal variance = 1 / (1 - phi^2) ≈ 5.26

    ar_coeff = 0.8_dp
    ar_state = 50.0_dp   ! start far from equilibrium

    do i = 1, n_total
       call random_number(r1)
       call random_number(r2)
       r1 = max(r1, 1.0e-15_dp)
       u1 = sqrt(-2.0_dp * log(r1)) * cos(2.0_dp * 3.14159265358979_dp * r2)
       noise = u1 * 1.0_dp

       transient = 20.0_dp * exp(-real(i, dp) / 300.0_dp)
       ar_state  = ar_coeff * ar_state + noise + transient
       ac_data(i) = ar_state
    end do

    print *, "AR(1) process: X(t) = 0.8*X(t-1) + N(0,1) + 20*exp(-t/300)"
    print *, "Initial X(0) = 50.0 (far from equilibrium)"
    print *, "True steady-state mean = 0.0"
    print *, "AR(1) marginal std dev ≈ 2.294"
    print *, ""

    ! First show fixed batch sizes for comparison
    print *, "--- Fixed batch sizes for comparison ---"
    do batch_size = 1, 32
       call mser_analysis(n_total, ac_data, batch_size, correct_autocorr, &
            trunc_batches, trunc_obs, steady_state_mean, min_mser_val, status)
       print '(A, I5, A, I6, A, F10.5, A, ES12.5)', &
         "  Batch=", batch_size, &
         "  Trunc(obs)=", trunc_obs, &
         "  SS Mean=", steady_state_mean, &
         "  MSER=", min_mser_val
    end do
    print *, ""

  end block

end program mser_demo
