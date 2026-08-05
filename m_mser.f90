! Implementation of the batched MSER (Marginal Standard Error Rules) method
!
! Given a time series of length n, find the number of initial observations
! d to discard so that the marginal standard error of the sample mean
! computed from the retained portion series(d+1 : n) is minimised.
!
! The MSER statistic for a candidate discard count d is defined as:
!
!   MSER(d) = Var(series(d+1 : n)) / (n - d)
!
! where Var denotes the sample variance of the retained observations and
! (n - d) is their count.  This equals the squared standard error of the
! retained sample mean.  The optimal discard point d* is the value of d
! that minimises MSER(d).
module m_mser
  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  public :: mser_find_truncation, mser_analysis

contains

  ! Find optimal truncation point
  !
  ! Note: the implementation uses the population-variance form (dividing by
  ! (n-d) rather than (n-d-1)) in the numerator; since the two forms differ
  ! only by a monotone scaling factor, the minimising d is identical.
  !
  ! The "retain at least half" rule is enforced: d is searched over
  ! 0 <= d <= floor(n/2)
  subroutine mser_find_truncation(n, series, optimal_d, min_mser, steady_mean, mser_values)
    !> Input data size
    integer, intent(in)   :: n
    !> Input data
    real(dp), intent(in)  :: series(n)
    !> Discard the first optimal_d elements; use series(optimal_d+1 : n)
    integer,  intent(out) :: optimal_d
    !> The MSER value at the optimum
    real(dp), intent(out) :: min_mser
    !> The mean of the retained batches at the optimum
    real(dp), intent(out) :: steady_mean
    !> Full MSER curve, indexed 0 .. max_discard
    real(dp), allocatable, intent(out), optional :: mser_values(:)

    integer  :: discard, max_discard, m
    real(dp) :: current_mser, x, old_mean, inv_rm
    real(dp) :: w_mean, w_M2
    logical  :: store_mser

    max_discard = n / 2                       ! retain-at-least-half rule
    store_mser  = present(mser_values)

    if (store_mser) allocate(mser_values(0:max_discard))

    if (n < 2) then
       optimal_d = 0
       min_mser  = huge(1.0_dp)
       steady_mean = huge(1.0_dp)
       if (store_mser) mser_values(0) = huge(1.0_dp)
       return
    end if

    ! Use a backward pass of Welford's online algorithm for accuracy
    w_mean = 0.0_dp
    w_M2   = 0.0_dp
    m      = 0

    min_mser    = huge(1.0_dp)
    optimal_d   = 0
    steady_mean = 0.0_dp

    do discard = n - 1, 0, -1
       ! Add series(discard+1) so the accumulator covers series(discard+1 : n)
       x        = series(discard + 1)
       m        = m + 1
       inv_rm   = 1.0_dp / m
       old_mean = w_mean
       w_mean   = w_mean + (x - w_mean) * inv_rm
       w_M2     = w_M2   + (x - old_mean) * (x - w_mean)

       if (discard > max_discard) cycle

       if (m < 2) then
          if (store_mser) mser_values(discard) = huge(1.0_dp)
          cycle
       end if

       ! Population-variance MSER; minimum coincides with sample-variance form
       current_mser = w_M2 * inv_rm * inv_rm
       if (store_mser) mser_values(discard) = current_mser

       if (current_mser < min_mser) then
          steady_mean = w_mean
          min_mser    = current_mser
          optimal_d   = discard
       end if
    end do

  end subroutine mser_find_truncation

  ! MSER analysis with batch means
  subroutine mser_analysis(n, data, batch_size, correct_autocorr, &
       trunc_obs, steady_state_mean, min_mser_val, status)
    !> Size of input data
    integer,  intent(in)  :: n
    !> Input data
    real(dp), intent(in)  :: data(n)
    !> Batch size
    integer,  intent(in)  :: batch_size
    !> Correct for autocorrelated data
    logical, intent(in)   :: correct_autocorr
    !> trunc_batches * batch_size (observations to discard)
    integer,  intent(out) :: trunc_obs
    !> mean of retained raw observations
    real(dp), intent(out) :: steady_state_mean
    !> MSER value at the optimum
    real(dp), intent(out) :: min_mser_val
    !> 0 = OK, 1 = batch_size > n, 2 = all data truncated
    integer,  intent(out) :: status

    integer :: n_batches, i, start_idx, end_idx, n_retained
    integer :: trunc_batches
    real(dp), allocatable :: batch_means(:)
    real(dp) :: w_mean, num, denom, phi, steady_batch_mean

    status = 0
    n_batches = n / batch_size

    if (n_batches < 1) then
       trunc_batches     = 0
       trunc_obs         = 0
       steady_state_mean = 0.0_dp
       min_mser_val      = huge(1.0_dp)
       status            = 1
       return
    end if

    ! Compute batch means
    allocate(batch_means(n_batches))
    do i = 1, n_batches
       start_idx = (i - 1) * batch_size + 1
       end_idx   = i * batch_size
       batch_means(i) = sum(data(start_idx:end_idx)) / real(batch_size, dp)
    end do

    ! Find truncation on the batch-mean series
    call mser_find_truncation(n_batches, batch_means, trunc_batches, &
         min_mser_val, steady_batch_mean)

    trunc_obs = trunc_batches * batch_size

    ! Steady-state mean: single Welford pass over retained raw data
    n_retained = n - trunc_obs
    if (n_retained < 1) then
       steady_state_mean = 0.0_dp
       status            = 2
       return
    end if

    w_mean = 0.0_dp
    do i = trunc_obs + 1, n
       w_mean = w_mean + (data(i) - w_mean) / real(i - trunc_obs, dp)
    end do
    steady_state_mean = w_mean

    if (correct_autocorr) then
       ! Correct for lag-1 autocorrelation of batch means

       ! Check for enough batches to estimate autocorrelation
       if (n_batches - trunc_batches < 3) return

       num = 0.0_dp
       denom = (batch_means(trunc_batches + 1) - steady_batch_mean)**2
       do i = trunc_batches + 2, n_batches
          denom = denom + (batch_means(i) - steady_batch_mean)**2
          num   = num   + (batch_means(i) - steady_batch_mean) * &
               (batch_means(i-1) - steady_batch_mean)
       end do

       if (abs(denom) <= 0) return
       phi = num / denom

       ! Avoid pathological corrections
       phi = max(-0.999_dp, min(0.999_dp, phi))

       ! Corrected MSER: account for lag-1 autocorrelation
       min_mser_val = min_mser_val * (1.0_dp + phi) / (1.0_dp - phi)
    end if

  end subroutine mser_analysis

end module m_mser
