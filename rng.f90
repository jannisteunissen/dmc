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
