module m_settings
  use iso_fortran_env, only: real32, real64

  implicit none
  public

#ifdef FLOAT_32
  integer, parameter :: fp = real32
#else
  integer, parameter :: fp = real64
#endif

  ! Defaults if not set by -D flags
#ifndef NDIM
#define NDIM 3
#endif
#ifndef NPART
#define NPART 1
#endif

  integer, parameter :: n_dim       = NDIM
  integer, parameter :: n_particles = NPART

#ifndef ATOMZ
#define ATOMZ 1.0_fp
#endif
  real(fp), parameter :: atom_z = ATOMZ

end module m_settings
