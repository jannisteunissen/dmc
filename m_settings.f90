module m_settings
  use iso_fortran_env, only: real32, real64

  implicit none
  public

  integer, parameter :: fp = real32

  ! Defaults if not set by -D flags
#ifndef NDIM
#define NDIM 3
#endif
#ifndef NPART
#define NPART 1
#endif

  integer, parameter :: n_dim       = NDIM
  integer, parameter :: n_particles = NPART

#ifndef ATOM_Z
#define ATOM_Z 1.0_fp
#endif
  real(fp), parameter :: atom_z = ATOM_Z

end module m_settings
