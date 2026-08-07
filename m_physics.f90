module m_physics
  use m_settings

  implicit none
  public

contains

  pure subroutine compute_potential(x, phi)
    !$acc routine seq
    real(fp), intent(in)  :: x(n_dim, n_particles)
    real(fp), intent(out) :: phi
    integer               :: n, m

    phi = 0.0_fp

#if defined(POT_ATOM)
    ! Electron-nucleus terms
    do n = 1, n_particles
       phi = phi - atom_z / norm2(x(:, n))
    end do

    ! Electron-electron terms
    do n = 1, n_particles
       do m = n+1, n_particles
          phi = phi + 1.0_fp / norm2(x(:, n) - x(:, m))
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

  pure integer function spin_region(x) result(r)
    !$acc routine seq
    real(fp), intent(in) :: x(n_dim, n_particles)
    real(fp)             :: s

#if defined(REGION_HO2D)
    s = norm2(x(:, 3)) - 2.0_fp/3.0_fp
#else if defined(REGION_NONE)
    s = 1.0_fp
#endif
    r = merge(1, -1, s >= 0.0_fp)
  end function spin_region

end module m_physics
