module test_particles

  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_basic

    use fortress_particles_t, only: fortress_particles

    integer :: nvar
    real(wp) :: mu(1)
    type(fortress_particles) :: particles
    

    nvar = 1
    particles = fortress_particles(npart=5, nvars=nvar)

    mu = particles%mean()

    call assert_equals(0.0_wp, mu(1))

    particles%particles(1,2) = 1.0_wp

    mu = particles%mean()
    call assert_equals(0.2_wp, mu(1))
    particles%weights(2) = 1000.0_wp
    call particles%normalize_weights(mu(1))
    call assert_equals(1000.8_wp, mu(1), 0.000001_wp)


    call particles%systematic_resampling(0.1_wp)
    mu = particles%mean()
    call assert_equals(1.0_wp, mu(1), 0.000001_wp)
    

  end subroutine test_basic

end module test_particles
