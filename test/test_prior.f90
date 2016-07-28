module test_prior
  use fruit
  implicit none

contains

  subroutine test_prior_read_from_file
    use fortress_prior_t, only: prior, new_prior

    type(prior) :: my_prior

    my_prior = new_prior('test/test_prior.txt')
    call assert_equals(1, my_prior%ptype(1))
    call assert_equals(0.5d0, my_prior%pmean(1))
    call assert_equals(0.2d0, my_prior%pstdd(1))
    call assert_equals(0.000001d0, my_prior%plower(1))
    call assert_equals(0.999999d0, my_prior%pupper(1))

  end subroutine test_prior_read_from_file

  subroutine test_prior_frank_er
    use fortress_prior_t, only: prior, new_prior

    type(prior) :: my_prior

    double precision :: parasim(13), logprior

    my_prior = new_prior('test/s5m13t1p2.txt')

    call assert_equals(13, my_prior%npara)

    parasim = [2.0d0, 0.25d0, 1.5d0, 1.0d0, 0.5d0, 0.90d0, 0.6d0, 0.44d0, 3.8d0, 0.3d0, 0.6d0, 0.80d0, 0.5d0]

    logprior = my_prior%logpdf(parasim)

    call assert_equals(-15.888976966d0, logprior, 0.00001d0)

  end subroutine test_prior_frank_er

end module test_prior
