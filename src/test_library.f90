program test_library

  use fortress, only: fortress_smc, fortress_random

  type(fortress_random) :: rn
  !rn = fortress_random()

  type(fortress_smc) :: test_smc
  !print*,'version ', FORTRESS_VERSION

  test_smc = fortress_smc()

end program test_library
