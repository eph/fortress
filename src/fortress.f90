!*******************************************************************************
!>
! author: Ed Herbst
!
!# Description
!  Framework for estimating Bayesian models in Fortran using MCMC and SMC.  
!
module fortress

  use fortress_bayesian_model_t, only: fortress_lgss_model
  use fortress_prior_t, only: fortress_abstract_prior
  use fortress_random_t, only: fortress_random
  use fortress_smc_t, only: fortress_smc

  !use fortress_utils, only: fortress_write_array

  implicit none


  character(len=144), parameter :: compilation_date = "${time.strftime('%Y-%m%-d')}$"
  

!*******************************************************************************  
end module fortress
!*******************************************************************************
