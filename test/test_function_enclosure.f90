program test_function_enclosure
  use, intrinsic :: iso_fortran_env, only: wp => real64
  
  use test_model_t, only : model

  implicit none 

  abstract interface 
     function likelihood(para, T)
       import :: wp 
       real(wp) :: likelihood
       real(wp), intent(in) :: para(2)
       integer, intent(in), optional :: T 
     end function likelihood
  end interface

  procedure(likelihood), pointer :: ll => null()
  type(model) :: m

  m = model()
  ll => m%lik!([0.2_wp, 0.5_wp])

  !print*,ll([0.2_wp,0.5_wp])
  
end program test_function_enclosure
