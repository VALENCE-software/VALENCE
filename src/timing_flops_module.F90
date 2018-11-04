module timing_flops
  use tools, only: dp
  implicit none
  integer(8) FLOP
  integer(8) count_determinants
  real(dp) kernel_time
  real(dp) initial_time
  real(dp) guess_time
  real(dp) final_time
end module timing_flops
