subroutine problem1()
    print *, "hello from zzyserver01"    
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(8), external :: f, composite_simpson, composite_tixing
    integer :: n_simpson, n_tixing

    print *, "Problem 1"
    ! call search(0.0_dp, 1.0_dp, n_simpson, composite_simpson, f, 1e-6_dp, 0.0_dp)
    ! call search(0.0_dp, 1.0_dp, n_tixing, composite_tixing, f, 1e-6_dp, 0.0_dp)
    ! print *, n_simpson, composite_simpson(0.0_dp, 1.0_dp, n_simpson, f, 0.0_dp)
    print *, composite_simpson(0.0_dp, 1.0_dp, n_simpson, f, 1e-6_dp, 0.0_dp), n_simpson 
    print *, composite_tixing(0.0_dp, 1.0_dp, n_tixing, f, 1e-6_dp, 0.0_dp), n_tixing 

end subroutine problem1
