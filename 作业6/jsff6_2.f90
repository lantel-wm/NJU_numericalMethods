subroutine problem2()
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(8), external :: g, composite_simpson, composite_tixing
    real(8) :: gamma_x = 10.0_dp, a = 0.0_dp, b = 100.0_dp
    integer :: n
    
    print *, "Problem 2" 
    call search(a, b, n, composite_simpson, g, 1e-6, gamma_x)
    ! print *, "n=", n, "a=", a, "b=", b, "gamma_x=", gamma_x, "gamma=", composite_simpson(a, b, n, g, gamma_x)
    print *, n, gamma_x, composite_simpson(a, b, n, g, gamma_x)
end subroutine problem2