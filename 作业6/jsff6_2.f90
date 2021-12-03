subroutine problem2()
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(8), external :: g, composite_simpson, composite_tixing
    real(8), dimension(5) :: gamma_xs = [1.0_dp, 5.0_dp, 10.0_dp, 2.333333_dp, 3.141593_dp]
    real(8) :: a = 0.0_dp, b = 60.0_dp
    integer :: n, i
    
    print *, "Problem 2" 
    ! call search(a, b, n, composite_simpson, g, 1e-6, gamma_x)
    ! print *, "n=", n, "a=", a, "b=", b, "gamma_x=", gamma_x, "gamma=", composite_simpson(a, b, n, g, gamma_x)
    do i = 1, 5
        ! b = gamma_xs(i) * 6.0_dp
        print *, composite_simpson(a, b, n, g, 1e-6_dp, gamma_xs(i)), n, gamma_xs(i)
    end do
end subroutine problem2