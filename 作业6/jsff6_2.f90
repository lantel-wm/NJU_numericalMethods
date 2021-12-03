subroutine problem2()
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(8), external :: g, u, composite_simpson, composite_tixing, gauss_raguel
    real(8), dimension(5) :: gamma_xs = [1.0_dp, 5.0_dp, 10.0_dp, 2.333333_dp, 3.141593_dp]
    real(8) :: a = 0.0_dp, b = 60.0_dp
    integer :: n, i
    
    print *, "Problem 2" 

    print *, "Composite Simpson:"
    print *, "\tresult\t\t\t     n\t\t  x"
    do i = 1, 5
        print *, composite_simpson(a, b, n, g, 1e-6_dp, gamma_xs(i)), n, gamma_xs(i)
    end do

    n = 5
    print *, "Gauss-Raguel:"
    print *, "\tresult\t\t\t     n\t\t  x"
    do i = 1, 5
        print *, gauss_raguel(n, u, gamma_xs(i)), n, gamma_xs(i)
    end do



end subroutine problem2