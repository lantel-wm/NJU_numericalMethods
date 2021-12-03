function f(gamma_x, x)
    ! return f(x)
    implicit none
    real(8) :: f, gamma_x, x
    f = 1 / sqrt(1 + x ** 3)
    return
end function f

function g(x, t)
    ! function used to calculate gamma function
    implicit none
    real(8) :: g, x, t
    g = t ** (x - 1) * exp(-t)
    return
end function g

function composite_tixing(a, b, n, func, eps, gamma_x)
    ! apply variable step trapezium composite quadrature
    ! parameters: a, b: integral interval boundray (a, b)
    !             n: number of small intervals
    !             func: integral function func(x)
    !             eps: precision
    !             gamma_x: parameter x in gamma function
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(8) :: composite_tixing
    real(8), external :: func
    real(8), intent(in) :: a, b, eps, gamma_x
    real(8) :: h, pre_h, res, pre_res, diff
    integer, intent(in out) :: n
    integer :: i

    n = 1
    pre_h = b - a;
    pre_res = 0.5 * (b - a) * (func(gamma_x, a) + func(gamma_x, b))
    diff = 1.0

    do while(abs(diff) > eps)
        res = 0_dp
        h = 0.5_dp * pre_h

        do i = 1, n
            res = res + func(gamma_x, a + (2.0_dp * dble(i) - 1) * h)
        end do

        res = 0.5_dp * pre_res + h * res
        diff = (res - pre_res) / 3.0_dp
        pre_h = h
        pre_res = res
        n = n * 2
    end do

    composite_tixing = res
    return
end function composite_tixing

function composite_simpson(a, b, n, func, eps, gamma_x)
    ! apply variable step simpson composite quadrature
    ! parameters: a, b: integral interval boundray (a, b)
    !             n: number of small intervals
    !             func: integral function func(x)
    !             gamma_x: parameter x in gamma function
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(8) :: composite_simpson
    real(8), external :: func
    real(8), intent(in) :: a, b, eps, gamma_x
    real(8) :: h, pre_h, res, pre_res, diff
    integer, intent(in out) :: n
    integer :: i

    n = 1
    pre_h = b - a;
    pre_res = (b - a) * (func(gamma_x, a) + 4.0_dp * func(gamma_x, 0.5_dp * (a + b)) + func(gamma_x, b)) / 6.0_dp
    diff = 1.0

    do while(abs(diff) > eps)
        h = 0.5_dp * pre_h
        res = 0.0_dp
        ! res = 2.0_dp * func(gamma_x, a + (2.0_dp * dble(n) - 0.5_dp) * h)
        do i = 1, 2 * n
            res = res + 2.0_dp * func(gamma_x, a + (dble(i) - 0.5_dp) * h)
        end do

        do i = 1, n
            res = res - func(gamma_x, a + (2.0_dp * dble(i) - 1.0_dp) * h)
        end do

        res = 0.5_dp * pre_res + h * res / 3.0_dp
        diff = (res - pre_res) / 15.0_dp
        pre_h = h
        pre_res = res
        n = n * 2
    end do

    composite_simpson = res
    return
end function composite_simpson

! subroutine search(a, b, n, int_method, func, eps, gamma_x)
!     ! search for n which satisfies error less than eps
!     ! parameters: n: to be calculated
!     !             int_method: integral method, simpson or trapezium composite quadrature
!     !             eps： precision
!     !             gamma_x: parameter x in gamma function
!     implicit none
!     integer, parameter :: dp = selected_real_kind(15)
!     integer, intent(in out) :: n
!     real(8), intent(in) :: a, b, eps, gamma_x
!     real(8), external :: int_method, func
!     real(8) :: res, res_pre, diff

!     n = 1
!     diff = 1.0
!     res_pre = int_method(a, b, n, func, gamma_x)
!     do while(abs(diff) > eps)
!         n = n + 1
!         res = int_method(a, b, n, func, gamma_x)
!         diff = res - res_pre
!         res_pre = res
!     end do

! end subroutine search
