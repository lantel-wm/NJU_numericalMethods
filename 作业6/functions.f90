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

function composite_simpson(a, b, n, func, gamma_x)
    ! apply simpson composite quadrature
    ! parameters: a, b: integral interval boundray (a, b)
    !             n: number of small intervals
    !             func: integral function func(x)
    !             gamma_x: parameter x in gamma function
    implicit none
    real(8) :: composite_simpson
    real(8), external :: func
    real(8), intent(in) :: a, b, gamma_x
    real(8) :: h
    integer, intent(in) :: n
    integer :: i

    h = (b - a) / dble(n)
    composite_simpson = func(gamma_x, a) + func(gamma_x, b) + 4 * func(gamma_x, a + 0.5 * h)
    do i = 1, n - 1
        composite_simpson = composite_simpson + 4 * func(gamma_x, a + (dble(i) + 0.5) * h) + 2 * func(gamma_x, a + dble(i) * h)
    end do
    composite_simpson = composite_simpson * h / 6
    return
end function composite_simpson

function composite_tixing(a, b, n, func, gamma_x)
    ! apply trapezium composite quadrature
    ! parameters: a, b: integral interval boundray (a, b)
    !             n: number of small intervals
    !             func: integral function func(x)
    !             gamma_x: parameter x in gamma function
    implicit none
    real(8) :: composite_tixing
    real(8), external :: func
    real(8), intent(in) :: a, b, gamma_x
    real(8) :: h
    integer, intent(in) :: n
    integer :: i

    h = (b - a) / dble(n)
    composite_tixing = func(gamma_x, a) + func(gamma_x, b)
    do i = 1, n - 1
        composite_tixing = composite_tixing + 2 * func(gamma_x, a + dble(i) * h)
    end do
    composite_tixing = composite_tixing * h / 2
    return
end function composite_tixing

subroutine search(a, b, n, int_method, func, eps, gamma_x)
    ! search for n which satisfies error less than eps
    ! parameters: n: to be calculated
    !             int_method: integral method, simpson or trapezium composite quadrature
    !             eps： precision
    !             gamma_x: parameter x in gamma function
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    integer, intent(in out) :: n
    real(8), intent(in) :: a, b, eps, gamma_x
    real(8), external :: int_method, func
    real(8) :: res, res_pre, diff

    n = 1
    diff = 1.0
    res_pre = int_method(a, b, n, func, gamma_x)
    do while(abs(diff) > eps)
        n = n + 1
        res = int_method(a, b, n, func, gamma_x)
        diff = res - res_pre
        res_pre = res
    end do

end subroutine search
