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

function u(x, t)
    ! function used in Gauss-Raguel integral
    implicit none
    real(8) :: u, x, t
    u = t ** (x - 1)
    return
end function u

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
    ! h: h_{2N}, pre_h: h_N, res: S_{2N}, pre_res: S_N, diff: the error of res
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
    !             n: the number of small intervals
    !             func: integral function func(x)
    !             gamma_x: parameter x in gamma function
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(8) :: composite_simpson
    real(8), external :: func
    real(8), intent(in) :: a, b, eps, gamma_x
    ! h: h_{2N}, pre_h: h_N, res: S_{2N}, pre_res: S_N, diff: the error of res
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

function gauss_raguel(n, func, gamma_x)
    ! apply Gauss-Raguel integral
    ! parameters: n: the number of small intervals
    ! func: integral function
    ! gamma_x: parameter x in gamma function
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    ! wi: coefficients, xi: nodes
    real(8) :: gauss_raguel, xi, wi, res
    real(8), external :: func
    real(8), intent(in) :: gamma_x
    integer, intent(in) :: n
    integer :: i
    character(2) :: str

    ! transfer integer to string
    write(str,"(i0)") n
    ! open ./nodes/datan.txt, the nth line of file datan.txt is "$(xi) $(wi)"
    open(1, file='./nodes/data' // trim(adjustl(str)) // '.txt', status='old')
    res = 0.0_dp
    do i = 1, n
        read(1, *) xi, wi
        res = res + wi * func(gamma_x, xi)
    end do
    close(1)

    gauss_raguel = res
    return
end function gauss_raguel