program jsff5
    ! homework5 of Numerical Methods
    ! arthor : zzy

    implicit none
    integer, parameter :: dp = SELECTED_REAL_KIND(15)

    call newton_iteration(1928332946.0_dp, 1e-6_dp)
    call bisection(1928332946.0_dp, 1e-6_dp)
    
end program jsff5

subroutine newton_iteration(a_double, eps)
    ! apply newton iteration
    ! parameters: a_double : double_precision of a
    !             eps : precision

    implicit none
    ! x_double, x_double_temp : variables used in iteration
    real(8) :: x_double, x_double_temp, a_double, eps
    ! x_float, x_float_temp : variables used in iteration
    real(4) :: x_float, x_float_temp, a_float
    ! cnt_double, cnt_float : record the iteration times
    integer(4) :: cnt_double = 0, cnt_float = 0

    ! initialize variables
    a_float = sngl(a_double)
    x_double = a_double
    x_float = a_float

    ! iteration of double precision variable
    x_double_temp = x_double - 0.5 * (x_double - a_double / x_double)
    do while(abs(x_double - x_double_temp) > eps)
        x_double = x_double_temp
        x_double_temp = x_double - 0.5 * (x_double - a_double / x_double)
        cnt_double = cnt_double + 1
    end do

    ! iteration of single precision variable
    x_float_temp = x_float - 0.5 * (x_float - a_float / x_float)
    do while(abs(x_float - x_float_temp) > eps)
        x_float = x_float_temp
        x_float_temp = x_float - 0.5 * (x_float - a_float / x_float)
        cnt_float = cnt_float + 1
    end do

    ! output the result
    print *, "Newton Iteration:"
    print *, "double:"
    print *, "a = ", a_double
    print *, "iteration times =", cnt_double
    print *, "x =", x_double
    print *, "sqrt(a) =", sqrt(a_double)
    print *, "error =", abs(x_double * x_double - a_double)
    print *, ""
    
    
    print *, "float:"
    print *, "a = ", a_float
    print *, "iteration times =", cnt_float
    print *, "x =", x_float
    print *, "sqrt(a) =", sqrt(a_float)
    print *, "error =", abs(x_float * x_float - a_float)
    print *, ""

end subroutine newton_iteration

subroutine bisection(a_double, eps)
    ! apply bisection method
    ! parameters: a_double : double_precision of a
    !             eps : precision

    implicit none
    ! l_double, r_double : lower and upper bound of current interval
    ! x_double : the solution
    real(8) :: l_double, r_double, x_double, a_double, eps
    ! l_float, r_float : lower and upper bound of current interval
    ! x_float : the solution
    real(4) :: l_float, r_float, x_float, a_float
    ! cnt_double, cnt_float : record the iteration times
    integer(4) :: cnt_double = 0, cnt_float = 0

    ! initialize variables
    a_float = sngl(a_double)

    l_double = 0
    r_double = a_double

    ! bisection
    do while(abs(r_double - l_double) > eps)
        x_double = l_double + 0.5 * (r_double - l_double)
        if (x_double * x_double > a_double) then
            r_double = x_double
        else
            l_double = x_double
        end if
        cnt_double = cnt_double + 1
    end do

    ! initialize variables
    l_float = 0
    r_float = a_float
    x_float = l_float + 0.5 * (r_float - l_float)

    ! bisection
    do while(abs(r_float - l_float) > eps)
        x_float = l_float + 0.5 * (r_float - l_float)
        if (x_float * x_float - a_float < 0) then
            if(abs(l_float - x_float) < eps) then
                exit
            end if
            l_float = x_float
        else
            if(abs(r_float - x_float) < eps) then
                exit
            end if
            r_float = x_float
        end if
        cnt_float = cnt_float + 1
    end do

    print *, "Bisection Method:"
    print *, "double:"
    print *, "a = ", a_double
    print *, "iteration times =", cnt_double
    print *, "x =", x_double
    print *, "sqrt(a) =", sqrt(a_double)
    print *, "error =", abs(x_double * x_double - a_double)
    print *, "l, r =", l_double, r_double
    print *, ""
    
    print *, "float:"
    print *, "a = ", a_float
    print *, "iteration times =", cnt_float
    print *, "x =", x_float
    print *, "sqrt(a) =", sqrt(a_float)
    print *, "error =", abs(x_float * x_float - a_float)
    print *, "l, r =", l_float, r_float

end subroutine bisection