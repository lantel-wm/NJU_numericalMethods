program jsff4_2
    ! homework4_2 of Numerical Methods
    ! arthor : zzy

    implicit none
    integer, parameter :: dp = SELECTED_REAL_KIND(15)
    real(8), dimension(30, 30) :: H, H_inv
    real(8), dimension(0:60, 0:30) :: C
    integer(4) :: ns(4) = [6, 8, 10, 15]
    integer(4) :: cur, n = 30

    call get_comb(C, n)
    call get_H(H, H_inv, C, n)

    do cur = 1, 4
        call solve(H(1: ns(cur), 1: ns(cur)), ns(cur))
    end do

end program jsff4_2

subroutine get_comb(C, n)
    ! calculate combinatorial numbers
    ! parameters: C(m, n) : ways of choose n items out of m items
    !             n : upper bound of C
    ! author : zzy

    implicit none
    real(8), dimension(0: 2 * n, 0: n) :: C
    integer(4) :: n, i, j

    C(0, 0) = 1
    do i = 1, 2 * n - 1
        C(i, 0) = 1
        do j = 1, min(i, n)
            C(i, j) = C(i - 1, j) + C(i - 1, j - 1)
        end do
    end do

end subroutine get_comb

subroutine get_H(H, H_inv, C, n)
    ! intitalize H, calculate H_inv and cond_inf
    ! parameters: H : Hilbert maxtrix
    !             H_inv : inverse matrix of Hilbert maxtrix
    !             C : combinatorial numbers
    !             n : upper bound of shape of H
    ! author : zzy

    implicit none
    real(8), dimension(n, n) :: H, H_inv
    real(8), dimension(0: 2 * n, 0: n), intent(in) :: C
    ! cond_inf : condition number (infinity)
    real(8) :: cond_inf
    integer(4) :: n, i, j, k

    ! initialize H by its definition
    do i = 1, n
        do j = 1, n
            H(j, i) = 1 / dble(i + j - 1);
        end do
    end do

    ! calculate H_inv and cond_inf
    do k = 1, n
        do i = 1, k
            do j = 1, k
                H_inv(i, j) = (i + j - 1) * C(k + i - 1, k - j) * C(k + j - 1, k - i) * C(i + j - 2, i - 1) ** 2
                if (mod(i + j, 2) == 1) then
                    H_inv(i, j) = -H_inv(i, j)
                end if
            end do
        end do
        cond_inf = maxval(sum(abs(H(1: k, 1: k)), 2)) * maxval(sum(abs(H_inv(1: k, 1: k)), 2))
        print *, "n =", k, "cond_inf =", cond_inf
    end do

end subroutine get_H

subroutine solve(H, n)
    ! solve the given problem in homework4
    ! parameters: H : Hilbert maxtrix
    !             n : shape of H
    ! author : zzy

    implicit none
    integer, parameter :: dp = SELECTED_REAL_KIND(15)
    real(8), dimension(n, n), intent(in) :: H
    ! L : lower triangular matrix, U : upper triangular matrix, D : diagonal matrix
    ! B_J : iteration matrix B of Jacobi iteration
    ! B_GS : iteration matrix B of Gauss-Seidel iteration
    real(8), dimension(n, n) :: L, U, D, B_J, B_GS
    ! x : the solution of equation Hx = Hx*
    real(8), dimension(n) :: x
    ! lam_J : the maximum absolute eigenvalue(aka spectral radius) of B_J
    ! lam_GS : the maximum absolute eigenvalue(aka spectral radius) of B_GS 
    real(8) :: lam_J, lam_GS
    integer(4) :: n, i, j, k

    ! initialize D, L, U
    do i = 1, n
        do j = 1, n
            D(j, i) = 0.0_dp
            L(j, i) = 0.0_dp
            U(j, i) = 0.0_dp
        end do
    end do

    ! H = L + U + D
    do i = 1, n
        D(i, i) = H(i, i)
        do j = 1, i - 1
            L(i, j) = H(i, j)
        end do
        do j = i + 1, n
            U(i, j) = H(i, j)
        end do
    end do

    ! B_J = -inv(D) * (L + U)
    B_J = -H
    do j = 1, n
        B_J(j, j) = 0
    end do

    ! B_GS = -inv(L + D) * U
    ! calc inv(L + D), saved in D
    L = L + D
    do i = 1, n
        D(i, i) = 1 / L(i, i)
        do k = i + 1, n
            D(k, i) = 0
            do j = i, k - 1
                D(k, i) = D(k, i) - L(k, j) * D(j, i)
            end do
            D(k, i) = D(k, i) / L(k, k)
        end do
    end do

    B_GS = -matmul(D, U)

    ! initialize lambda
    lam_J = 1e8_dp
    lam_GS = 1e8_dp
    
    ! calculate the maximum eigenvalue by power method
    call power_method(B_J, n, 1e-2_dp, lam_J)
    call power_method(B_GS, n, 1e-2_dp, lam_GS)

    print *, "n = ", n
    print *, "lam_J :", abs(lam_J)
    print *, "lam_GS :", abs(lam_GS)

    ! solve Hx = Hx* by Gauss-Seidel iteration 
    call gauss_seidel(H, x, n, 1e-2_dp)
    print *, "x =", x

end subroutine solve

subroutine power_method(A, n, eps, lambda)
    ! apply power method to calculate the largest eigenvalue and corresponding eigenvector
    ! parameters: A : the matrix to be calculated
    !             n : shape of A is (n, n)
    !             eps : precision
    !             lambda : eigenvalue

    implicit none
    real(8), dimension(n, n) :: A
    ! v : iteration vector
    real(8), dimension(n, 2) :: v 
    real(8) :: lambda, lam_temp = 0, eps
    integer(4) :: n, i, j
    
    do i = 1, n
        do j = 1, 2
            v(i, j) = i
        end do
    end do

    do while(abs(lambda - lam_temp) > eps)
        lambda = lam_temp
        v(:, 2) = v(:, 2) / maxval(v(:, 2))
        v(:, 1) = v(:, 2)
        v(:, 2) = matmul(A, v(:, 2))
        lam_temp = dot_product(v(:, 2), matmul(A, v(:, 2))) / dot_product(v(:, 2), v(:, 2))
    end do

end subroutine power_method

subroutine gauss_seidel(A, x, n, eps)
    ! apply Gauss-Seidel iteration to solve linear equtions
    ! parameters: A : coefficient matrix
    !             x : the solution
    !             n : shape of A is (n, n)
    !             eps : precision

    implicit none
    integer, parameter :: dp = SELECTED_REAL_KIND(15)
    real(8), dimension(n, n) :: A
    ! x_star : true solution, b : Ax = b, b = x_star * H
    real(8), dimension(n) :: x, x_star, b
    real(8) :: eps
    integer(4) :: n, i, j, cnt = 0

    do i = 1, n
        x(i) = 0.0_dp
        x_star(i) = 1.0_dp
    end do
    b = matmul(x_star, A)

    ! implement Gauss-Seidel iteration
    do while(maxval(abs(x - x_star)) > eps)
        do i = 1, n
            x(i) = b(i)
            do j = 1, i - 1
                x(i) = x(i) - A(i, j) * x(j)
            end do
            do j = i + 1, n
                x(i) = x(i) - A(i, j) * x(j)
            end do
            x(i) = x(i) / A(i, i)
        end do
        cnt = cnt + 1
    end do
    print *, "cnt:", cnt

end subroutine gauss_seidel

subroutine print_matrix(A, m, n)
    ! debug function, print a matrix
	! parameters: A : matrix to be printed
    !             (m, n) : shape of matrix
	! author: zzy

    implicit none
    integer(4) :: m, n, i
    real(8), dimension(m, n) :: A

    do i = 1, m
        print *, A(i, 1: n)
    end do

end subroutine print_matrix