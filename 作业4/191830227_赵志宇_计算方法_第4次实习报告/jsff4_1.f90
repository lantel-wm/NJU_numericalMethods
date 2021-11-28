program jsff4_1
    ! homework4_1 of Numerical Methods
    ! arthor : zzy

    implicit none
    ! dp : set presion for literal
    integer, parameter :: dp = SELECTED_REAL_KIND(15)
    ! X : x coordinates, Y : y coordinates
    real(8), dimension(10) :: X = [1.02_dp, 0.95_dp, 0.87_dp, 0.77_dp, 0.67_dp, 0.56_dp, 0.44_dp, 0.3_dp, 0.16_dp, 0.01_dp]
    real(8), dimension(10) :: Y = [0.39_dp, 0.32_dp, 0.27_dp, 0.22_dp, 0.18_dp, 0.15_dp, 0.13_dp, 0.12_dp, 0.13_dp, 0.15_dp]
    real(8), dimension(10) :: del_X = [-0.0029_dp, 0.0007_dp, -0.0082_dp, -0.0038_dp, -0.0041_dp, &
                                        0.0026_dp, -0.0001_dp, -0.0058_dp, -0.0005_dp, -0.0034_dp]
    real(8), dimension(10) :: del_Y = [-0.0033_dp, 0.0043_dp, 0.0006_dp, 0.002_dp, 0.0044_dp, &
                                        0.0009_dp, 0.0028_dp, 0.0034_dp, 0.0059_dp, 0.0024_dp]
    ! x_train : character variables in linear regression
    real(8), dimension(10, 5) :: x_train
    ! y_train : target variables in linear regression
    real(8), dimension(10) :: y_train
    ! i : loop variable, n : the number of samples, m : the number of chracters
    ! 10 points and 5 characters(1, x, y, x*y, y^2) are used in linear regression
    integer(4) :: i, n = 10, m = 5
    
    do i = 1, n
        x_train(i, 1) = 1
        x_train(i, 2) = X(i)
        x_train(i, 3) = Y(i)
        x_train(i, 4) = X(i) * Y(i)
        x_train(i, 5) = Y(i) * Y(i)
        y_train(i) = X(i) * X(i)
    end do

    call linear_regression(x_train, y_train, n, m)

    do i = 1, n
        x_train(i, 1) = 1
        x_train(i, 2) = X(i) + del_X(i)
        x_train(i, 3) = Y(i) + del_Y(i)
        x_train(i, 4) = (X(i) + del_X(i)) * (Y(i) + del_Y(i))
        x_train(i, 5) = (Y(i) + del_Y(i)) * (Y(i) + del_Y(i))
        y_train(i) = (X(i) + del_X(i)) * (X(i) + del_X(i))
    end do

    call linear_regression(x_train, y_train, n, m)

end program jsff4_1

subroutine linear_regression(A, y, n, m)
    ! apply linear regression algorithm
	! parameters: A : matrix of character variables, shape is (n, m)
    !             y : vector of target variables
    !             n : the number of (x, y) 
    !             m : the number of characters
	! author: zzy

    implicit none
    integer(4), intent(in) :: n, m
    real(8), dimension(n, m) :: A
    ! B : agumented matrix
    real(8), dimension(m, m + 1) :: B, B0
    real(8), dimension(n) :: y
    ! theta : solution of ATA*b == ATy 
    real(8), dimension(m) :: theta
    ! y_mean : mean value of y, Syy : variance of y, Q : sum of squared error (SSE), R : multiple correlation coefficient
    real(8) :: y_mean, Syy, Q, R
    integer(4) :: i

    ! intialize agumented matrix
    B(1: m, 1: m) = matmul(transpose(A), A)
    B(:, m + 1) = matmul(transpose(A), y)
    B0 = B

    ! solve the equation ATA*b == ATy
    call LU_factoriation(B0, theta, m)

    print *, 'b :', theta

    ! calculate y_mean, Syy, Q, R
    y_mean = 0
    Syy = 0
    Q = 0

    y_mean = sum(y) / dble(n)

    do i = 1, n
        Syy = Syy + (y(i) - y_mean) ** 2
        Q = Q + (y(i) - dot_product(theta, A(i, :))) ** 2
    end do

    R = sqrt((Syy - Q) / Syy)

    ! print *, 'y_mean :', y_mean
    ! print *, 'Syy :', Syy
    ! print *, 'Q :', Q
    ! print *, 'R :', R

end subroutine linear_regression

subroutine LU_factoriation(A, theta, n)
    ! apply LU factoriation, calculate inverse matrix of A(1:n, 1:n)
	! parameters: A : agumented matrix
    !             theta : solution of linear equations
    !             n : the length of theta is n
	! author: zzy

    implicit none
    integer(4), intent(in) :: n
    real(8), intent(in out), dimension(n, n + 1) :: A
    ! LU combines the matrix L and U (PA = LU)
    real(8), dimension(n, n) :: LU, L_inv, U_inv
    ! A(1:n, 1:n) * theta = A(:, n+1) 
    real(8), dimension(n), intent(in out) :: theta
    ! L * zeta = P * A(:, n+1)
    ! U * theta = zeta
    real(8), dimension(n) :: zeta
    ! temp : intermediate varible for vector swap
    real(8), dimension(n + 1) :: temp
    ! cond_inf : conditional number
    real(8) :: cond_inf
    ! i, j, k, r : loop varibles
    integer(4) :: i, j, k, r
    ! p : save the output of maxloc
    integer(4) :: p(1)

    do r = 1, n - 1
        ! find column pivot, and swap the rows
        p = maxloc(abs(A(r: n, r)))
        if (p(1) > r) then
            temp = A(p(1), :)
            A(p(1), :) = A(r, :)
            A(r, :) = temp
        end if

        ! calculate row r of U
        do j = r, n
            LU(r, j) = A(r, j)
            do k = 1, r - 1
                LU(r, j) = LU(r, j) - LU(r, k) * LU(k, j)
            end do
        end do

        ! calculate column r of L
        do i = r + 1, n
            LU(i, r) = A(i, r)
            do k = 1, r - 1
                LU(i, r) = LU(i, r) - LU(i, k) * LU(k, r)
            end do
            LU(i, r) = LU(i, r) / LU(r, r)
        end do
    end do

    ! calculate U(n, n)
    LU(n, n) = A(n, n)
    do k = 1, n - 1
        LU(n, n) = LU(n, n) - LU(n, k) * LU(k, n)
    end do

    ! solve L * zeta = P * A(:, n+1)
    zeta(1) = A(1, n + 1)
    do r = 2, n
        zeta(r) = A(r, n + 1)
        do j = 1, r - 1
            zeta(r) = zeta(r) - LU(r, j) * zeta(j)
        end do
    end do

    ! solve U * theta = zeta
    theta(n) = zeta(n) / LU(n, n)
    do r = n - 1, 1, -1
        theta(r) = zeta(r)
        do j = r + 1, n
            theta(r) = theta(r) - LU(r, j) * theta(j)
        end do
        theta(r) = theta(r) / LU(r, r)
    end do

    ! calc inv(U) and inv(L)
    do i = 1, n
        U_inv(i, i) = 1 / LU(i, i)
        do k = i - 1, 1, -1
            U_inv(k, i) = 0
            do j = k + 1, i
                U_inv(k, i) = U_inv(k, i) - LU(k, j) * U_inv(j, i)
            end do
            U_inv(k, i) = U_inv(k, i) / LU(k, k)
        end do
    end do

    do i = 1, n
        L_inv(i, i) = 1
        do k = i + 1, n
            L_inv(k, i) = 0
            do j = i, k - 1
                L_inv(k, i) = L_inv(k, i) - LU(k, j) * L_inv(j, i)
            end do
        end do
    end do

    cond_inf = maxval(sum(abs(A(1: n, 1: n)), 2)) * maxval(sum(abs(matmul(U_inv, L_inv)), 2))
    print *, "cond_inf :", cond_inf

end subroutine LU_factoriation

subroutine print_matrix(A, m, n)
    ! debug function, print a matrix
	! parameters: A : matrix to be printed
    !             (m, n) : shape of matrix
	! author: zzy

    implicit none
    integer(4) :: m, n, i
    real(8), dimension(m, n) :: A

    do i = 1, m
        print *, A(i, :)
    end do

end subroutine print_matrix