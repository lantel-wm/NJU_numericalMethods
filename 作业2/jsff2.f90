program jsff2
    ! homework2 of Numerical Methods
    ! arthor : zzy

    implicit none
    ! X : x coordinates, Y : y coordinates
    ! X(1) = X(21), Y(1) = Y(21) inorder to apply periodical boundary conditions
    real(8), dimension(21) :: X = [-0.20, 0.01, 0.16, 0.30, 0.44, 0.56, 0.67, 0.77, 0.87, 0.95,&
                                    0.99, 0.93, 0.85, 0.73, 0.59, 0.42, 0.29, 0.16, 0.05, -0.11, -0.20]
    real(8), dimension(21) :: Y = [0.22, 0.15, 0.13, 0.12, 0.13, 0.15, 0.18, 0.22, 0.27, 0.32,&
                                    0.39, 0.4, 0.41, 0.42, 0.43, 0.42, 0.41, 0.4, 0.36, 0.32, 0.22]
    ! u : the parameter of the curve X = X(u), Y = Y(u)
    real(8), dimension(21) :: u
    ! x_train : character variables in linear regression
    real(8), dimension(20, 5) :: x_train
    ! y_train : target variables in linear regression
    real(8), dimension(20) :: y_train
    ! i : loop variable, n : the number of samples, m : the number of chracters
    ! 20 + 1 points are used in cubic spline
    ! 20 points and 4 characters(x, y, x*y, y^2) are used in linear regression
    integer(4) :: i, n = 20, m = 4 

    ! initialize u
    do i = 1, n + 1
        u(i) = dble(i)
    end do
    print *, 'Mx : '
    call cubic_spline(u, X, n)
    print *, 'My : '
    call cubic_spline(u, Y, n)

    ! initialize x_train and y_train
    do i = 1, n
        x_train(i, 1) = 1
        x_train(i, 2) = X(i)
        x_train(i, 3) = Y(i)
        x_train(i, 4) = X(i) * Y(i)
        x_train(i, 5) = Y(i) * Y(i)
        y_train(i) = X(i) * X(i)
    end do

    call linear_regression(x_train, y_train, n, m)

end program jsff2

subroutine cubic_spline(x, y, n)
    ! apply cubic spline interpolation algorithm
	! parameters: x, y : coordinates of the points to be interpolated
    !             n : the number of the points to be interpolated
	! author: zzy

    implicit none
    integer(4), intent(in) :: n
    integer(4) :: i
    real(8), intent(in), dimension(n + 1) :: x
    real(8), intent(in), dimension(n + 1) :: y
    ! B : augmented matrix, B's shape is (n + 1, n + 2) since there are (n + 1) points 
    real(8), dimension(n + 1, n + 2) :: B
    ! M : second derivative of spline functions in each interval
    real(8), dimension(n + 1) :: M
    ! h : h(i) = x(i) - x(i - 1)
    real(8), dimension(n) :: h

    ! calculate h
    do i = 2, n + 1
        h(i) = x(i) - x(i - 1)
    end do

    ! calculate B according to three-moment method and periodical boundary condition

    ! M(1) == M(n + 1), periodical boundary condition 
    B(1, 1) = 1
    B(1, n + 1) = -1

    do i = 2, n
        ! alpha(i) * M(i - 1) + 2 * M(i) + (1 - alpha(i)) * M(i + 2) == beta(i)
        B(i, i - 1) = h(i) / (h(i) + h(i + 1))
        B(i, i) = 2.0
        B(i, i + 1) = h(i + 1) / (h(i) + h(i + 1))
        B(i, n + 2) = 6 / (h(i) + h(i + 1)) * ((y(i + 1) - y(i)) / h(i + 1) - (y(i) - y(i - 1)) / h(i))
    end do

    ! periodical boundary condition
    B(n + 1, 2) = h(2) / (h(2) + h(n + 1))
    B(n + 1, n) = -h(n) / (h(2) + h(n + 1))
    B(n + 1, n + 1) = 2
    B(n + 1, n + 2) = 6 / (h(2) + h(n + 1)) * ((y(2) - y(1)) / h(2) - (y(n + 1) - y(n)) / h(n + 1))

    call gauss_elimination(B, M, n)

    print *, M

end subroutine cubic_spline

subroutine linear_regression(A, y, n, m)
    ! apply linear regression algorithm
	! parameters: A : matrix of character variables, shape is (n, m + 1)
    !             y : vector of target variables
    !             n : the number of (x, y) 
    !             m : the number of characters
	! author: zzy

    implicit none
    integer(4), intent(in) :: n, m
    real(8), dimension(n, m + 1) :: A
    ! ATA : result of transpose(A) * A
    real(8), dimension(m + 1, m + 1) :: ATA
    ! B : agumented matrix
    real(8), dimension(m + 1, m + 2) :: B
    real(8), dimension(n) :: y
    ! ATy : result of transpose(A) * y
    real(8), dimension(m + 1) :: ATy
    ! theta : solution of ATA*b == ATy 
    real(8), dimension(m + 1) :: theta
    ! y_mean : mean value of y, Syy : variance of y, Q : sum of squared error (SSE), R : multiple correlation coefficient
    real(8) :: y_mean, Syy, Q, R
    integer(4) :: i, j

    ! call print_matrix(A, n, m + 1)
    ! call print_matrix(y, m + 1, 1)

    ATA = matmul(transpose(A), A)
    ATy = matmul(transpose(A), y)
    
    ! call print_matrix(ATA, m + 1, m + 1)
    ! call print_matrix(ATy, m + 1, 1)

    ! intialize agumented matrix
    do i = 1, m + 1
        do j = 1, m + 1
            B(i, j) = ATA(i, j)
        end do
    end do

    do i = 1, m + 1
        B(i, 6) = ATy(i)
    end do

    ! call print_matrix(B, m + 1, m + 2)

    ! solve the equation ATA*b == ATy
    call gauss_elimination(B, theta, m)

    print *, 'b :', theta

    ! calculate y_mean, Syy, Q, R
    y_mean = 0
    Syy = 0
    Q = 0

    do i = 1, n
        y_mean = y_mean + y(i)
    end do
    y_mean = y_mean / dble(n)

    do i = 1, n
        Syy = Syy + (y(i) - y_mean) ** 2
        Q = Q + (y(i) - dot_product(theta, A(i, :))) ** 2
    end do

    R = sqrt((Syy - Q) / Syy)

    print *, 'y_mean :', y_mean
    print *, 'Syy :', Syy
    print *, 'Q :', Q
    print *, 'R :', R

end subroutine linear_regression

subroutine gauss_elimination(B, theta, n)
    ! apply gauss elimination algorithm
	! parameters: B : agumented matrix
    !             theta : solution of linear equations
    !             n : the length of theta is (n + 1)
	! author: zzy

    implicit none
    integer(4), intent(in) :: n
    real(8), intent(in out), dimension(n + 1, n + 2) :: B
    real(8), intent(in out), dimension(n + 1) :: theta
    integer(4) :: i, j, k
    
    ! use elementary transformation to transform B into upper triangular matrix
    do i = 1, n + 1 ! ii : rows
	    do j = i + 1, n + 2 ! j : columns
            B(i, j) = B(i, j) / B(i, i)
        end do
		B(i, i) = 1
		do j = i + 1, n + 1 ! j : rows
			do k = i + 1, n + 2 ! k : columns
				B(j, k) = B(j, k) - B(j, i) * B(i, k) 
            end do
			B(j, i) = 0
        end do
    end do

    ! solve theta by transform B(1:n+1, 1:n+1) into diagonal matrix 
	do i = n + 1, 1, -1
		do j = i + 1, n + 1
            B(i, n + 2) = B(i, n + 2) - theta(j) * B(i, j)
        end do
		theta(i) = B(i, n + 2);
    end do

end subroutine gauss_elimination

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