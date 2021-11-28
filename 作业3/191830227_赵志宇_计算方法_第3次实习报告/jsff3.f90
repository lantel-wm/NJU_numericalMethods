program jsff3
    ! homework3 of Numerical Methods
	! arthor : zzy

    implicit none
    ! A0 : record the primary value of A
    ! A : the matrix to be iterated in QR method
    ! P : the rows of P are eigenvectors
    ! Q : a orthogonal matrix
    ! R : a upper triangular matrix
    ! E : identity matrix
    real(8), dimension(5, 5) :: A0, A, P, Q, R, E
    ! v : initial vector of inverse power method
    real(8), dimension(5, 1) :: v
    ! eps : calculation precision
    real(8), parameter :: eps = 1e-8
    ! S : quadratic sum of A(i, j), where 1 <= i <= n, j < i
    real(8) :: S
    integer(4) :: i, j, n = 5

    ! initialize A, A0, E
    A = reshape([11, -3,  -8,  1,   8, &
                 -6,  5,  12,  6, -18, &
                  4, -2,  -3, -2,   8, &
                -10,  4,  12,  3, -14, &
                 -4,  1,   4, -1,  -1] &
                 , shape(A))
    A0 = A
    do i = 1, n
        do j = 1, n
            if (i == j) then
                E(i, j) = 1.
            else
                E(i, j) = 0.
            end if
        end do
    end do

    ! QR method
    do while(sqrt(S(A, n)) > eps)
        call schmidt_orthogonalization(A, Q, R, n)
        A = matmul(R, Q)
    end do

    ! output eigenvalues
    ! print *, 'eigenvalues:'
    ! do i = 1, n
    !     print "(f8.5)", A(i, i)
    ! end do

    ! calculate eigenvectors by inverse power method
    v = reshape([1., 0., 0., 0., 0.], shape(v))
    call inverse_power_method(A0, v, A(1, 1), n, eps)
    P(:, 1:1) = v

    v = reshape([1., 0., 0., 0., -1.], shape(v))
    call inverse_power_method(A0, v, A(2, 2), n, eps)
    P(:, 2:2) = v

    v = reshape([1., 0., 0., 0., 0.], shape(v))
    call inverse_power_method(A0, v, A(3, 3), n, eps)
    P(:, 3:3) = v

    v = reshape([1., 0., 0., 0., 0.], shape(v))
    call inverse_power_method(A0, v, A(4, 4), n, eps)
    P(:, 4:4) = v

    v = reshape([1., 0., 0., 0., -1.], shape(v))
    call inverse_power_method(A0, v, A(5, 5), n, eps)
    P(:, 5:5) = v

    print *, 'eigenvectors:'
    call print_matrix(P, n, n)
    print *, ' '
    ! call print_matrix(matmul(A0, P), n, n)

    call power_method(A0, n, eps)

end program jsff3

function S(A, n)
    ! calculate quadratic sum of A(i, j), where 1 <= i <= n, j < i
    ! parameters : A : input matrix
    !              n : shape of A is (n, n)
    implicit none
    real(8), dimension(5, 5), intent(in out) :: A
    real(8) :: S
    integer :: i, n

    S = 0
    do i = 2, n
        S = S + dot_product(A(i, 1: i - 1), A(i, 1: i - 1))
    end do

end function S

subroutine schmidt_orthogonalization(A, Q, R, n)
    ! apply schmidt orthogonalization to A
    ! parameters: A : the matrix to be orthogonalize
    !             Q : orthogonal matrix
    !             R : upper triangular matrix
    !             n : shape of A, Q, R is (n, n)
    implicit none
    real(8), dimension(5, 5), intent(in out) :: A, Q, R
    ! beta : orthogonal vectors
    real(8), dimension(5, 5) :: beta
    integer :: i, j, n

    do i = 1, n
        beta(:, i) = A(:, i)
        do j = 1, i - 1
            beta(:, i) = beta(:, i) - dot_product(A(:, i), Q(:, j)) * Q(:, j)
        end do
        Q(:, i) = beta(:, i) / sqrt(dot_product(beta(:, i), beta(:, i)))
        R(i, i) = sqrt(dot_product(beta(:, i), beta(:, i)))
        do j = i + 1, n
            R(i, j) = dot_product(A(:, j), Q(:, i))
        end do
    end do

end subroutine

subroutine power_method(A, n, eps)
    ! apply power method to calculate the largest eigenvalue and corresponding eigenvector
    ! parameters: A : the matrix to be calculated
    !             n : shape of A is (n, n)
    !             eps : precision
    implicit none
    real(8), dimension(n, n) :: A
    ! v : iteration vector
    real(8), dimension(n, 2) :: v 
    ! lambda : initial eigenvalue
    real(8) :: lambda = -1e8, lam_temp, eps
    integer(4) :: n
    
    ! (1)
    v = reshape([1., 2., 3., 4., 5., 1., 2., 3., 4., 5.], shape(v))
    
    ! (2)
    ! v = reshape([-1.00000004e+00,  3.33333333e-01,  1.00000000e+00, -4.21368997e-08, -9.99999958e-01, &
                !  -1.00000004e+00,  3.33333333e-01,  1.00000000e+00, -4.21368997e-08, -9.99999958e-01], shape(v))
    
    ! (3)
    ! v = reshape([1., 0., 0., -2., -1., 1., 0., 0., -2., -1.], shape(v))
    
    print *, 'v0 =', v(:, 1) 
    do while(abs(lambda - lam_temp) > eps)
        ! lambda = v(1, 2) / v(1, 1)
        lambda = lam_temp
        ! print *, lambda, ','
        v(:, 2) = v(:, 2) / maxval(v(:, 2))
        v(:, 1) = v(:, 2)
        v(:, 2) = matmul(A, v(:, 2))
        lam_temp = dot_product(v(:, 2), matmul(A, v(:, 2))) / dot_product(v(:, 2), v(:, 2))
    end do
    print *, 'eigenvalue:'
    print *, lambda
    print *, 'eigenvector:'
    print *, v(:, 2)
    ! print *, matmul(A, v(:, 2))
end subroutine power_method

subroutine inverse_power_method(A, v, l0, n, eps)
    ! apply inverse power method to calculate the eigenvector of given eigenvalue l0
    ! parameters: A : the matrix to be calculated
    !             v : the iteration vector
    !             l0 : the given eigenvalue 
    !             n : shape of A is (n, n)
    !             eps : precision
    implicit none
    real(8), dimension(n, n), intent(in) :: A
    real(8), dimension(n, n + 1) :: B
    ! v1, v : the iteration vectors
    real(8), dimension(n, 1) :: v1, v 
    ! lambda : initial eigenvalue
    real(8) :: l0, eps
    integer(4) :: n, i, j
    
    v1 = v

    do i = 1, 5
        B(: , 1: n) = A
        B(: , n + 1: n + 1) = v1
        do j = 1, n
            B(j, j) = B(j, j) - l0
        end do
        ! print *, lambda
        v1 = v
        call gauss_elimination(B, v, n)
        v = v / maxval(v)

    end do
    ! print *, lambda
    ! print *, v
end subroutine inverse_power_method

subroutine gauss_elimination(A, theta, n)
    ! apply gauss elimination algorithm
	! parameters: B : agumented matrix
    !             theta : solution of linear equations
    !             n : the length of theta is (n + 1)
	! author: zzy

    implicit none
    integer(4), intent(in) :: n
    real(8), intent(in out), dimension(n, n + 1) :: A
    real(8), intent(in out), dimension(n) :: theta
    integer(4) :: i, j, k
    
    ! use elementary transformation to transform B into upper triangular matrix
    do i = 1, n ! ii : rows
	    do j = i + 1, n + 1 ! j : columns
            A(i, j) = A(i, j) / A(i, i)
        end do
		A(i, i) = 1
		do j = i + 1, n ! j : rows
			do k = i + 1, n + 1 ! k : columns
				A(j, k) = A(j, k) - A(j, i) * A(i, k) 
            end do
			A(j, i) = 0
        end do
    end do

    ! solve theta by transform B(1:n+1, 1:n+1) into diagonal matrix 
	do i = n, 1, -1
		do j = i + 1, n
            A(i, n + 1) = A(i, n + 1) - theta(j) * A(i, j)
        end do
		theta(i) = A(i, n + 1);
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