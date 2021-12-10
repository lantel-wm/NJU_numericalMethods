subroutine read_uv(U, V)
    implicit none
    real(8) :: U(0:18, 0:18), V(0:18, 0:18)

    open(1, file='u.txt', status='old')
    read(1, *) U
    close(1)
    
    open(2, file='v.txt', status='old')
    read(2, *) V
    close(2)

end subroutine

subroutine read_grid(grid_X, grid_Y)
    implicit none
    real(8) :: grid_X(0:18, 0:18), grid_Y(0:18, 0:18)
    
    open(1, file='grid.txt', status='old')
    read(1, *) grid_X, grid_Y
    close(1)

end subroutine read_grid 

subroutine calc_div(D, U, V, grid_X, grid_Y)
    implicit none
    real(8) :: U(0:18, 0:18), V(0:18, 0:18)
    real(8) :: D(0:18, 0:18)
    real(8) :: grid_X(0:18, 0:18), grid_Y(0:18, 0:18)
    integer :: i, j

    do i = 1, 17
        do j = 1, 17
            D(i, j) = (U(i + 1, j) - U(i - 1, j)) / (grid_X(i + 1, j) - grid_X(i - 1, j))&
                    + (V(i, j + 1) - V(i, j - 1)) / (grid_Y(i, j + 1) - grid_Y(i, j - 1))
        end do
    end do

end subroutine calc_div 

subroutine solve_equation(D, phi, h, eps)
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(8), intent(in) :: D(0:18, 0:18), h, eps
    real(8),intent(in out) ::  phi(0:18, 0:18)
    real(8) :: R(0:18, 0:18)
    real(8) :: alpha = 1.6_dp, diff = 1.0_dp
    integer i, j
    
    do i = 0, 18
        do j = 0, 18
            phi(i, j) = 0.0_dp
        end do
    end do
    do while(diff > eps)
        diff = 1e-7_dp
        do i = 1, 17
            do j = 1, 17
                R(i, j) = (phi(i + 1, j) + phi(i, j + 1) + phi(i - 1, j) + phi(i, j - 1) - 4 * phi(i, j)) / (h * h) - D(i, j)
                phi(i, j) = phi(i, j) + 0.25_dp * alpha * R(i, j)
                diff = max(diff, abs(0.25_dp * alpha * R(i, j)))
            end do
        end do
    end do
    
end subroutine solve_equation
