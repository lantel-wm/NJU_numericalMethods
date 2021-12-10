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

subroutine solve_equation()
end subroutine solve_equation
