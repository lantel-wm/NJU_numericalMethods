subroutine step1()
    implicit none
    real(8) :: U(0:18, 0:18), V(0:18, 0:18)
    real(8) :: D(0:18, 0:18)
    real(8) :: phi(0:18, 0:18)
    real(8) :: grid_X(0:18, 0:18), grid_Y(0:18, 0:18)

    call read_uv(U, V)
    call read_grid(grid_X, grid_Y)
    call calc_div(D, U, V, grid_X, grid_Y)

end subroutine step1:w

