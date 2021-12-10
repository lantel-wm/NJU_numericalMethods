subroutine step1()
    implicit none
    real(8) :: U(17, 17), V(17, 17)
    real(8) :: grid_X(17, 17), grid_Y(17, 17)

    call read_uv(U, V)
    print *, U
    call read_grid(grid_X, grid_Y)

end subroutine
