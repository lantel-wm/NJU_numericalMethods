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
    real(8) :: grid_X(17, 17), grid_Y(17, 17)
    
    open(1, file='gird.txt', status='old')
    read(1, *) grid_X, grid_Y
    close(1)

end subroutine
