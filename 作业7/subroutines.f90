subroutine read_uv(U, V)
    implicit none
    real(8) :: U(17, 17), V(17, 17)

    open(1, file='u.txt', status='old')
    read(1, *) U
    close(1)
    
    open(2, file='v.txt', status='old')
    read(2, *) V
    close(2)
end subroutine
