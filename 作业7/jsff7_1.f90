subroutine step1()
    implicit none
    real(8) :: U(17, 17), V(17, 17)

    call read_uv(U, V)
end subroutine
