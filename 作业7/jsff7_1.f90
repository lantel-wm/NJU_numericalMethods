subroutine step1()
    implicit none
    real(8) :: U(17, 17), V(17, 17)

    call read_uv(U, V)
    print *, U(1, 3)
    print *, V(3, 1)
end subroutine
