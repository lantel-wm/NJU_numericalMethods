program jsff7
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(8) :: U(0:18, 0:18), V(0:18, 0:18)
    real(8) :: U_p(17, 17), V_p(17, 17)
    real(8) :: D(17, 17)
    real(8) :: phi(0:18, 0:18)

    call read_uv(U, V)
    ! call read_grid(grid_X, grid_Y)
    call calc_div(D, U, V, 0.25_dp)
    print *, D
    call solve_equation(D, phi, 0.25_dp, 1e-7_dp)
    print *, phi
    call calc_uv(phi,U_p, V_p, 0.25_dp)
end program jsff7
