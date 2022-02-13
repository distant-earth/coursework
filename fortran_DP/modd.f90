module modd

contains

! Правые части уравнений через истинную аномалию:
subroutine func_f(N, t, y, res, RPAR, IPAR)
implicit none
integer N, IPAR
real(8) t, y(N), res(N), a, ecc, omega, f, n_mean, RPAR, mass
real(8), parameter :: G = 32.633120664043254, alpha = 1e-5, n_exp = 0.0d0 !yr
    a = y(1)
    ecc = y(2)
    omega = y(3)
    f = y(4)
    mass = y(5)
    n_mean= sqrt(G * mass / a**3)
    res(5) = - alpha * mass**n_exp
    res(1) = a * ((1 + ecc**2 + 2 * ecc * cos(f)) / (ecc**2 - 1)) * (res(5) / mass)
    res(2) = - (ecc + cos(f)) * (res(5) / mass)
    res(3) = - (sin(f) / ecc) * (res(5) / mass)
    res(4) = - res(3) + n_mean * (1 + ecc * cos(f))**2 / (1 - ecc**2)**(1.5)
end subroutine func_f


! Правые части уравнений через эксцентрическую аномалию:
subroutine func_E(N, t, y, res, RPAR, IPAR)
implicit none
integer N, IPAR
real(8) t, y(N), res(N), a, ecc, omega, E, n_mean, RPAR, mass
real(8), parameter :: G = 32.633120664043254, alpha = 1e-5, n_exp = 0.0d0 !yr
    a = y(1)
    ecc = y(2)
    omega = y(3)
    E = y(4)
    mass = y(5)
    n_mean= sqrt(G * mass / a**3)
    res(5) = - alpha * mass**n_exp
    res(1) = - a * ((1 + ecc * cos(E)) / (1 - ecc * cos(E))) * (res(5) / mass)
    res(2) = ((1 - ecc**2) * cos(E) / (ecc * cos(E) - 1)) * (res(5) / mass)
    res(3) = (sqrt(1 - ecc**2) * sin(E) / (ecc * (ecc * cos(E) - 1))) * (res(5) / mass)
    res(4) = - res(3)/ sqrt(1 - ecc**2) + n_mean / (1 - ecc * cos(E))
end subroutine func_E


subroutine solout(NR, t_old, t, y, N, con, icomp, ND, RPAR, IPAR, IRTRN)
integer IPAR, NR, N, IRTRN, ND
real(8) RPAR, t_old, t, y(N), con, icomp

if (mod(NR, 10000) == 0) write(IPAR, *) NR, t, y

end subroutine solout

end module modd
