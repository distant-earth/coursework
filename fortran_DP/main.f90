program main
use :: modd
implicit none

integer, parameter :: LWORK = 75, LIWORK = 20
real(8), parameter :: pi = 4 * atan(1.0d0), a0 = 0.006d0, ecc0 = 1e-2, omega0 = 0.2d0, f0 = 1.6d0, &
                      E0 = acos((ecc0 + cos(f0)) / (1 + ecc0 * cos(f0))), &
                      RTOL(5) = 0, ATOL(5) = (/1e-6, 1e-6, 1e-6, 1e-6, 1e-6/), &
                      m0 = 1.0d0
real(8) :: y(5), WORK(LWORK) = 0, t0 = 0, time = 70000 !yr
integer :: IDID, IWORK(LIWORK) = 0, IPAR
character(1) choice

IPAR = 1
IWORK(1) = 2147483647

write(*,*) 'f / E'
read (*,*) choice
select case(choice)
	case('f')
		open(unit = IPAR, file = 'f_RESULT')
		y(1) = a0
		y(2) = ecc0
		y(3) = omega0
		y(4) = f0
		y(5) = m0
		call DOP853(5, func_f, t0, y, time, RTOL, ATOL, 1, solout, 1, WORK, LWORK, IWORK, LIWORK, 0.0, IPAR, IDID)
		write(*,*) time, y
		write(*,*) IDID
	case('E')
		open(unit = IPAR, file = 'E_RESULT')
		y(1) = a0
		y(2) = ecc0
		y(3) = omega0
		y(4) = E0
		y(5) = m0
		call DOP853(5, func_E, t0, y, time, RTOL, ATOL, 1, solout, 1, WORK, LWORK, IWORK, LIWORK, 0.0, IPAR, IDID)
		write(*,*) time, y
		write(*,*) IDID
	case default
		write(*,*) 'Неверный ввод!'
endselect

end
