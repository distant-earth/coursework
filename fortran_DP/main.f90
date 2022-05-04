program main
use :: modd
implicit none

integer, parameter :: LWORK = 75, LIWORK = 20
real(8), parameter :: pi = 4 * atan(1.0d0), a0 = 0.121d0, ecc0 = 0.71d0, omega0 = 237 * pi / 180, f0 = 0.0d0, &
                      E0 = acos((ecc0 + cos(f0)) / (1 + ecc0 * cos(f0))), &
                      RTOL(5) = 0.00001, ATOL(5) = (/1e-8, 1e-8, 1e-8, 1e-8, 1e-8/), &
                      m0 = 646.8717053028447d0 ! M_jup
                      !0.6175d0 ! M_sun
real(8) :: y(5), WORK(LWORK) = 0.D0, t0 = 0, time = 1e+8 ! yr
integer :: IDID, IWORK(LIWORK) = 0, IPAR
character(1) choice

IPAR = 1
IWORK(1) = 2147483647
!WORK(7) = 0.1 ! initial step size
!WORK(3) = 1 ! ограничение снизу на отношение нового шага к старому

write(*,*) time
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
