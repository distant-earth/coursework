program main
use :: modd
implicit none

integer, parameter :: LWORK = 75, LIWORK = 20
real(8), parameter :: pi = 4 * atan(1.0d0), a0 = 0.121d0, ecc0 = 0.71, omega0 = 237 * pi / 180, f0 = 0.0d0, &
                      E0 = acos((ecc0 + cos(f0)) / (1 + ecc0 * cos(f0))), &
                      RTOL(5) = 1e-12, ATOL(5) = (/1e-15, 1e-15, 1e-15, 1e-15, 1e-15/), &
                      m0 =	0.6175d0 ! M_sun
                      ! 646.8717053028447d0 ! M_jup
real(8) :: y(5), WORK(LWORK) = 0.D0, t0 = 0, time = 2e+7 ! yr
integer :: IDID, IWORK(LIWORK) = 0, IPAR
character(1) choice

!=================================================
! Параметры, влияющие на работу интегратора:

IWORK(1) = 2147483647 ! максимальное допустимое кол-во шагов
! Поставлено максимальное возможное число, при значении по умолчанию
! интегратор выдает ошибку "Недостаточно шагов".

! WORK(3) = 0.333d0 ! ограничение снизу на отношение нового шага к старому
! WORK(4) = 6.d0 ! ограничение сверху на отношение нового шага к старому
! WORK(7) = 0.d0 ! начальный размер шага
! В методе DOPRI шаг подбирается автоматически. Если результат выглядит 
! неадекватно, желательно сначала поискать ошибку в шапке main.f90.
! Например, проверить RTOL и ATOL.
!=================================================

IPAR = 1 ! индикатор устройства вывода
! Вывод осуществляется в файл f_RESULT или E_RESULT (см. ниже).

write(*,*) 'Время интегрирования (в годах):', time
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
