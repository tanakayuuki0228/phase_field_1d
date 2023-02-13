program test
    use,intrinsic :: iso_fortran_env
    implicit none
    character(8) :: date
    character(10) :: time
    character(5) :: zone
    integer :: value(8)
    call date_and_time(date,time,zone,value)
    write(*,*) 'date ',date
    write(*,*) 'time ',time
end program test