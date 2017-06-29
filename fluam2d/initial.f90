program initial
  integer, parameter :: n=32000
  integer, parameter :: ng=n/3
  integer :: p
  real :: x, y, Lx, Ly
  
  Lx=453.74818583181210299
  Ly=Lx
  
  call RANDOM_SEED()
  
  write(*,*) n
  do p=1,n
    call random_number(x)
    call random_number(y)
    x=x*Lx
    if(p<ng) then ! Green
       y=y*Ly/3+Ly/3
    else if((p>=ng).and.(p<n-ng)) then ! Red bottom band
      y=y*Ly/3    
    else
      y=y*Ly/3+2*Ly/3 ! Red top band    
    end if
    ! fluam centers the unit cell around (0,0)
    y=y-Ly/2
    x=x-Lx/2
    write(*,*) x, y, 0.0
  end do
    
end program

