program initial
  ! Creates a random suspension of ideal particles with half green and half red, phase separated
  integer, parameter :: n=4096
  integer, parameter :: ng=n/2
  integer :: p
  real :: x, y, z, Lx, Ly
  
  Lx=73.5375 ! Periodic dimension
  Ly=Lx ! Periodic dimension
  z=5.0 ! Height above wall
  
  write(*,*) n
  do p=1,n
    call random_number(x)
    call random_number(y)
    x=x*Lx
    if(p<ng) then ! Green
       y=y*Ly/2
    else 
       y=y*Ly/2+Ly/2 ! Red top band    
    end if
    ! fluam centers the unit cell around (0,0)
    y=y-Ly/2
    x=x-Lx/2
    ! Note that we write quaternions here even though rotation is not tracked explicitly
    write(*,*) x, y, z, " 1.0 0.0 0.0 0.0"
  end do
    
end program
