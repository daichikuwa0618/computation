! 4th Order Runge-Kutta program
! created by Daichi Hayashi 29 Sep. 2019.

program runge_kutta_4th
  double precision x0,y0,x,y,h,k1,k2,k3,k4,f,g
  x0 = 0.d0
  y0 = 1.d0
  h  = 0.1d0
  x  = x0
  y  = y0

  write(6,100) x,y,g(x)
  do i = 1, 1000
     k1 = h*f(x,y)
     k2 = h*f(x+h/2.0d0,y+k1/2.0d0)
     k3 = h*f(x+h/2.0d0,y+k2/2.0d0)
     k4 = h*f(x+h,y+k3)
     x  = x + h
     y  = y + (k1+2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
     write(6,*) x,y,g(x)
  end do

100 format (f10.7, 1x, f10.7, 1x, f10.7)
end program runge_kutta_4th

function f(x,y)
  ! 1th differential eqation,
  double precision f,x,y
  f = x - y
end function f

function g(x)
  ! exact function
  double precision g,x
  g = 2.0d0*dexp(-x) + x - 1.0d0
end function g
