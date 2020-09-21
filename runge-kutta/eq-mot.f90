! 4th Order Runge-Kutta program
! created by Daichi Hayashi 29 Sep. 2019.

program eq_of_mot
  implicit none
  double precision t0,x0,v0,h,t,x,v,m,k,g,k11,k12,k21,k22,k31,k32,k41,k42,f1,f2
  integer i
  t0 = 0.d0
  x0 = 1.d0
  v0 = -0.15d0
  h  = 0.05d0
  t  = t0
  x  = x0
  v  = v0
  k  = 0.15d0

  open(200,file='out.dat')
  write(200,100) t,x,g(t,k)
  do i = 1, 300
     k11 = h*f1(v)
     k12 = h*f2(x,v,k)

     k21 = h*f1(v+h*k12/2.0d0)
     k22 = h*f2(x+h*k11/2.0d0,v+h*k12/2.0d0,k)

     k31 = h*f1(v+h*k22/2.0d0)
     k32 = h*f2(x+h*k21/2.0d0,v+h*k22/2.0d0,k)

     k41 = h*f1(v+h*k32)
     k42 = h*f2(x+h*k31,v+h*k32,k)

     t  = t + h
     x  = x + (k11 + 2.0d0*k21 + 2.0d0*k31 + k41)/6.0d0
     v  = v + (k12 + 2.0d0*k22 + 2.0d0*k32 + k42)/6.0d0
     write(200,100) t,x,g(t,k)
  end do

  close(200)
100 format (10e12.4)
end program eq_of_mot

function f1(v)
  ! 1th differential eqation,
  double precision f1,v
  f1 = v
end function f1

function f2(x,v,k)
  ! 1th differential eqation,
  double precision f2,x,k
  f2 = -2.0d0*k*v-x
end function f2

function g(t,k)
  ! exact function
  double precision g,k,m,t
  g = dexp(-k*t)*dcos(dsqrt(1.0d0-k**2d0)*t)
end function g
