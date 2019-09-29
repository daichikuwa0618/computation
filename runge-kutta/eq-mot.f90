! 4th Order Runge-Kutta program
! created by Daichi Hayashi 29 Sep. 2019.

program eq_of_mot
  double precision t0,x0,v0,h,t,x,v,m,k,g,k11,k12,k21,k22,k31,k32,k41,k42,f1,f2
  t0 = 0.d0
  x0 = 1.d0
  v0 = 0.d0
  h  = 0.01d0
  t  = t0
  x  = x0
  v  = v0
  k  = 1.0d0
  m  = 1.0d0

  open(200,file='out.dat')
  write(200,*) t,x,g(t,m,k)
  do i = 1, 1000
     k11 = h*f1(v)
     k12 = h*f2(x,m,k)

     k21 = h*f1(v+k12/2.0d0)
     k22 = h*f2(x+k11/2.0d0,m,k)

     k31 = h*f1(v+k22/2.0d0)
     k32 = h*f2(x+k21/2.0d0,m,k)

     k41 = h*f1(v+k32)
     k42 = h*f2(x+k31,m,k)

     t  = t + h
     x  = x + (k11+2.0d0*k21 + 2.0d0*k31 + k41)/6.0d0
     v  = v + (k12+2.0d0*k22 + 2.0d0*k32 + k42)/6.0d0
     write(200,*) t,x,g(t,m,k)
  end do

  close(200)
100 format (f10.7, 1x, f10.7, 1x, f10.7)
end program eq_of_mot

function f1(v)
  ! 1th differential eqation,
  double precision f1,v
  f1 = v
end function f1

function f2(x,m,k)
  ! 1th differential eqation,
  double precision f2,x,m,k
  f2 = -k/m*x
end function f2

function g(t,m,k)
  ! exact function
  double precision g,k,m,t
  g = dcos(dsqrt(k/m)*t)
end function g
