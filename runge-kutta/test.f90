program main
  implicit none
  integer::i,N,step,Nmax
  double precision::x,h
  double precision,allocatable::y(:)
  double precision,external::grk

  !==============
  h=1.d-3
  Nmax=20000
  step=200
  N=2   ! +-------+
  !==============

  x=0d0
  allocate(y(1:N))
  y(1)=1d0
  y(2)=-0.15d0   ! +-------+

  write(10,*)x,y(1),y(2)
  do i=1,Nmax
     call rk4(grk,N,x,h,y)
     if(mod(i,step).eq.0)then
        write(10,*)x,y(1),y(2)
     endif
  enddo

  stop
end program main

function grk(N,x,y,s)
  implicit none
  integer,intent(in)::N,s
  double precision,intent(in)::x,y(1:N)
  double precision::grk

  grk=0d0
  if(s.eq.1)then
     grk=y(2)
  elseif(s.eq.2)then        ! +-------+
     grk=-0.3d0*y(2)-y(1)   ! +-------+
  else
     write(6,*)"****Error at grk"; stop
  endif

  return
end function grk

!===================================

subroutine rk4(grk,N,x,h,y)
  implicit none
  integer,intent(in)::N
  double precision,intent(in)::h
  double precision,intent(inout)::x,y(1:N)
  double precision,external::grk
  !   sikinote(http://slpr.sakura.ne/jp/qp)
  ! N  --> number of differential equation
  ! x  --> variable
  ! h  --> step size
  ! y  --> solution
  integer::i,j
  double precision,allocatable::tmp(:),K(:,:)
  double precision::tx,c(1:4)

  c(1:4)=(/0d0,0.5d0,0.5d0,1d0/)
  allocate(tmp(1:N),K(1:N,0:4))
  tmp=0d0; K=0d0

  do j=1,4
     tx = x + c(j)*h
     do i=1,N
        tmp(i) = y(i) + K(i,j-1)*c(j)
     enddo
     do i=1,N
        K(i,j) = h*grk(N,tx,tmp,i)
     enddo
  enddo

  x=x+h
  do i=1,N
     y(i)=y(i)+(K(i,1)+K(i,4))/6d0+(K(i,2)+K(i,3))/3d0
  enddo

  return
end subroutine rk4
