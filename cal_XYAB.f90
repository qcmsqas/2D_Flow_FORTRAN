subroutine cal_XYAB(Nx, Ny, Xf, Yf, Xv, Yu, XBu, XBv, YBu, YBv, XAu, YAv)
  implicit none
  integer*8 :: Nx, Ny, i
  real*8, dimension(1:Nx+1)				:: xf, XBu, XAu
  real*8, dimension(1:Ny+1)				:: yf, YBv, YAv
  real*8, dimension(1:Nx+2)				:: xv, XBv
  real*8, dimension(1:Ny+2)				:: yu, YBu
  XBu = 0.d0
  XAu = 0.d0
  XBv = 0.d0
  YBu = 0.d0
  YBv = 0.d0
  YAv = 0.d0
  do i=2,Nx
    XBu(i)=1.d0/(Xf(i+1)-Xf(i-1))
  enddo
  do i=2,Nx+1
    XBv(i)=1.d0/(Xv(i+1)-Xv(i-1))
  enddo
  do i=2,Ny+1
    YBu(i)=1.d0/(Yu(i+1)-Yu(i-1))
  enddo
  do i=2,Ny
    YBv(i)=1.d0/(Yf(i+1)-Yf(i-1))
  enddo
  do i=2,Nx
    XAu(i)=1.d0/(Xv(i+1)-Xv(i))
  enddo
  do i=2,Ny
    YAv(i)=1.d0/(Yu(i+1)-Yu(i))
  enddo

end subroutine cal_XYAB