subroutine cal_uv_source(Nx, Ny, ux0, dt, u, v, XBu, XBv, YBu, YBv, Su, Sv)
  implicit none

  integer*8 :: Nx, Ny, i, j, jj, nnx
  real*8  :: dt, ux0, q
  real*8, dimension(1:Nx+1,1:Ny+2)		:: u
  real*8, dimension(1:Nx+2,1:Ny+1)		:: v
  real*8, dimension(1:Nx+1)				:: XBu
  real*8, dimension(1:Ny+1)				:: YBv
  real*8, dimension(1:Nx+2)				:: XBv
  real*8, dimension(1:Ny+2)				:: YBu
  real*8, dimension(1:(Nx+1),1:(Ny+2))	:: Su
  real*8, dimension(1:(Nx+2),1:(Ny+1))	:: Sv
  q=1.d0/dt
  !首先应用边界条件
  do i=1,Nx+1
    u(i,Ny+2)=2.d0*ux0-u(i,Ny+1)
    u(i,1)=-u(i,2)
  enddo
  do j=1,Ny+1
    v(Nx+2,j)=-v(Nx+1,j)
    v(1,j)=-v(2,j)
  enddo
  !计算u方程源矩阵
  Su=0.d0
  !nnx=Nx+1
  !!计算内部情况
  do i=2,Nx
    do jj=2,Ny+1
      j=jj-1
      Su(i,jj)=-u(i,jj)*XBu(i)*(u(i+1,jj)-u(i-1,jj)) & 
             & -0.25d0*(v(i,jj-1)+v(i,jj)+v(i+1,jj-1)+v(i+1,jj))*YBu(jj)*(u(i,jj+1)-u(i,jj-1)) &
             & +q*u(i,jj)
    enddo
  enddo
  !!计算上下边界
  do i=1,Nx+1
    !!!下边界
    Su(i,1)=0.d0
    !!!上边界
    Su(i,Ny+2)=2.d0*ux0
  enddo
  !!计算左右边界
  do jj=2,Ny+1
    j=jj-1
    Su(1,jj)=0.d0
    Su(Nx+1,jj)=0.d0
  enddo
  
  
  !计算v方程源矩阵
  Sv=0.d0
  !nnx=Nx+2
  !!计算内部情况
  do i=2,Nx+1
    do jj=2,Ny
      j=jj-1
      Sv(i,jj)=-0.25d0*(u(i-1,jj)+u(i-1,jj+1)+u(i,jj)+u(i,jj+1))*XBv(i)*(v(i+1,jj)-v(i-1,jj)) & 
             & -v(i,jj)*YBv(jj)*(v(i,jj+1)-v(i,jj-1)) &
             & +q*v(i,jj)
    enddo
  enddo
  !!计算上下边界
  do i=2,Nx+1
    !!!下边界
    Sv(i,1)=0.d0
    !!!上边界
    Sv(i,Ny+1)=0.d0
  enddo
  !!计算左右边界
  do jj=1,Ny+1
    j=jj-1
    Sv(1,jj)=0.d0
    Sv(Nx+2,jj)=0.d0
  enddo
end subroutine cal_uv_source