subroutine cal_p_source(Nx, Ny, vs, us, xf, yf, Sp)

  implicit none

  integer*8 :: Nx, Ny, i, j
  real*8, dimension(1:Nx+1)				:: xf
  real*8, dimension(1:Ny+1)				:: yf
  real*8, dimension(1:Nx,1:Ny)			:: Sp
  real*8, dimension(1:(Nx+1),1:(Ny+2))	:: us
  real*8, dimension(1:(Nx+2),1:(Ny+1))	:: vs
  
  do i=1,Nx
    do j=1,Ny
	  Sp(i,j)=(us(i+1,j+1)-us(i,j+1))/(xf(i+1)-xf(i))+(vs(i+1,j+1)-vs(i+1,j))/(yf(j+1)-yf(j))
    enddo
  enddo

end subroutine cal_p_source