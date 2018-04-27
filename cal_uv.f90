subroutine cal_uv(Nx, Ny, XAu, YAv, us, vs, u, v, fi)
  implicit none
  integer*8 :: Nx, Ny, nnx, i, j, jj
  real*8, dimension(1:Nx+1,1:Ny+2)		:: u
  real*8, dimension(1:Nx+2,1:Ny+1)		:: v
  real*8, dimension(1:Nx+1)				:: XAu
  real*8, dimension(1:Ny+1)				:: YAv
  real*8, dimension(1:(Nx+1),1:(Ny+2))	:: us
  real*8, dimension(1:(Nx+2),1:(Ny+1))	:: vs
  real*8, dimension(1:Nx*Ny)			:: fi
  
  nnx=Nx+1
  do i=2,Nx
    do jj=2,Ny+1
      j=jj-1
      u(i,jj)=us(i,jj)-(fi(i+(j-1)*Nx)-fi(i-1+(j-1)*Nx))*XAu(i)
    enddo
  enddo
  nnx=Nx+2
  do i=2,Nx+1
    do jj=2,Ny
      j=jj-1
      v(i,jj)=vs(i,jj)-(fi(i-1+j*Nx)-fi(i-1+(j-1)*Nx))*YAv(jj)
    enddo
  enddo
	
end subroutine cal_uv