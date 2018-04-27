subroutine save_VTK_uvp(u, v, p, fi, dt, Ro, xc, yc, it, Nx, Ny)
  implicit none

  integer*8 :: it, i, j
  integer*8 :: Nx, Ny
  real*8    :: dt, Ro, rodt
  real*8, dimension(1:Nx+1,1:Ny+2)		:: u
  real*8, dimension(1:Nx+2,1:Ny+1)		:: v
  real*8, dimension(1:Nx)				:: xc
  real*8, dimension(1:Ny)				:: yc
  real*8, dimension(0:Nx-1,0:Ny-1)		:: p, uout, vout
  real*8, dimension(1:Nx*Ny)			:: fi
  
  
  uout = (u(1:Nx,2:Ny+1)+u(2:Nx+1,2:Ny+1))*0.5d0
  vout = (v(2:Nx+1,1:Ny)+v(2:Nx+1,2:Ny+1))*0.5d0
  p = 0.d0
  rodt=ro/dt;
  do i=1,Nx
    do j=1,Ny
      p(i-1,j-1)=fi(i+(j-1)*Nx)*rodt;
    enddo
  enddo
  call save_VTK_trainee(uout, vout, p, xc, yc, it, Nx-1, Ny-1)

end subroutine save_VTK_uvp