!
!######################################################################=$
!######################################################################=$
!
!                         PROGRAMME 2D CAVITY FLOW
!
!######################################################################=$
!######################################################################=$
!
!C========+=========+=========+=========+=========+=========+=========+=$
!C PROGRAM: cavity_flow
!C TYPE   : Main program
!C PURPOSE: Solve 2D cavity flow using Finite Volumes
!C========+=========+=========+=========+=========+=========+=========+=$
!
program cavity_flow
  implicit none
  !
  !定义变量
  !
  real*8, dimension(:)    , allocatable :: xc, yc, xf, yf, xv, yu, us, vs 
  real*8, dimension(:,:)  , allocatable :: u, v, p
  
  integer*8 :: Nx, Ny
  integer*8 :: imaille
  real*8    :: dt, tmax, Tf, ux0, Ro, rmu, Xmin, Xmax, Ymin, Ymax
  real*8    :: pi
  pi = dacos(-1.d0)
  !
  !读取参数
  !
  call read_input(Nx, Ny, dt, ux0, Ro, rmu, tmax, Xmin, Xmax, Ymin, Ymax, imaille)
  !
  !初始化变量
  !
  allocate(p(1:Nx,1:Ny), u(1:Nx+1,1:Ny+2), v(1:Nx+2,1:Ny+1))
  allocate(us(1:(Nx+1)*(Ny+2)), vs(1:(Nx+2)*(Ny+1)))
  allocate(xc(1:Nx), yc(1:Ny), xf(1:Nx+1), yf(1:Ny+1), xv(1:Nx+2), yu(1:Ny+2))
  !
  !设定网格
  !
  call set_mesh(xc, yc, xf, yf, xv, yu, Nx, Ny, Xmin, Xmax, Ymin, Ymax, pi, imaille)
  
  call solve_flow(Nx, Ny, dt, tmax, ux0, Ro, rmu, xc, yc, xf, yf, xv, yu, p, u, v, us, vs)
  
  deallocate(p, u, v, us, vs, xc, yc, xf, yf, xv, yu)
  
end program cavity_flow
subroutine read_input(Nx, Ny, dt, ux0, Ro, rmu, tmax, Xmin, Xmax, Ymin, Ymax, imaille)
  implicit none
  !
  integer*8 :: Nx, Ny
  integer*8 :: imaille, ivisu, pbtype
  real*8    :: dt, tmax, Xmin, Xmax, Ymin, Ymax, Ro, rmu, ux0

  write(*,*)
  write(*,*) "  ### Loading parameters ###"
  open(unit=18, file="parameters.dat", status="old")
  read(18,*) Nx
  read(18,*) Ny
  read(18,*) Xmin
  read(18,*) Xmax
  read(18,*) Ymin
  read(18,*) Ymax
  read(18,*) dt
  read(18,*) tmax
  read(18,*) Ro
  read(18,*) rmu
  read(18,*) ux0
  read(18,*) imaille
  close(18)
  write(*,*) "  ### Parameters loading : done ###"
  write(*,*) Nx, Ny, Xmin, Xmax, Ymin, Ymax, dt, tmax, Ro, rmu, ux0, imaille

end subroutine read_input

subroutine set_mesh(xc, yc, xf, yf, xv, yu, Nx, Ny, Xmin, Xmax, Ymin, Ymax, pi, imaille)
  implicit none
  !
  integer*8                 :: i, j, Nx, Ny, imaille
  real*8, dimension(1:Nx)   :: xc,yc
  real*8, dimension(1:Nx+1) :: xf,yf
  real*8, dimension(1:Nx+2) :: xv,yu
  real*8                    :: dx, dy, pi, Xmin, Xmax, Ymin, Ymax
  !
  ! Definition du maillage homogene :
  !----------------------------------
  if (imaille==0) then 
     dx=(Xmax-Xmin)/Nx
     dy=(Ymax-Ymin)/Ny
     do i=1,Nx
        xc(i)=Xmin+i*dx-0.5d0*dx;
     enddo
     do j=1,Ny
        yc(j)=Ymin+j*dy-0.5d0*dy;
     enddo
  endif
  !
  ! Definition du maillage inhomogene :
  !------------------------------------
  if (imaille==1) then
     do i=1,Nx
        xc(i)=(1.d0-dcos((i-0.5d0)*pi/Nx))/2.d0
     enddo
     do j=1,Ny
        yc(j)=(1.d0-dcos((j-0.5d0)*pi/Ny))/2.d0	 
     enddo
  endif

  if ((imaille .ne. 0).and.(imaille .ne. 1)) then
     write(*,*) "   ***   Mauvais choix de maillage   ***   "
     stop
  endif

  ! Definition des faces des volumes de controle :
  ! ----------------------------------------------
  do i=1,Nx
     xv(i+1)=xc(i)
  enddo
  do j=1,Ny
     yu(j+1)=yc(j)
  enddo
  xv(1)=2*xmin-xc(1)
  xv(Nx+2)=2*xmax-xc(Nx)
  yu(1)=2*ymin-yc(1)
  yu(Nx+2)=2*ymax-yc(Nx)
  do i=1,Nx+1
     xf(i)=(xv(i)+xv(i+1))/2.d0
  enddo
  do j=1,Ny+1
     yf(j)=(yu(j)+yu(j+1))/2.d0
  enddo
  write(*,*) "  ###   Mesh Setted   ###"
end subroutine set_mesh

subroutine solve_flow(Nx, Ny, dt, tmax, ux0, Ro, rmu, xc, yc, xf, yf, xv, yu, p, u, v, us, vs)
  use mod_inverse
  use mod_thomas
  implicit none
  !
  integer*8 :: it, Nt
  integer*8 :: Nx, Ny
  real*8    :: dt, tmax, ux0, Ro, rmu
  real*8, dimension(1:Nx+1,1:Ny+2)		:: u
  real*8, dimension(1:Nx+2,1:Ny+1)		:: v
  real*8, dimension(1:Nx)				:: xc
  real*8, dimension(1:Ny)				:: yc
  real*8, dimension(1:Nx+1)				:: xf, XBu, XAu
  real*8, dimension(1:Ny+1)				:: yf, YBv, YAv
  real*8, dimension(1:Nx+2)				:: xv, XBv
  real*8, dimension(1:Ny+2)				:: yu, YBu
  real*8, dimension(1:Nx,1:Ny)			:: p, Sp
  real*8, dimension(1:(Nx+1),1:(Ny+2))	:: us, Su
  real*8, dimension(1:(Nx+2),1:(Ny+1))	:: vs, Sv
  real*8, dimension(1:Nx*Ny)			:: fi
  real*8, dimension(0:Nx,0:Nx,0:Ny+1)		:: MuL,MuD,MuU
  real*8, dimension(0:Nx+1,0:Nx+1,0:Ny)		:: MvL,MvD,MvU
  real*8, dimension(0:Nx-1,0:Nx-1,0:Ny-1)	:: MpL,MpD,MpU
  real*8, dimension(1:(Nx+1)*(Ny+2),1:(Nx+1)*(Ny+2))	:: Mu
  
  
  write(*,*) "  ###   pre-calculating   ###"
  
  call cal_uvp_matrix(Nx, Ny, rmu, dt, Xc, Yc, Xv, Yu, Xf, Yf, MuL, MuD, MuU, MvL, MvD, MvU, MpL, MpD, MpU)
  
  
  call cal_XYAB(Nx, Ny, Xf, Yf, Xv, Yu, XBu, XBv, YBu, YBv, XAu, YAv)
  u=0.d0
  v=0.d0
  Nt=int(tmax/dt)
  do it=0,Nt
    call cal_uv_source(Nx, Ny, ux0, dt, u, v, XBu, XBv, YBu, YBv, Su, Sv)
!    write(*,*) "Su=["
!    write(*,*) Su
!    write(*,*) "]';Sv=["
!    write(*,*) Sv
	call RESOL_THOMAS(Nx,Ny+1,MuD,MuL,MuU,Su,us)
	call RESOL_THOMAS(Nx+1,Ny,MvD,MvL,MvU,Sv,vs)
 !   write(*,*) "]';us=["
 !   write(*,*) us
 !   write(*,*) "]';vs=["
 !   write(*,*) vs
 !   write(*,*) "]';"
    
	call cal_p_source(Nx, Ny, vs, us, xf, yf, Sp)
	call RESOL_THOMAS(Nx-1,Ny-1,MpD,MpL,MpU,Sp,fi)
	call cal_uv(Nx, Ny, XAu, YAv, us, vs, u, v, fi)

	if (mod(it,int(Nt/100)) .eq. 0) then
	  call save_VTK_uvp(u, v, p, fi, dt, Ro, xc, yc, it, Nx, Ny)
      write(*,*) "   ***   time=", it*dt, " calculated   ***   "
	end if
  enddo
end subroutine solve_flow