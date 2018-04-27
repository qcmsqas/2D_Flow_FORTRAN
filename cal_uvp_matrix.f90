subroutine cal_uvp_matrix(Nx, Ny, rmu, dt, Xc, Yc, Xv, Yu, Xf, Yf, MuL, MuD, MuU, MvL, MvD, MvU, MpL, MpD, MpU)
  use mod_inverse
  use mod_thomas
  implicit none
  
  integer*8 :: Nx, Ny, nnx, i, j, jj, ip, iq, jp, jq,stx,edx,sty,edy
  real*8  :: rmu, dt, q
  real*8, dimension(1:Nx)				:: xc
  real*8, dimension(1:Ny)				:: yc
  real*8, dimension(1:Nx+1)				:: xf
  real*8, dimension(1:Ny+1)				:: yf
  real*8, dimension(1:Nx+2)				:: xv
  real*8, dimension(1:Ny+2)				:: yu
  real*8, dimension(1:(Nx+1)*(Ny+2),1:(Nx+1)*(Ny+2))	:: Mu
  real*8, dimension(1:(Nx+2)*(Ny+1),1:(Nx+2)*(Ny+1))	:: Mv
  real*8, dimension(1:(Nx  )*(Ny  ),1:(Nx  )*(Ny  ))	:: Mp
!  real*8, dimension(1:(Nx  )*(Ny  ),1:(Nx+1)*(Ny+2))	:: Mpu, MMpu
!  real*8, dimension(1:(Nx  )*(Ny  ),1:(Nx+2)*(Ny+1))	:: Mpv, MMpv
  real*8, dimension(0:Nx  ,0:Nx  ,0:Ny+1)		:: MuL,MuD,MuU
  real*8, dimension(0:Nx+1,0:Nx+1,0:Ny  )		:: MvL,MvD,MvU
  real*8, dimension(0:Nx-1,0:Nx-1,0:Ny-1)		:: MpL,MpD,MpU

  
  q=1.d0/dt
  !!获得u更新us的泊松方程的方程矩阵。
  Mu=0.d0
  nnx=Nx+1
  !!考虑内部点
  do i=2,Nx
    do jj=2,Ny+1
      j=jj-1
      Mu(i+j*nnx,i+j*nnx)=Mu(i+j*nnx,i+j*nnx)+(q &
                       & +rmu/(Xf(i)-Xf(i-1))/(Xc(i)-Xc(i-1))&
                       & +rmu/(Xf(i+1)-Xf(i))/(Xc(i)-Xc(i-1))&
                       & +rmu/(Yu(jj)-Yu(jj-1))/(Yf(jj)-Yf(jj-1))&
                       & +rmu/(Yu(jj+1)-Yu(jj))/(Yf(jj)-Yf(jj-1))&
                       & )
      Mu(i+j*nnx,(i-1)+j*nnx)=Mu(i+j*nnx,(i-1)+j*nnx)-rmu/(Xf(i)-Xf(i-1))/(Xc(i)-Xc(i-1))
      Mu(i+j*nnx,(i+1)+j*nnx)=Mu(i+j*nnx,(i+1)+j*nnx)-rmu/(Xf(i+1)-Xf(i))/(Xc(i)-Xc(i-1))
      Mu(i+j*nnx,i+(j-1)*nnx)=Mu(i+j*nnx,i+(j-1)*nnx)-rmu/(Yu(jj)-Yu(jj-1))/(Yf(jj)-Yf(jj-1))
      Mu(i+j*nnx,i+(j+1)*nnx)=Mu(i+j*nnx,i+(j+1)*nnx)-rmu/(Yu(jj+1)-Yu(jj))/(Yf(jj)-Yf(jj-1))
    enddo
  enddo
  !!考虑上下边界
  do i=1,Nx+1
    !!!下边界
    Mu(i,i)=Mu(i,i)+1.d0
    Mu(i,i+nnx)=Mu(i,i+nnx)+1.d0
    !!!上边界
    Mu(i+(Ny+1)*nnx,i+(Ny+1)*nnx)=Mu(i+(Ny+1)*nnx,i+(Ny+1)*nnx)+1.d0
    Mu(i+(Ny+1)*nnx,i+(Ny)*nnx)=Mu(i+(Ny+1)*nnx,i+(Ny)*nnx)+1.d0
  enddo
  !!考虑左右边界
  do jj=2,Ny+1
    j=jj-1
    Mu(1+j*nnx,1+j*nnx)=Mu(1+j*nnx,1+j*nnx)+1.d0
    Mu(Nx+1+j*nnx,Nx+1+j*nnx)=Mu(Nx+1+j*nnx,Nx+1+j*nnx)+1.d0
  enddo
  
  !!将矩阵转换为分块三对角表示
  do i=0,Ny+1
	stx=i*(Nx+1)+1
	sty=i*(Nx+1)+1
	edx=i*(Nx+1)+Nx+1
	edy=i*(Nx+1)+Nx+1
    MuD(:,:,i)=Mu(stx:edx,sty:edy)
  enddo
  MuL(:,:,0)=0;
  do i=0,Ny
	stx=(i+1)*(Nx+1)+1
	sty=i*(Nx+1)+1
	edx=(i+1)*(Nx+1)+Nx+1
	edy=i*(Nx+1)+Nx+1
	MuL(:,:,i+1)=Mu(stx:edx,sty:edy)
  enddo
  MuU(:,:,Ny+1)=0;
  do i=0,Ny
	stx=i*(Nx+1)+1
	sty=(i+1)*(Nx+1)+1
	edx=i*(Nx+1)+Nx+1
	edy=(i+1)*(Nx+1)+Nx+1
	MuU(:,:,i)=Mu(stx:edx,sty:edy)
  enddo
  !Mu=Mu^(-1)
  !call inverse(Mu,(Nx+1)*(Ny+2))
  !write(*,*) "Mu=["
  !write(*,*) reshape(Mu,(/(Nx+1)**2*(Ny+2)**2,1/))
  !write(*,*) "]';MuD=["
  !write(*,*) reshape(MuD,(/(Nx+1)**2*(Ny+2),1/))
!  write(*,*) "]';MuL=["
 ! write(*,*) reshape(MuL,(/(Nx+1)**2*(Ny+2),1/))
  !write(*,*) "]';MuU=["
  !write(*,*) reshape(MuU,(/(Nx+1)**2*(Ny+2),1/))
  call INIT_THOMAS(Nx,Ny+1,MuD,MuL,MuU)
  !write(*,*) "]';MD=["
  !write(*,*) reshape(MuD,(/(Nx+1)**2*(Ny+2),1/))
!  write(*,*) "]';ML=["
!  write(*,*) reshape(MuL,(/(Nx+1)**2*(Ny+2),1/))
!  write(*,*) "]';MU=["
!  write(*,*) reshape(MuU,(/(Nx+1)**2*(Ny+2),1/))
!  write(*,*) "]';"
  
  
  !获得v更新vs的泊松方程的方程矩阵
  Mv=0.d0
  nnx=Nx+2
  !!考虑内部点
  do i=2,Nx+1
    do jj=2,Ny
      j=jj-1
      Mv(i+j*nnx,i+j*nnx)=Mv(i+j*nnx,i+j*nnx)+(q &
                       & +rmu/(Xv(i)-Xv(i-1))/(Xf(i)-Xf(i-1)) &
                       & +rmu/(Xv(i+1)-Xv(i))/(Xf(i)-Xf(i-1)) &
                       & +rmu/(Yf(jj)-Yf(jj-1))/(Yc(jj)-Yc(jj-1)) &
                       & +rmu/(Yf(jj+1)-Yf(jj))/(Yc(jj)-Yc(jj-1)) &
                       & )
      Mv(i+j*nnx,(i-1)+j*nnx)=Mv(i+j*nnx,(i-1)+j*nnx)-rmu/(Xv(i)-Xv(i-1))/(Xf(i)-Xf(i-1))
      Mv(i+j*nnx,(i+1)+j*nnx)=Mv(i+j*nnx,(i+1)+j*nnx)-rmu/(Xv(i+1)-Xv(i))/(Xf(i)-Xf(i-1))
      Mv(i+j*nnx,i+(j-1)*nnx)=Mv(i+j*nnx,i+(j-1)*nnx)-rmu/(Yf(jj)-Yf(jj-1))/(Yc(jj)-Yc(jj-1))
      Mv(i+j*nnx,i+(j+1)*nnx)=Mv(i+j*nnx,i+(j+1)*nnx)-rmu/(Yf(jj+1)-Yf(jj))/(Yc(jj)-Yc(jj-1))
    enddo
  enddo
  !!考虑左右边界
  do jj=1,Ny+1
    j=jj-1
    !!!左边界
    Mv(1+j*nnx,1+j*nnx)=Mv(1+j*nnx,1+j*nnx)+1.d0
    Mv(1+j*nnx,2+j*nnx)=Mv(1+j*nnx,2+j*nnx)+1.d0
    !!!!右边界
    Mv(Nx+2+j*nnx,Nx+2+j*nnx)=Mv(Nx+2+j*nnx,Nx+2+j*nnx)+1.d0
    Mv(Nx+2+j*nnx,Nx+1+j*nnx)=Mv(Nx+2+j*nnx,Nx+1+j*nnx)+1.d0
  enddo
  !!考虑上下边界
  do i=2,Nx+1
    !!!下边界
    Mv(i,i)=Mv(i,i)+1.d0
    !!!上边界
    Mv(i+Ny*nnx,i+Ny*nnx)=Mv(i+Ny*nnx,i+Ny*nnx)+1.d0
  enddo
  
  
  do i=0,Ny
	stx=i*(Nx+2)+1
	sty=i*(Nx+2)+1
	edx=i*(Nx+2)+Nx+2
	edy=i*(Nx+2)+Nx+2
    MvD(:,:,i)=Mv(stx:edx,sty:edy)
  enddo
  MvL(:,:,0)=0;
  do i=0,Ny-1
	stx=(i+1)*(Nx+2)+1
	sty=i*(Nx+2)+1
	edx=(i+1)*(Nx+2)+Nx+2
	edy=i*(Nx+2)+Nx+2
	MvL(:,:,i+1)=Mv(stx:edx,sty:edy)
  enddo
  MvU(:,:,Ny)=0;
  do i=0,Ny-1
	stx=i*(Nx+2)+1
	sty=(i+1)*(Nx+2)+1
	edx=i*(Nx+2)+Nx+2
	edy=(i+1)*(Nx+2)+Nx+2
	MvU(:,:,i)=Mv(stx:edx,sty:edy)
  enddo
  !Mv=Mv^(-1)
  !call inverse(Mv,(Nx+2)*(Ny+1))
  call INIT_THOMAS(Nx+1,Ny,MvD,MvL,MvU)
  
  
  !计算p的泊松方程的方程矩阵
  Mp=0.d0
  !!考虑内部点
  do i=1,Nx
    do jj=1,Ny
      j=jj-1
      !disp([num2str(i) ' ' num2str(jj)])
      !处理四种边界
      ip=1
	  iq=1
	  jp=1
	  jq=1
      !必须有一个狄利克雷边界，否则会不稳定
      if ((i .eq. 1) .and. (jj .eq. 1)) then
        Mp(i+j*Nx,i+j*Nx)=Mp(i+j*Nx,i+j*Nx)+4.d0/(Xv(2)-Xv(1))/(Yu(2)-Yu(1))
        continue
      endif
      if (i .eq. 1) then
        ip=0
      endif
      if (i .eq. Nx) then
        iq=0
      endif
      if (jj .eq. 1) then
        jp=0
      endif
      if (jj .eq. Ny) then
        jq=0
      endif
      Mp(i+j*Nx,i+j*Nx)=Mp(i+j*Nx,i+j*Nx)-1.d0/(Xv(i+1)-Xv(i))/(Xf(i+1)-Xf(i)) &
                         -1.d0/(Xv(i+2)-Xv(i+1))/(Xf(i+1)-Xf(i)) &
                         -1.d0/(Yu(jj+1)-Yu(jj))/(Yf(jj+1)-Yf(jj)) &
                         -1.d0/(Yu(jj+2)-Yu(jj+1))/(Yf(jj+1)-Yf(jj))
      Mp(i+j*Nx,(i-ip)+j*Nx)=Mp(i+j*Nx,(i-ip)+j*Nx)+1.d0/(Xv(i+1)-Xv(i))/(Xf(i+1)-Xf(i))
      Mp(i+j*Nx,(i+iq)+j*Nx)=Mp(i+j*Nx,(i+iq)+j*Nx)+1.d0/(Xv(i+2)-Xv(i+1))/(Xf(i+1)-Xf(i))
      Mp(i+j*Nx,i+(j-jp)*Nx)=Mp(i+j*Nx,i+(j-jp)*Nx)+1.d0/(Yu(jj+1)-Yu(jj))/(Yf(jj+1)-Yf(jj))
      Mp(i+j*Nx,i+(j+jq)*Nx)=Mp(i+j*Nx,i+(j+jq)*Nx)+1.d0/(Yu(jj+2)-Yu(jj+1))/(Yf(jj+1)-Yf(jj))
    enddo
  enddo
  !Mp=Mp^(-1)
  !call inverse(Mp,Nx*Ny)
  
  do i=0,Ny-1
	stx=i*(Nx)+1
	sty=i*(Nx)+1
	edx=i*(Nx)+Nx
	edy=i*(Nx)+Nx
    MpD(:,:,i)=Mp(stx:edx,sty:edy)
  enddo
  MpL(:,:,0)=0;
  do i=0,Ny-2
	stx=(i+1)*(Nx)+1
	sty=i*(Nx)+1
	edx=(i+1)*(Nx)+Nx
	edy=i*(Nx)+Nx
	MpL(:,:,i+1)=Mp(stx:edx,sty:edy)
  enddo
  MpU(:,:,Ny-1)=0;
  do i=0,Ny-2
	stx=i*(Nx)+1
	sty=(i+1)*(Nx)+1
	edx=i*(Nx)+Nx
	edy=(i+1)*(Nx)+Nx
	MpU(:,:,i)=Mp(stx:edx,sty:edy)
  enddo

  call INIT_THOMAS(Nx-1,Ny-1,MpD,MpL,MpU)
  
  !MMpu=0.d0;
  !MMpv=0.d0;
  !do i=1,Nx
  !  do jj=1,Ny
  !    !为了方程稳定性，必须有一个狄利克雷边条
  !    if ((i .eq. 1) .and. (jj .eq. 1)) then
  !      continue
  !    endif
  !    j=jj-1
  !    MMpu(i+j*Nx,i+(j+1)*(Nx+1))=-1.d0/(Xf(i+1)-Xf(i))
  !    MMpu(i+j*Nx,i+1+(j+1)*(Nx+1))=1.d0/(Xf(i+1)-Xf(i))
  !    MMpv(i+j*Nx,i+1+j*(Nx+2))=-1.d0/(Yf(jj+1)-Yf(jj))
  !    MMpv(i+j*Nx,i+1+(j+1)*(Nx+2))=1.d0/(Yf(jj+1)-Yf(jj))
  !  enddo
  !enddo
  !Mpu=matmul(Mp,MMpu)
  !Mpv=matmul(Mp,MMpv)
  
  
end subroutine cal_uvp_matrix