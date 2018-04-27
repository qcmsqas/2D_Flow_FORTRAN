!
!C========+=========+=========+=========+=========+=========+=========+=$
!C PROGRAM: save_VTK_TRAINEE 
!C TYPE   :
!C PURPOSE: save fields with index itsave with VTK format
!C========+=========+=========+=========+=========+=========+=========+=$
!
subroutine save_VTK_trainee(velx, vely, press, x, y, itsave, Nx, Ny)
  implicit none
  !
  integer*8 Nx, Ny, i, j, itsave
  real*8, dimension(0:Nx,0:Ny) :: velx, vely, press
  real*8, dimension(0:Nx)      :: x
  real*8, dimension(0:Ny)      :: y

  character*5 numarch
  character*30 nomfichier
  !
  write(*,*)
  write(numarch,'(I5.5)') itsave
  write(*,*)' ## saving fields VTK n=',numarch,' ##'
  nomfichier=adjustl('champs_'//numarch//'.vtk')
  open(unit=20, file=nomfichier, status='unknown')
  !
  write(20,'(A26)') '# vtk DataFile Version 2.0'
  write(20,'(A19)') 'Lid driven velocity'
  write(20,'(A5)') 'ASCII'
  write(20,'(A24)') 'DATASET RECTILINEAR_GRID'
  write(20,'((A10),3((I5)))') 'DIMENSIONS',Nx+1,Ny+1,1
  write(20,'((A13),(I5),X,(A5))') 'X_COORDINATES',Nx+1,'float'
  do i=0,Nx
     write(20,'(e13.6)') x(i)
  enddo
  write(20,'((A13),(I5),X,(A5))') 'Y_COORDINATES',Ny+1,'float'
  do j=0,Ny
     write(20,'(e13.6)') y(j)
  enddo
  write(20,'((A13),(I2),X,(A5))') 'Z_COORDINATES',1,'float'
  write(20,'(e13.6)') 0.
  write(20,'((A10),X,(I5))') 'POINT_DATA',(Nx+1)*(Ny+1)
  write(20,*) 'VECTORS Velocity float '
  do j=0,Ny
    do i=0,Nx
        write(20,'(3(e13.6,X))') velx(i,j),vely(i,j),0.
    enddo
  enddo
  write(20,*) 'SCALARS Pression float 1'
  write(20,*) 'LOOKUP_TABLE default'
  do j=0,Ny
    do i=0,Nx
        write(20,'(e13.6)') press(i,j)
    enddo
  enddo
  close(20)
  !
end subroutine save_VTK_trainee

!C========+=========+=========+=========+=========+=========+=========+=$
!C  save_VTK_CHAUFFEE
!C========+=========+=========+=========+=========+=========+=========+=$
subroutine save_VTK_chauffee(temp, velx, vely, press, x, y, itsave, Nx, Ny)
  implicit none
  !
  integer*8 Nx, Ny, i, j, itsave
  real*8, dimension(0:Nx,0:Ny) :: temp, velx, vely, press
  real*8, dimension(0:Nx)      :: x
  real*8, dimension(0:Ny)      :: y

  character*5 numarch
  character*30 nomfichier
  !
  write(*,*)
  write(numarch,'(I5.5)') itsave
  write(*,*)' ## saving fields VTK n=',numarch,' ##'
  nomfichier=adjustl('champs_'//numarch//'.vtk')
  open(unit=20, file=nomfichier, status='unknown')
  !
  write(20,'(A26)') '# vtk DataFile Version 2.0'
  write(20,'(A19)') 'Lid driven velocity'
  write(20,'(A5)') 'ASCII'
  write(20,'(A24)') 'DATASET RECTILINEAR_GRID'
  write(20,'((A10),3((I5)))') 'DIMENSIONS',Nx+1,Ny+1,1
  write(20,'((A13),(I5),X,(A5))') 'X_COORDINATES',Nx+1,'float'
  do i=0,Nx
     write(20,'(e13.6)') x(i)
  enddo
  write(20,'((A13),(I5),X,(A5))') 'Y_COORDINATES',Ny+1,'float'
  do j=0,Ny
     write(20,'(e13.6)') y(j)
  enddo
  write(20,'((A13),(I2),X,(A5))') 'Z_COORDINATES',1,'float'
  write(20,'(e13.6)') 0.
  write(20,'((A10),X,(I5))') 'POINT_DATA',(Nx+1)*(Ny+1)
  write(20,*) 'SCALARS Temperature float 1'
  write(20,*) 'LOOKUP_TABLE default'
  do j=0,Ny
    do i=0,Nx
        write(20,'(e13.6)') temp(i,j)
    enddo
  enddo
  write(20,*) 'VECTORS Velocity float '
  do j=0,Ny
    do i=0,Nx
        write(20,'(3(e13.6,X))') velx(i,j),vely(i,j),0.
    enddo
  enddo
  write(20,*) 'SCALARS Pressure float 1'
  write(20,*) 'LOOKUP_TABLE default'
  do j=0,Ny
    do i=0,Nx
        write(20,'(e13.6)') press(i,j)
    enddo
  enddo
  close(20)
  !

end subroutine save_VTK_chauffee

!Subroutine pour la visualisation des opérateurs
!------------------------------------------------------------------------------
subroutine visualisation(D,L,U,Nx,Ny,matname)
  use ISO_FORTRAN_ENV

  integer :: Nx, Ny,j,i
  real*8, dimension (0:Nx,0:Nx,0:Ny) :: D, L, U
  character*5 matname
  integer ::  unit_descriptor

!  unit_descriptor=55          !ecriture dans un fichier fort.55
  unit_descriptor=OUTPUT_UNIT !ecriture dans la console

  write(unit_descriptor,*) '==========='
  write(unit_descriptor,*) 'Operateur D',matname
  write(unit_descriptor,*) '==========='
  do j=0,Ny
     do i=0,Nx
        write(unit_descriptor,1000) D(i,:,j)
     enddo
     write(unit_descriptor,*)
  enddo

  write(unit_descriptor,*) '==========='
  write(unit_descriptor,*) 'Operateur L',matname
  write(unit_descriptor,*) '==========='
  do j=0,Ny
     do i=0,Nx
        write(unit_descriptor,1000) L(i,:,j)
     enddo
     write(unit_descriptor,*)
  enddo

  write(unit_descriptor,*) '==========='
  write(unit_descriptor,*) 'Operateur U',matname
  write(unit_descriptor,*) '==========='
  do j=0,Ny
     do i=0,Nx
        write(unit_descriptor,1000) U(i,:,j)
     enddo
     write(unit_descriptor,*)
  enddo

! 1000 format(5(1x,E16.9))  ! si ecriture dans un fichier fort.55
 1000 format(7(1x,f8.4)) !si ecriture dans la console

end subroutine visualisation

!Subroutine pour la visualisation des termes sources
!------------------------------------------------------------------------------
subroutine visualisation_source (S,Nx,Ny,srcname)
  use ISO_FORTRAN_ENV

  integer :: Nx, Ny,j,i
  real*8, dimension (0:Nx,0:Ny) :: S
  character*3 srcname
  integer ::  unit_descriptor

!  unit_descriptor=44          !ecriture dans un fichier fort.44
  unit_descriptor=OUTPUT_UNIT !ecriture dans la console
 
  write(unit_descriptor,*) '==========='
  write(unit_descriptor,*) 'Terme source : ',srcname
  write(unit_descriptor,*) '==========='
  do j=1,Ny-1
  do i=1,Nx-1
     write(unit_descriptor,1000) S(i,j)
  enddo
  enddo

! 1000 format((1x,E16.9)) !si ecriture dans un fichier fort.44
 1000 format(7(1x,f8.4)) !si ecriture dans la console

end subroutine visualisation_source


