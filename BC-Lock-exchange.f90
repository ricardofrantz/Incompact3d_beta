module flow_type
  use decomp_2d, only : mytype
  real(mytype) :: pfront
  integer :: iprocessing

end module flow_type

subroutine ft_parameter(arg)

  USE param
  USE variables
  USE flow_type
  USE complex_geometry
  USE decomp_2d, only : nrank
  implicit none

  logical,intent(in) :: arg
  integer :: is
  character :: a

  iscalar = 1

  open(10,file='BC-Lock-exchange.prm',status='unknown',form='formatted')
  read (10,*) a ! 
  read (10,*) a ! INCOMPACT 3D computational parameters
  read (10,*) a !
  read (10,*) nx
  read (10,*) ny
  read (10,*) nz
  read (10,*) nphi
  read (10,*) p_row
  read (10,*) p_col
  if (arg) then
     close(10)
     return
  endif
  read (10,*) a ! 
  read (10,*) a ! INCOMPACT 3D Flow parameters
  read (10,*) a !
  read (10,*) xlx
  read (10,*) yly
  read (10,*) zlz
  read (10,*) re
  read (10,*) pfront  
  do is=1,nphi
     read (10,*) a
     read (10,*) ri(is)
     read (10,*) nsc(is)
     read (10,*) uset(is)
     read (10,*) cp(is)
  enddo
  read (10,*) a
  read (10,*) angle
  read (10,*) u1
  read (10,*) u2
  read (10,*) noise
  read (10,*) iin
  read (10,*) dt
  read (10,*) a !
  read (10,*) a ! INCOMPACT3D Flow configuration
  read (10,*) a !
  read (10,*) ifirst
  read (10,*) ilast
  read (10,*) nscheme
  read (10,*) cont_phi
  read (10,*) a !velocity
  read (10,*) nclx1
  read (10,*) nclxn
  read (10,*) ncly1
  read (10,*) nclyn
  read (10,*) nclz1
  read (10,*) nclzn
  read (10,*) a !scalar
  read (10,*) nclxS1
  read (10,*) nclxSn
  read (10,*) nclyS1
  read (10,*) nclySn
  read (10,*) nclzS1
  read (10,*) nclzSn
  read (10,*) a !
  read (10,*) a ! INCOMPACT 3D File parameters
  read (10,*) a !
  read (10,*) ilit
  read (10,*) isave
  read (10,*) imodulo
  read (10,*) iprocessing
  read (10,*) a !
  read (10,*) a ! NUMERICAL DISSIPATION
  read (10,*) a
  read (10,*) jLES
  read (10,*) fpi2
  close(10)

  if (nrank==0) then
     print *,'=====================Lock-exchange========================='
     write(*,"(' pfront             : ',F15.8)") pfront
     write(*,"(' iprocessing        : ',I8)") iprocessing
     print *,'==========================================================='
  endif
  return
end subroutine ft_parameter
!********************************************************************
subroutine boundary_conditions (ux1,uy1,uz1,phi1,ep1)

  USE param
  USE variables
  USE decomp_2d
  USE MPI

  implicit none

  integer  :: i,j,k,is
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi2


  do is=1, nphi
     if (uset(is) .eq. 0.) cycle
     call transpose_x_to_y(phi1(:,:,:,is),phi2)

     do k=1,ysize(3)
        do i=1,ysize(1)

           !phi2(i,yend(2),k)= phi2(i,yend(2)-1,k) / (1.+uset(is)*dy*nsc(is)/xnu) !Robin on top BC

           if ( phi2(i,2,k) .gt. phi2(i,1,k) ) then

              phi2(i,1,k)= phi2(i,1,k) + ((uset(is)*gdt(itr))/dy)*(phi2(i,2,k)-phi2(i,1,k)) !Deposit on bottom BC
           else
              phi2(i,1,k)= phi2(i,2,k)! dc/dn=0
           endif

        enddo
     enddo
     call transpose_y_to_x(phi2,phi1(:,:,:,is))
  enddo

  return
end subroutine boundary_conditions
!********************************************************************
subroutine init (ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)

  USE flow_type
  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1,phis1,phiss1

  real(mytype) :: y,r,um,r1,r2,r3,x,z,h
  integer :: k,j,i,ijk,fh,ierror,ii,is,code
  integer (kind=MPI_OFFSET_KIND) :: disp

  if (iscalar==1) then

     !lock-exchange

     do is=1,nphi
        do k=1,xsize(3)
           do j=1,xsize(2)
              do i=1,xsize(1)
                 x=real(i+xstart(1)-1-1,mytype)*dx-pfront
                 phi1(i,j,k,is)=(0.5-0.5*tanh((nsc(is)/xnu)**(0.5)*x))*cp(is)
              enddo
           enddo
        enddo
        do ijk = 1,xsize(1)*xsize(2)*xsize(3)
           if (phi1(ijk,1,1,is).gt.cp(is)) phi1(ijk,1,1,is) =  cp(is)
            if (phi1(ijk,1,1,is).lt.zero) phi1(ijk,1,1,is) =  zero
        enddo
     enddo

     do is=1,nphi
        do ijk=1,xsize(1)*xsize(2)*xsize(3)
           phis1(ijk,1,1,is)=phi1(ijk,1,1,is)
           phiss1(ijk,1,1,is)=phis1(ijk,1,1,is)
        enddo
     enddo

  endif

  ux1=zero; uy1=zero; uz1=zero

  if (iin.ne.0) then
     call system_clock(count=code)
     if (iin.eq.2) code=0
     call random_seed(size = ii)
     call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

     call random_number(ux1)
     call random_number(uy1)
     call random_number(uz1)

     !lock-exchange
     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)        
              x=real(i-1,mytype)*dx-pfront
              um=exp(-twentyfive*x*x)*noise
              ux1(i,j,k)=um*(two*ux1(i,j,k)-one)
              uy1(i,j,k)=um*(two*uy1(i,j,k)-one)
              uz1(i,j,k)=um*(two*uz1(i,j,k)-one)
           enddo
        enddo
     enddo
  endif


  !INIT FOR G AND U=MEAN FLOW + NOISE
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           ux1(i,j,k)=ux1(i,j,k)!+bxx1(j,k)
           uy1(i,j,k)=uy1(i,j,k)!+bxy1(j,k)
           uz1(i,j,k)=uz1(i,j,k)!+bxz1(j,k)
           gx1(i,j,k)=ux1(i,j,k)
           gy1(i,j,k)=uy1(i,j,k)
           gz1(i,j,k)=uz1(i,j,k)
           hx1(i,j,k)=gx1(i,j,k)
           hy1(i,j,k)=gy1(i,j,k)
           hz1(i,j,k)=gz1(i,j,k)
        enddo
     enddo
  enddo

#ifdef DEBG 
  if (nrank .eq. 0) print *,'# init end ok'
#endif

  return
end subroutine init
!*******************************************************************
module post_processing

  USE decomp_2d
  USE variables
  USE param
  USE flow_type

  implicit none
  !
  real(mytype), save, allocatable, dimension(:,:,:) :: vol1
  real(mytype), save, allocatable, dimension(:,:) :: area2
  !
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character
  !
  !probes
  integer :: nprobes
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes

contains

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: dxdydz, dxdz, x, xprobes, yprobes, zprobes
    integer :: i,j,k,code
    character :: a

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post start'
#endif

    call alloc_x(vol1, opt_global=.true.)
    vol1 = zero

    !X PENCILS !Utilizar para integral volumétrica dentro do domínio físico (método de Simpson)
    dxdydz=dx*dy*dz
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          do i=xstart(1),xend(1)
             vol1(i,j,k)=dxdydz
             if (i .eq. 1 .or. i .eq. nx) vol1(i,j,k) = vol1(i,j,k) * five/twelve
             if (j .eq. 1 .or. j .eq. ny) vol1(i,j,k) = vol1(i,j,k) * five/twelve
             if (k .eq. 1 .or. k .eq. nz) vol1(i,j,k) = vol1(i,j,k) * five/twelve
             if (i .eq. 2 .or. i .eq. nx-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
             if (j .eq. 2 .or. j .eq. ny-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
             if (k .eq. 2 .or. k .eq. nz-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
          end do
       end do
    end do

    !Y PENCILS
    allocate(area2(ystart(1):yend(1),ystart(3):yend(3)))
    dxdz=dx*dz
    area2=zero
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
          area2(i,k)=dxdz
          if (i .eq. 1 .or. i .eq. nx) area2(i,k) = area2(i,k)/two
          if (k .eq. 1 .or. k .eq. nz)  area2(i,k) = area2(i,k)/two
       end do
    end do

    !probes
    !WORK X-PENCILS
    open(10,file='probes.prm',status='unknown',form='formatted')
    read (10,*) nprobes
    read (10,*) a
    if (nprobes .gt. 0) then
       allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
       rankprobes(:)=0
       do i=1, nprobes
          read (10,*) xprobes, yprobes, zprobes
          !x
          if (nclx) then
             nxprobes(i)=int(xprobes/dx)
          else
             nxprobes(i)=int(xprobes/dx+1)
          end if
          !y
          if (ncly) then
             nyprobes(i)=int(yprobes/dy)
          else
             nyprobes(i)=int(yprobes/dy+1)
          end if
          !z
          if (nclz) then
             nzprobes(i)=int(zprobes/dz)
          else
             nzprobes(i)=int(zprobes/dz+1)
          end if
          if       (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
             if    (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
                if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
                   rankprobes(i)=1
                endif
             endif
          endif
       enddo
    endif
    close(10)

#ifdef DEBG 
    if (nrank .eq. 0) print *,'# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine postprocessing(ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1

    real(mytype),dimension(ysize(1),ysize(2),ysize(3),nphi) :: phi2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3),nphi) :: phi3
    real(mytype),dimension(ysize(1),ysize(3),nphi) :: dep2
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phisum1
    real(mytype),dimension(zsize(1),zsize(2),nphi) :: phim3
    integer :: i,j,k,is

    real(mytype) :: mp(nphi),dms(nphi),xf(1:2,1:3),xf2d(1:2,1:2)

    if (mod(itime,iprocessing).ne.0) return

    mp=zero; dms=zero; xf=zero; xf2d=zero

    phisum1=zero
    do is=1, nphi
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                phisum1(i,j,k)=phisum1(i,j,k)+phi1(i,j,k,is)
             enddo
          end do
       end do
       call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))
       call transpose_y_to_z(phi2(:,:,:,is),phi3(:,:,:,is))
       call mean_plane_z(phi3(:,:,:,is),zsize(1),zsize(2),zsize(3),phim3(:,:,is))
    enddo

    call budget(ux1,uy1,uz1,phi1,vol1)
    call dep(phi1,dep2)
    call suspended(phi1,vol1,mp)
    call depositrate (dep2,dms)
    call front(phisum1,xf)
    call front2d(phim3,xf2d)

    if (nrank .eq. 0) then
       FS = 1+nphi+nphi+3+2 !Number of columns
       write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
       FS = FS*14+1  !Line width
       open(67,file='./out/statistics',status='unknown',form='formatted',&
            access='direct',recl=FS)
       write(67,fileformat,rec=itime/iprocessing+1) t,& !1
            mp,&                                    !nphi
            dms,&                                   !nphi
            xf(1,:),&                               !3
            xf2d(1,:),&                             !2
            NL                                      !+1
       close(67)
    end if

  end subroutine postprocessing
  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,phi1) !By Felipe Schuch

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),nphi) :: phi1

    integer :: i
    character(len=30) :: filename
    FS = 1+3+nphi !Number of columns
    write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
    FS = FS*14+1  !Line width

    do i=1, nprobes
       if (rankprobes(i) .eq. 1) then
          write(filename,"('./out/probe',I4.4)") i
          open(67,file=trim(filename),status='unknown',form='formatted'&
               ,access='direct',recl=FS)
          write(67,fileformat,rec=itime) t,&                         !1
               ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
               uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
               uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
               phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !nphi
               NL                                                    !+1
          close(67)
       endif
    enddo

  end subroutine write_probes
  !############################################################################
  subroutine suspended(phi1,vol1,mp1)

    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: vol1
    real(mytype),intent(out) :: mp1(1:nphi)

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1
    real(mytype) :: mp(1:nphi)
    integer :: is,code

    mp=zero; mp1=zero

    do is=1, nphi
       temp1 = phi1(:,:,:,is)*vol1(:,:,:)
       mp(is)= sum(temp1)
    end do

    call MPI_REDUCE(mp,mp1,nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    return
  end subroutine suspended
  !############################################################################
  subroutine front ( phisum1, xp )

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: phisum1
    real(mytype),intent(out) :: xp(1:2,1:3)

    real(mytype) :: xp1(1:2),y
    integer :: is, i, j ,k, code

    xp(2,:) = real(nrank,mytype)
    xp(1,:)=zero
    xp1=zero
    kloop: do k=xstart(3),xend(3) 
       jloop: do j=xstart(2),xend(2)
          iloop: do i=xend(1), xstart(1), -1
             if ( phisum1(i,j,k) .ge. 0.01_mytype) then
                xp(1,1:3) = (/ real (i-1,mytype)*dx, real (j-1,mytype)*dy, real(k-1,mytype)*dz /)
                exit kloop
             end if
          end do iloop
       end do jloop
    end do kloop

    call MPI_ALLREDUCE(xp(:,1),xp1,1,real2_type,MPI_MAXLOC,MPI_COMM_WORLD,code)
    call MPI_Bcast(xp(1,:), 3,real_type, int(xp1(2)), MPI_COMM_WORLD,code)


  end subroutine front
  !############################################################################
  subroutine front2d ( phim3, xp )

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(zstart(1):zend(1),zstart(2):zend(2),1:nphi) :: phim3
    real(mytype),intent(out) :: xp(1:2,1:2)
    real(mytype) :: xp1(1:2),y
    integer :: is, i, j, code

    xp(2,:) = real(nrank,mytype)
    xp(1,:) = zero
    xp1=zero
    jloop: do j=zstart(2),zend(2)
       y = real( j - 1, mytype )*dy
       iloop: do i=zend(1), zstart(1), -1
          if ( sum(phim3(i,j,:)) .ge. 0.01_mytype) then
             xp(1,1:2) = (/ real (i-1,mytype)*dx, real (j-1,mytype)*dy /)
             exit jloop
          end if
       end do iloop
    end do jloop

    call MPI_ALLREDUCE(xp(:,1),xp1,1,real2_type,MPI_MAXLOC,MPI_COMM_WORLD,code)
    call MPI_Bcast(xp(1,:), 2,real_type, int(xp1(2)), MPI_COMM_WORLD,code)


  end subroutine front2d
  !############################################################################
  subroutine dep(phi1,dep2)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),nphi) :: phi1
    real(mytype),intent(out),dimension(ystart(1):yend(1),ystart(3):yend(3),nphi) :: dep2

    real(mytype),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),nphi) :: tempdep2
    real(mytype),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),nphi) :: phi2

    integer :: i, j, k, is
    character(len=30) :: filename

    dep2 = zero
    do is=1, nphi
       if (uset(is) .eq. 0.) cycle
       call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))

       tempdep2=zero

       do k=ystart(3),yend(3) 
          do i=ystart(1),yend(1)
             tempdep2(i,1,k,is) = phi2(i,1,k,is)*uset(is)
             dep2(i,k,is) = tempdep2(i,1,k,is)
             !dep2(i,k,is) = dep2(i,k,is) +  phi2(i,1,k,is)*uset(is)
          end do
       end do

       write(filename,"('./out/dep',I1.1,I4.4)") is,itime/iprocessing
       call decomp_2d_write_plane(2,tempdep2(:,:,:,is),2,1,filename)
    enddo

  end subroutine dep
  !############################################################################
  subroutine depositrate ( dep2, dms1)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(ystart(1):yend(1),ystart(3):yend(3),nphi) :: dep2
    real(mytype),intent(out) :: dms1(nphi)

    real(mytype) :: dms(nphi)
    integer :: i,k,is,code

    dms=zero; dms1=zero
    do is=1, nphi
       if (uset(is) .eq. 0.) cycle
       do k=ystart(3),yend(3)
          do i=ystart(1),yend(1)
             dms(is)=dms(is)+dep2(i,k,is)*area2(i,k)
          end do
       end do
    enddo

    call MPI_REDUCE(dms,dms1,nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

  end subroutine depositrate
  !############################################################################
  subroutine budget(ux1,uy1,uz1,phi1,vol1)

    USE decomp_2d
    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,vol1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: diss1,dphiy1, dphixx1, dphiyy1, dphizz1, ddphi1, di1, temp1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi2, dphiy2, dphixx2, dphiyy2, dphizz2, ddphi2, di2, vol2, temp2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: phi3, dphiy3, dphixx3, dphiyy3, dphizz3, ddphi3, di3, temp3

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3
    real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A

    real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu
  
    real(8) :: ek,ek1,dek,dek1,ep,ep1,dep,dep1,xvol
    integer :: ijk,i,j,k,l,m,is,code
    character(len=30) :: filename

    ek=zero;ek1=zero;dek=zero;dek1=zero;ep=zero;ep1=zero;dep=zero;dep1=zero;diss1=zero

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    !y-derivatives
    call transpose_x_to_y(ux1,td2)
    call transpose_x_to_y(uy1,te2)
    call transpose_x_to_y(uz1,tf2)
    call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    !!z-derivatives
    call transpose_y_to_z(td2,td3)
    call transpose_y_to_z(te2,te3)
    call transpose_y_to_z(tf2,tf3)
    call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)

    A(:,:,:,:,:)=zero
    A(1,1,:,:,:)=ta1(:,:,:)
    A(2,1,:,:,:)=tb1(:,:,:)
    A(3,1,:,:,:)=tc1(:,:,:)
    A(1,2,:,:,:)=td1(:,:,:)
    A(2,2,:,:,:)=te1(:,:,:)
    A(3,2,:,:,:)=tf1(:,:,:)
    A(1,3,:,:,:)=tg1(:,:,:)
    A(2,3,:,:,:)=th1(:,:,:)
    A(3,3,:,:,:)=ti1(:,:,:)

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             do m=1,3
                do l=1,3
                   diss1(i,j,k)=diss1(i,j,k)+two*xnu*half*half*(A(l,m,i,j,k)+A(m,l,i,j,k))**two
                enddo
             enddo
          enddo
       enddo
    enddo

    do ijk=1,xsize(1)*xsize(2)*xsize(3)
       xvol=real(vol1(ijk,1,1),8)
       ek = ek + half * xvol * (ux1(ijk,1,1)**two+uy1(ijk,1,1)**two+uz1(ijk,1,1)**two)
       dek = dek + xvol * diss1(ijk,1,1)
    enddo

    call transpose_x_to_y(vol1,vol2)

    if (ivirt==2) then
       ilag=0
    endif
    do is=1, nphi
       if (ri(is) .eq. 0.) cycle
       call derxxS (dphixx1,phi1(:,:,:,is),di1,sx,sfxpS,ssxpS,swxpS,xsize(1),xsize(2),xsize(3),1)

       call transpose_x_to_y(dphixx1,dphixx2)

       call transpose_x_to_y(phi1(:,:,:,is),phi2)

       call deryS (dphiy2,phi2,di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),1)

       call deryyS (dphiyy2,phi2,di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1)

       call transpose_y_to_z(phi2,phi3)

       call derzzS (dphizz3,phi3,di3,sz,sfzpS,sszpS,swzpS,zsize(1),zsize(2),zsize(3),1)

       call transpose_z_to_y(dphizz3,dphizz2)

       do ijk=1,ysize(1)*ysize(2)*ysize(3)
          ddphi2(ijk,1,1)=dphixx2(ijk,1,1)+dphiyy2(ijk,1,1)+dphizz2(ijk,1,1)
       enddo

       do k=1,ysize(3)
          do j=1,ysize(2)
             do i=1,ysize(1)
                xvol=real(vol2(i,j,k),8)
                ep=ep + xvol * (phi2(i,j,k)*(j-1)*dy)
                dep=dep - xvol * (ddphi2(i,j,k)*xnu/nsc(is)+uset(is)*dphiy2(i,j,k))*(j-1)*dy
             enddo
          enddo
       enddo
    enddo
    if (ivirt==2) then
       ilag=1
    endif

    call MPI_REDUCE(ek,ek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(dek,dek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(ep,ep1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(dep,dep1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)

    if (nrank .eq. 0) then
       open(67,file='./out/budget',status='unknown',form='formatted',&
            access='direct',recl=71) !71=5*14+1
       write(67,"(5E14.6,A)",rec=itime/iprocessing+1) t,ek1,dek1,ep1,dep1,NL
       close(67)
    end if

    if (mod(itime,imodulo).eq.0) then
     !if (save_diss.eq.1) then
        uvisu=zero
        call fine_to_coarseV(1,diss1,uvisu)
        write(filename,"('./data/diss',I4.4)") itime/imodulo
        call decomp_2d_write_one(1,uvisu,filename,2)
     !endif

     !if (save_dissm.eq.1) then
        call transpose_x_to_y (diss1,temp2)
        call transpose_y_to_z (temp2,temp3)
        call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
        write(filename,"('./data/dissm',I4.4)") itime/imodulo
        call decomp_2d_write_plane(3,temp3,3,1,filename)
     !endif
    endif

  end subroutine budget
end module post_processing
