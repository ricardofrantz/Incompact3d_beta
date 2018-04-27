module flow_type
  use decomp_2d, only : mytype

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

  open(10,file='BC-ekman.prm',status='old',form='formatted')
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
  read (10,*) ro
  read (10,*) noise
  read (10,*) dt
  do is=1,nphi
    read (10,*) a
    read (10,*) ri(is)
    read (10,*) nsc(is)
    read (10,*) uset(is)
    read (10,*) cp(is)
  enddo
  read (10,*) a !
  read (10,*) a ! INCOMPACT3D Flow configuration
  read (10,*) a !
  read (10,*) iin
  read (10,*) ifirst
  read (10,*) ilast
  read (10,*) nscheme
  read (10,*) istret
  read (10,*) beta
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
  read (10,*) itest
  read (10,*) a !
  read (10,*) a ! NUMERICAL DISSIPATION
  read (10,*) a !
  read (10,*) jLES
  read (10,*) fpi2

  if (nrank==0) then
    print *,'==================Ekman flow==================='
    print *,'==============================================='
  endif
  return
end subroutine ft_parameter
!********************************************************************
subroutine init (ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uz2

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1,phis1,phiss1

  real(mytype) :: y,r,r3,x,z,h,ct
  real(mytype) :: cx0,cy0,cz0,hg,lg
  real(mytype) :: a_a,a_0,a_1,b_0,b_1,b_2,b_3,c_c,c_0,c_1,c_2,c_3,d_0,d_1,d_2,d_3

  integer :: k,j,i,ijk,fh,ierror,ii,is,code
  integer (kind=MPI_OFFSET_KIND) :: disp
  real(mytype), dimension(ysize(2)) :: um
  integer, dimension (:), allocatable :: seed
  real(mytype), save, allocatable, dimension(:) :: uansorge, vansorge
  allocate(uansorge(ysize(2)),vansorge(ysize(2)))

  a_a =  1.94150000
  a_0 = -1.65020716
  a_1 =  0.26934755
  b_0 =  1.36058818e5
  b_1 = -0.212789514
  b_2 =  59.14127000
  b_3 =  0.122560653

  c_c =  1.94150000
  c_0 =  0.69517727
  c_1 =  2.03288356
  c_2 = -1.35471899
  c_3 =  0.18865185
  d_0 = -2.55958895
  d_1 = -0.26217732
  d_2 =  6.09783536
  d_3 = -0.08174588

  do j=1,ysize(2)
    if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
    if (istret.ne.0) y=yp(j+xstart(2)-1)
    if (y.lt.c_c) then
      um(j) = c_0*y**(c_1-one)*exp(c_2*y**(c_1-one) + c_3*y**(c_1))
    else
      um(j) = d_0*exp(d_1*(y+d_2))*sin(d_3*(y+d_2))
    endif
  enddo

  um = um / maxval(um)


  if (iscalar==1) then

    do is=1,nphi

      call random_number(phi1(:,:,:,is))

      do k=1,xsize(3)
        do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)
          do i=1,xsize(1)
            phi1(i,j,k,is) = um(j)*noise*(2.*phi1(i,j,k,is)-1.) + erf(y*uset(is))
          enddo
        enddo
      enddo

    enddo

    !do not delete this
    phis1=phi1
    phiss1=phis1

  endif

  ux1=zero;uy1=zero;uz1=zero

  if (iin==1) then

    call system_clock(count=code)
    if (iin.eq.2) code=0
    call random_seed(size = ii)
    call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

    call random_number(ux1)
    call random_number(uy1)
    call random_number(uz1)

    if (ncly1==2) then !ensure no noise near the wall !!! done in the hard way
      if (xstart(2)==1) then
        do k=1,xsize(3)
          do i=1,xsize(1)
            ux1(i,1,k)=zero
            uy1(i,1,k)=zero
            uz1(i,1,k)=zero
          enddo
        enddo
      endif
    endif

    !espiral de Ekman
    !do k=1,xsize(3)
    ! do j=1,xsize(2)
    !    if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
    !    if (istret.ne.0) y=yp(j+xstart(2)-1)
    !       do i=1,xsize(1)  !pfront aqui é deltaE
    !       !Ug = 1. e consideramos a coordenada vertical como y
    !       ux1(i,j,k) = ux1(i,j,k) + one*(one-exp(-y/delta_Ekman)*cos(y/delta_Ekman))
    !       uz1(i,j,k) = uz1(i,j,k) + one*exp(-y/delta_Ekman)*sin(y/delta_Ekman)
    !  enddo
    ! enddo
    !enddo

    !ansorge com fits
    do k=1,xsize(3)
      do j=1,xsize(2)
        if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
        if (istret.ne.0) y=yp(j+xstart(2)-1)
        do i=1,xsize(1)

          if (y.lt.a_a) then
            ux1(i,j,k) = noise*um(j)*(two*ux1(i,j,k)-one) + one-exp(a_0*y + a_1*y**2)
          else
            ux1(i,j,k) = noise*um(j)*(two*ux1(i,j,k)-one) + one-b_0*exp(b_1*(y+b_2))*cos(b_3*(y+b_2))
          endif

          uy1(i,j,k)=noise*um(j)*(two*uy1(i,j,k)-one)

          if (y.lt.c_c) then
            uz1(i,j,k) = noise*um(j)*(two*uz1(i,j,k)-one) + c_0*y**(c_1-one)*exp(c_2*y**(c_1-one) + c_3*y**(c_1))
          else
            uz1(i,j,k) = noise*um(j)*(two*uz1(i,j,k)-one) + d_0*exp(d_1*(y+d_2))*sin(d_3*(y+d_2))
          endif

          !if (y.le. 2.16                  ) then
          !  ux1(i,j,k) = ux1(i,j,k) + 0.02144   + 1.15392*y   - 0.79181*y**3   + 0.46653*y**4   - 0.08034*y**5
          !elseif (y.gt. 2.16 .AND. y.le. 16.3  ) then
          !  ux1(i,j,k) = ux1(i,j,k) + 8.051e-01 + 5.072e-02*y - 7.359e-04*y**3 + 6.008e-05*y**4 - 1.453e-06*y**5
          !elseif (y.gt. 16.3 .AND. y.le. 28.24 ) then
          !  ux1(i,j,k) = ux1(i,j,k) + 1.098e+00 - 6.238e-03*y + 3.791e-06*y**3 + 7.596e-08*y**4 - 3.060e-09*y**5
          !elseif (y.gt. 28.24                 ) then
          !  ux1(i,j,k) = ux1(i,j,k) + 9.945e-01 + 2.481e-04*y - 2.923e-09*y**3 - 2.695e-09*y**4 + 3.375e-11*y**5
          !endif

          !if (y.le. 2.4                  ) then
          !  uz1(i,j,k) = uz1(i,j,k) + -0.0001992 + 0.6344476*y - 0.6641268*y**2 + 0.2702919*y**3 - 0.0305002*y**5 + 0.0061358*y**6
          !elseif (y.gt. 2.4 .AND. y.le. 16.3 ) then
          !  uz1(i,j,k) = uz1(i,j,k) + 2.679e-01  - 4.611e-02*y + 3.672e-03*y**2 - 1.654e-04*y**3 + 4.148e-07*y**5 - 1.326e-08*y**6
          !elseif (y.gt. 16.3                 ) then
          !  uz1(i,j,k) = uz1(i,j,k) + 8.523e-02  - 1.210e-02*y + 6.228e-04*y**2 - 1.208e-05*y**3 + 2.617e-09*y**5 - 2.208e-11*y**6
          !endif

        enddo
      enddo
    enddo

  endif

  !  if (iin==2) then !import fields as initial condition
  !    if (nrank==0) print *,'Reading external files"'


  !    call system_clock(count=code)
  !    if (iin.eq.2) code=0
  !    call random_seed(size = ii)
  !    call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

  !    call random_number(ux1)
  !    call random_number(uy1)
  !    call random_number(uz1)

  !    !modulation of the random noise + initial velocity profile
  !    do k=1,xsize(3)
  !      do j=1,xsize(2)
  !        if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy- 0.8
  !        if (istret.ne.0) y=yp(j+xstart(2)-1)- 0.8
  !        um=exp(-zptwo*y*y)
  !        do i=1,xsize(1)
  !          ux1(i,j,k)=noise*um(j)*(two*ux1(i,j,k)-one)
  !          uy1(i,j,k)=noise*um(j)*(two*uy1(i,j,k)-one)
  !          uz1(i,j,k)=noise*um(j)*(two*uz1(i,j,k)-one)
  !        enddo
  !      enddo
  !    enddo

  !    if (ncly1==2) then !ensure no noise near the wall !!! done in the hard way
  !      if (xstart(2)==1) then
  !        do k=1,xsize(3)
  !          do i=1,xsize(1)
  !            ux1(i,1,k)=zero
  !            uy1(i,1,k)=zero
  !            uz1(i,1,k)=zero
  !          enddo
  !        enddo
  !      endif
  !    endif

  !    call transpose_x_to_y(ux1,ux2)
  !    call transpose_x_to_y(uz1,uz2)

  !    open(10,file='ansorge.txt',status='unknown',form='formatted')
  !    do i=1,ysize(2)
  !      read(10,*) uansorge(i) , vansorge(i)
  !      !read(10,*) ux2(:,i,:), uz2(:,i,:)
  !    enddo

  !    do i=1,ysize(2)
  !      ux2(:,i,:) = ux2(:,i,:) + uansorge(i)
  !      uz2(:,i,:) = uz2(:,i,:) - vansorge(i)
  !    enddo

  !    call transpose_y_to_x(ux2,ux1)
  !    call transpose_y_to_x(uz2,uz1)

  !    !call decomp_2d_read_one(1,ux1,'ux1_init.dat')
  !    !call decomp_2d_read_one(1,uy1,'uy1_init.dat')
  !    !call decomp_2d_read_one(1,uz1,'uz1_init.dat')

  !  endif

  if (ncly1==2) then !ensure no noise near the wall !!! done in the hard way
    if (xstart(2)==1) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          ux1(i,1,k)=zero
          uy1(i,1,k)=zero
          uz1(i,1,k)=zero
        enddo
      enddo
    endif
  endif

  !INIT FOR G AND U=MEAN FLOW + NOISE
  do k=1,xsize(3)
    do j=1,xsize(2)
      do i=1,xsize(1)
        ux1(i,j,k)=ux1(i,j,k)
        uy1(i,j,k)=uy1(i,j,k)
        uz1(i,j,k)=uz1(i,j,k)
        gx1(i,j,k)=ux1(i,j,k)
        gy1(i,j,k)=uy1(i,j,k)
        gz1(i,j,k)=uz1(i,j,k)
        hx1(i,j,k)=gx1(i,j,k)
        hy1(i,j,k)=gy1(i,j,k)
        hz1(i,j,k)=gz1(i,j,k)
      enddo
    enddo
  enddo

  return
end subroutine init
!********************************************************************
subroutine boundary_conditions (ux1,uy1,uz1,phi1,ep1)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  integer  :: i,j,k,is
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi2

  byx1=zero
  byy1=zero
  byz1=zero

  !  do is=1, nphi
  !     call transpose_x_to_y(phi1(:,:,:,is),phi2)

  !     do k=1,ysize(3)
  !        do i=1,ysize(1)

  !           !phi2(i,yend(2),k)= phi2(i,yend(2)-1,k) / (1.+uset(is)*dy*nsc(is)/xnu) !Robin on top BC

  !           if ( phi2(i,2,k) .gt. phi2(i,1,k) ) then
  !              phi2(i,1,k)= phi2(i,1,k) + ((uset(is)*gdt(itr))/dy)*(phi2(i,2,k)-phi2(i,1,k)) !Deposit on bottom BC
  !           else
  !              phi2(i,1,k)= phi2(i,2,k)! dc/dn=0
  !           endif

  !        enddo
  !     enddo
  !     call transpose_y_to_x(phi2,phi1(:,:,:,is))
  !  enddo

  if (nclyS1.eq.2) then
    do is=1,nphi
      if (xstart(2)==1) then
        do k=1,xsize(3)
          do i=1,xsize(1)
            !if (nrank==0) write(*,*) 'antes',phi1(i,1,k,is)
            phi1(i,1,k,is)=zero
            !if (nrank==0) write(*,*) 'depois',phi1(i,1,k,is)
          enddo
        enddo
      endif
    enddo
  endif

  if (nclySn.eq.2) then
    do is=1,nphi
      if (xend(2)==ny) then
        do k=1,xsize(3)
          do i=1,xsize(1)
            phi1(i,xsize(2),k,is)=cp(is)
          enddo
        enddo
      endif
    enddo
  endif

  return

end subroutine boundary_conditions
#ifdef POST
!********************************************************************
module post_processing

  USE decomp_2d
  USE variables
  USE param
  USE flow_type

  implicit none
  !
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character
  !
  !probes
  integer, save :: nprobes, ntimes1, ntimes2
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes
  !
  !real(mytype),save,allocatable,dimension(:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

  real(mytype), save, allocatable, dimension(:,:,:) :: vol1
  real(mytype), save, allocatable, dimension(:,:) :: area2

contains

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: x, xprobes, yprobes, zprobes, dxdydz, dxdz
    integer :: i,j,k,code
    character :: a

    !allocate(usum(ysize(2)),vsum(ysize(2)),wsum(ysize(2)))
    !allocate(uusum(ysize(2)),uvsum(ysize(2)),uwsum(ysize(2)))
    !allocate(vvsum(ysize(2)),vwsum(ysize(2)),wwsum(ysize(2)))

    call alloc_x(vol1, opt_global=.true.)

    !usum=zero;vsum=zero;wsum=zero
    !uusum=zero;uvsum=zero;uwsum=zero
    !vvsum=zero;vwsum=zero;wwsum=zero
    !ntimes1 = 0
    !ntimes2 = 0
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

    !    nprobes  = 0
    !    !probes
    !    !WORK X-PENCILS
    !    open(10,file='probes.prm',status='unknown',form='formatted')
    !    read (10,*) nprobes
    !    read (10,*) a
    !    if (nprobes .gt. 0) then
    !       allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
    !       rankprobes(:)=0
    !       do i=1, nprobes
    !          read (10,*) xprobes, yprobes, zprobes
    !          !x
    !          if (nclx) then
    !             nxprobes(i)=int(xprobes/dx)
    !          else
    !             nxprobes(i)=int(xprobes/dx+1)
    !          end if
    !          !y
    !          if (ncly) then
    !             nyprobes(i)=int(yprobes/dy)
    !          else
    !             nyprobes(i)=int(yprobes/dy+1)
    !          end if
    !          !z
    !          if (nclz) then
    !             nzprobes(i)=int(zprobes/dz)
    !          else
    !             nzprobes(i)=int(zprobes/dz+1)
    !          end if
    !          if       (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
    !             if    (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
    !                if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
    !                   rankprobes(i)=1
    !                endif
    !             endif
    !          endif
    !       enddo
    !    endif
    !    close(10)

    nprobes = ny - 2
    allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
    rankprobes(:)=0

    do i=1, nprobes-1
      nyprobes(i) = i+1
    enddo
    nxprobes(:) = nx / 2
    nzprobes(:) = nz / 2

    do i=1, nprobes
      if       (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
        if    (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
          if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
            rankprobes(i)=1
          endif
        endif
      endif
    enddo

  end subroutine init_post
  !############################################################################
  subroutine postprocessing(ux1,uy1,uz1,phi1,ep1,pre1,diss1)

    USE param
    USE variables
    USE decomp_2d
    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, pre1, ep1
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3)) :: diss1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2 , uy2,  uz2,  phi2
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1p, uy1p, uz1p, phi1p, pre1p
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2p, uy2p, uz2p, phi2p, pre2p
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,temp1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,di2,temp2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,temp3
    real(mytype),dimension(ysize(2)) :: sumxz, sumxz1
    real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu
    real(mytype) :: mp(nphi)
    real(8) :: tstart
    integer :: i,j,k,is,code
    character(len=1),parameter :: NL=char(10)
    character(len=60) :: filename
    character(len=300) :: fileformat

    tstart=MPI_WTIME()

    if (mod(itime,iprocessing).eq.0) then

      mp=zero

      call budget(ux1,uy1,uz1,phi1,vol1)
      call suspended(phi1,vol1,mp)

      if (nrank .eq. 0) then
        FS = 1+nphi !Number of columns
        write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
        FS = FS*14+1  !Line width
        open(67,file='./out/statistics',status='unknown',form='formatted',&
        access='direct',recl=FS)
        write(67,fileformat,rec=itime/iprocessing+1) t,& !1
        mp,&                                    !nphi
        NL                                      !+1
        close(67)
      end if

    endif

    ta1=zero;tb1=zero;tc1=zero;td1=zero;te1=zero;tf1=zero;tg1=zero;th1=zero;ti1=zero
    ta2=zero;tb2=zero;tc2=zero;td2=zero;te2=zero;tf2=zero;tg2=zero;th2=zero;ti2=zero
    ta3=zero;tb3=zero;tc3=zero;td3=zero;te3=zero;tf3=zero;tg3=zero;th3=zero;ti3=zero

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

    call dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1) !sij sij

    call transpose_x_to_y(diss1,temp2)
    sumxz=zero; sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sum(temp2(:,j,:))
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/diss')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    !!!!!!!!!!!!!!!!!!!!!! u* or ufriction

    call transpose_x_to_y(sqrt(xnu*sqrt(td1**2+ tf1**2)),temp2)
    sumxz=zero; sumxz1=zero
    sumxz(1) = sum(temp2(:,1,:))
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 2
    write(filename,"('./out/ustar')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=(2*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(1),NL
    close(67)


    !!!!!!!!!!!!!!!!!!!!!! ux

    call transpose_x_to_y(ux1,ux2)
    sumxz=zero; sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sum(ux2(:,j,:))
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/ux')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    !!!!!!!!!!!!!!!!!!!!!! ux'

    do k=1,ysize(3)
      do j=1,ysize(2)
        do i=1,ysize(1)
          ux2p(i,j,k) = ux2(i,j,k) - sumxz1(j)
        enddo
      enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!! ux' ux'

    sumxz=zero; sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sum(ux2p(:,j,:)**2)
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/uxp2')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    !!!!!!!!!!!!!!!!!!!!!! uy

    call transpose_x_to_y(uy1,uy2)
    sumxz=zero; sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sum(uy2(:,j,:))
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/uy')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    !!!!!!!!!!!!!!!!!!!!!! uy'

    do k=1,ysize(3)
      do j=1,ysize(2)
        do i=1,ysize(1)
          uy2p(i,j,k) = uy2(i,j,k) - sumxz1(j)
        enddo
      enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!! uy' uy'

    sumxz=zero; sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sum(uy2p(:,j,:)**2)
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/uyp2')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    !!!!!!!!!!!!!!!!!!!!!! uz

    call transpose_x_to_y(uz1,uz2)
    sumxz=zero; sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sum(uz2(:,j,:))
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/uz')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    !!!!!!!!!!!!!!!!!!!!!! uz'

    do k=1,ysize(3)
      do j=1,ysize(2)
        do i=1,ysize(1)
          uz2p(i,j,k) = uz2(i,j,k) - sumxz1(j)
        enddo
      enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!! uz' uz'

    sumxz=zero; sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sum(uz2p(:,j,:)**2)
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/uzp2')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    !!!!!!!!!!!!!!!!!!!!!! tke

    sumxz=zero; sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sum(ux2p(:,j,:)**2 + uy2p(:,j,:)**2 + uz2p(:,j,:)**2)*half
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/tke')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    if (mod(itime,imodulo).eq.0) then
      uvisu=zero
      call fine_to_coarseV(1,ux1p**2+uy1p**2+uz1p**2,uvisu)
      write(filename,"('./data/tke',I4.4)") itime/imodulo
      call decomp_2d_write_one(1,uvisu,filename,2)
    endif

    !!!!!!!!!!!!!!!!!!!!!! turbulent intrnesity

    sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sqrt(sumxz(j))/sum(sqrt(ux2(:,j,:)**2 + uy2(:,j,:)**2 + uz2(:,j,:)**2))
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/tint')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    do is = 1, nphi

      !!!!!!!!!!!!!!!!!!!!!! phi
      call transpose_x_to_y(phi1(:,:,:,is),phi2)
      sumxz=zero; sumxz1=zero
      do j=1,ysize(2)
        sumxz(j) = sum(phi2(:,j,:))
      enddo
      call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      sumxz1=sumxz1/real(nx*nz,mytype)
      write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
      write(filename,"('./out/phi',I1.1)") is
      open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
      write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
      close(67)

      !!!!!!!!!!!!!!!!!!!!!! phi'

      do k=1,ysize(3)
        do j=1,ysize(2)
          do i=1,ysize(1)
            phi2p(i,j,k) = phi2(i,j,k) - sumxz1(j)
          enddo
        enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!! phi' phi'

      sumxz=zero; sumxz1=zero
      do j=1,ysize(2)
        sumxz(j) = sum(phi2p(:,j,:)**2)*half*ri(is)
      enddo
      call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      sumxz1=sumxz1/real(nx*nz,mytype)
      write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
      write(filename,"('./out/phi',I1.1,'p')") is
      open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
      write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
      close(67)

      !!!!!!!!!!!!!!!!!!!!!! phi' u' or buoyancy flux

      sumxz=zero; sumxz1=zero
      do j=1,ysize(2)
        sumxz(j) = sum(phi2p(:,j,:)*ux2p(:,j,:))*ri(is)
      enddo
      call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      sumxz1=sumxz1/real(nx*nz,mytype)
      write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
      write(filename,"('./out/phi',I1.1,'p_uxp')") is
      open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
      write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
      close(67)

      !!!!!!!!!!!!!!!!!!!!!! phi' v' or buoyancy flux

      sumxz=zero; sumxz1=zero
      do j=1,ysize(2)
        sumxz(j) = sum(phi2p(:,j,:)*uy2p(:,j,:))*ri(is)
      enddo
      call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      sumxz1=sumxz1/real(nx*nz,mytype)
      write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
      write(filename,"('./out/phi',I1.1,'p_uyp')") is
      open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
      write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
      close(67)

      !!!!!!!!!!!!!!!!!!!!!! phi' w' or buoyancy flux

      sumxz=zero; sumxz1=zero
      do j=1,ysize(2)
        sumxz(j) = sum(phi2p(:,j,:)*uz2p(:,j,:))*ri(is)
      enddo
      call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      sumxz1=sumxz1/real(nx*nz,mytype)
      write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
      write(filename,"('./out/phi',I1.1,'p_uzp')") is
      open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
      write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
      close(67)

    enddo

    ta1=zero;tb1=zero;tc1=zero;td1=zero;te1=zero;tf1=zero;tg1=zero;th1=zero;ti1=zero
    ta2=zero;tb2=zero;tc2=zero;td2=zero;te2=zero;tf2=zero;tg2=zero;th2=zero;ti2=zero
    ta3=zero;tb3=zero;tc3=zero;td3=zero;te3=zero;tf3=zero;tg3=zero;th3=zero;ti3=zero

    call derx (ta1,ux1p,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tb1,uy1p,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (tc1,uz1p,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    !y-derivatives
    call transpose_x_to_y(ux1p,td2)
    call transpose_x_to_y(uy1p,te2)
    call transpose_x_to_y(uz1p,tf2)
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

    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

    !!!!!!!!!!!!!!!!!!!!!! diss sij sij

    call dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1) !sij sij

    call transpose_x_to_y(diss1,temp2)
    sumxz=zero; sumxz1=zero
    do j=1,ysize(2)
      sumxz(j) = sum(temp2(:,j,:))
    enddo
    call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    sumxz1=sumxz1/real(nx*nz,mytype)
    write(fileformat, '( "(",I4,"(E16.8),A)" )' ) 1+ny
    write(filename,"('./out/dissp')")
    open(67,file=trim(filename),status='unknown',form='formatted',access='direct',recl=((ny+1)*16+1))
    write(67,fileformat,rec=itime+1) t,sumxz1(:),NL
    close(67)

    if (mod(itime,imodulo).eq.0) then
      uvisu=zero
      call fine_to_coarseV(1,diss1,uvisu)
      write(filename,"('./data/dissp',I4.4)") itime/imodulo
      call decomp_2d_write_one(1,uvisu,filename,2)
    endif

    if (nrank.eq.0) print *,'Time in statistics (s)', real(MPI_WTIME()-tstart,4)

    return

  end subroutine postprocessing
  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,pre1,diss1,phi1) !By Felipe Schuch

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1, pre1, diss1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),nphi) :: phi1

    integer :: i
    character(len=30) :: filename
    FS = 1+3+2+nphi !Number of columns
    write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
    FS = FS*14+1  !Line width

    do i=1, nprobes
      if (rankprobes(i) .eq. 1) then
        write(filename,"('./probes/probe',I4.4)") i+2
        open(67,file=trim(filename),status='unknown',form='formatted'&
        ,access='direct',recl=FS)
        write(67,fileformat,rec=itime) t,&                         !1
        ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
        uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
        uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
        pre1(nxprobes(i),nyprobes(i),nzprobes(i)),&           !5
        diss1(nxprobes(i),nyprobes(i),nzprobes(i)),&          !6
        phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !nphi
        NL                                                    !+1
        close(67)
      endif
    enddo

  end subroutine write_probes
  !############################################################################


  subroutine dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1
    real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A
    integer :: k,j,i,m,l

    !INSTANTANEOUS DISSIPATION RATE
    diss1=0._mytype
    A(:,:,:,:,:)=0._mytype
    A(1,1,:,:,:)=ta1(:,:,:)!du/dx=ta1
    A(2,1,:,:,:)=tb1(:,:,:)!dv/dx=tb1
    A(3,1,:,:,:)=tc1(:,:,:)!dw/dx=tc1
    A(1,2,:,:,:)=td1(:,:,:)!du/dy=td1
    A(2,2,:,:,:)=te1(:,:,:)!dv/dy=te1
    A(3,2,:,:,:)=tf1(:,:,:)!dw/dy=tf1
    A(1,3,:,:,:)=tg1(:,:,:)!du/dz=tg1
    A(2,3,:,:,:)=th1(:,:,:)!dv/dz=th1
    A(3,3,:,:,:)=ti1(:,:,:)!dw/dz=ti1
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          do m=1,3
            do l=1,3
              diss1(i,j,k)=diss1(i,j,k)+two*xnu*half*half*(A(l,m,i,j,k)+A(m,l,i,j,k))**2
            enddo
          enddo
        enddo
      enddo
    enddo

    return

  end subroutine dissipation
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

  end subroutine budget

end module post_processing
#endif
