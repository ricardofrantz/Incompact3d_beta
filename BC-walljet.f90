module flow_type
  use decomp_2d, only : mytype
  integer :: initstats1,initstats2

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

  open(10,file='BC-walljet.prm',status='unknown',form='formatted')
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
  read (10,*) iscalar
  read (10,*) cont_phi
  do is=1,nphi
    read (10,*) a
    read (10,*) ri(is)
    read (10,*) nsc(is)
    read (10,*) uset(is)
    read (10,*) cp(is)
  enddo
  read (10,*) noise1
  read (10,*) dt
  read (10,*) u1
  read (10,*) u2
  read (10,*) a !
  read (10,*) a ! INCOMPACT3D Flow configuration
  read (10,*) a !
  read (10,*) jLES
  read (10,*) iin
  read (10,*) itype
  read (10,*) iturb
  read (10,*) ifirst
  read (10,*) ilast
  read (10,*) nscheme
  read (10,*) istret
  read (10,*) beta
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
  read (10,*) initstats1
  read (10,*) initstats2
  read (10,*) a !
  read (10,*) a ! NUMERICAL DISSIPATION
  read (10,*) a !
  read (10,*) fpi2
  read (10,*) a !
  read (10,*) a ! NUMERICAL TRIPPING
  read (10,*) a !
  read (10,*) itrip
  read (10,*) A_trip
  read (10,*) xs_tr
  read (10,*) ys_tr
  read (10,*) zs_tr
  read (10,*) ts_tr
  read (10,*) x0_tr
  close(10)

  if (nrank==0) then
    print *,'======================= Wall jet =========================='
    write(*,"(' initstats1         : ',I15)") initstats1
    write(*,"(' initstats2         : ',I15)") initstats2
    print *,'==========================================================='
  endif

  if (itrip==1) then
    if (nrank==0) then
      print *,'======================= Tripping ======================='
      write(*,"(' A_trip, x0_tr, ts_tr       : (',F6.2,',',F6.2,',',F6.2,')')") A_trip, x0_tr, ts_tr
      write(*,"(' xs_tr , ys_tr, zs_tr       : (',F6.2,',',F6.2,',',F6.2,')')") xs_tr, ys_tr, zs_tr
      print *,'========================================================'
    endif
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
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1,phis1,phiss1

  real(mytype) :: y,r,um,r3,x,z,h,ct
  real(mytype) :: cx0,cy0,cz0,hg,lg
  real(mytype) :: a,b,c,d,eta,profile,u
  integer :: k,j,i,ijk,fh,ierror,ii,is,code
  integer (kind=MPI_OFFSET_KIND) :: disp

  if (iscalar==1) then

    phi1 = zero !change as much as you want

    !do not delete this
    phis1=phi1
    phiss1=phis1

  endif

  ux1=zero;uy1=zero;uz1=zero

  !  if (iin.ne.0) then
  !     call system_clock(count=code)
  !     if (iin.eq.2) code=0
  !     call random_seed(size = ii)
  !     call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

  !     call random_number(ux1)
  !     call random_number(uy1)
  !     call random_number(uz1)
  !  endif

  if (iin.eq.1) then !prescribed profile

    if (itype==0) then

      if (nrank==0) write(*,*) "Glauert Wall Jet!"
      a = 1.77875559_mytype
      b = -1.797565198_mytype
      c = -0.40051025_mytype
      d = 1.97454561_mytype
      eta= 2.1266172089_mytype

      do k=1,xsize(3)
        do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)
          do i=1,xsize(1)

            profile=b*c*d*((y*eta)**(d-1.))*exp(c*(y*eta)**d)
            !coflow
            if (y.LT.0.519444444444) then
              u=profile
            else
              u=0.1_mytype*(1._mytype-profile)+profile
            endif

            ux1(i,j,k)=ux1(i,j,k)+ u

          enddo
        enddo
      enddo
    endif

    if (itype==1) then
      if (nrank==0) write(*,*) "Hyperbolic tangent Wall Jet!"
      do k=1,xsize(3)
        do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)
          do i=1,xsize(1)

            ux1(i,j,k)=ux1(i,j,k)+ (0.55+0.45*tanh(10*(1.-y)))*(tanh(10.*y))

          enddo
        enddo
      enddo
    endif

  endif

  !INIT FOR G AND U=MEAN FLOW + NOISE
  do k=1,xsize(3)
    do j=1,xsize(2)
      do i=1,xsize(1)
        ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
        uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
        uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
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
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1

  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx

  real(mytype) :: u,h,y,z,de,gapy,gapz,nondaz,nonday
  real(mytype) :: flut,temp,profile,basin,randn
  real(mytype) :: uxmax,uxmin,uxmax1,uxmin1
  real(mytype) :: a,b,c,d,e,eta,u_in,Av,v_w,y0,vv,wv,Ast,tau,ust
  real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,cz

  integer  :: i,k,j,is,ii,code

  udx=1./dx; udy=1./dy; udz=1./dz; uddx=0.5/dx; uddy=0.5/dy; uddz=0.5/dz

  uxmax=-1609.
  uxmin=1609.

  ux1 = ux1*(one-ep1)
  uy1 = uy1*(one-ep1)
  uz1 = uz1*(one-ep1)

  h = 1.         !entrance height
  de = 20.

  !!INFLOW

  if (itype==0) then !glauer jet

    a = 1.77875559_mytype
    b = -1.797565198_mytype
    c = -0.40051025_mytype
    d = 1.97454561_mytype
    eta= 2.1266172089_mytype

    do is=1,nphi
      do k=1,xsize(3)
        do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)

          !channel
          !u=5.*y*exp(-2.5*y**2)
          !u=4.65858*y*exp(-3.88762*y**2)
          !u=(0.55+0.45*tanh(10*(h-y)))*(tanh(10.*y))

          !coflow
          !if (y.LT.SQRT(1/(2*3.88762))) then
          !    u=profile
          !else
          !   u=0.01*(1-profile)+profile
          !endif

          !TBL
          ! u=1.-exp(-4.0*y) !1.5*y - 0.5*y**3
          !u = 1. - exp( -0.3*(4.91*y) - 0.1*(4.91*y)**2 )

          profile=b*c*d*((y*eta)**(d-1.))*exp(c*(y*eta)**d)

          !coflow
          if (y.LT.0.519444444444) then
            u=profile
          else
            u=0.1_mytype*(1._mytype-profile)+profile
          endif

          !TBL
          ! u=1.-exp(-4.0*y) !1.5*y - 0.5*y**3
          !u = 1. - exp( -0.3*(4.91*y) - 0.1*(4.91*y)**2 )

          bxx1(j,k) = u!b*c*d*((y*eta)**(d-1.))*exp(c*(y*eta)**d)

          !ahlman theoretical
          phi1(1,j,k,is)=(half + half * tanh(10. * (h - y))) * cp(is)


        enddo !z
      enddo !y
    enddo !nphi

  endif


  if (itype==1) then !hyperbolic tangent
    do is=1,nphi
      do k=1,xsize(3)
        do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)

          bxx1(j,k) = (0.55+0.45*tanh(10*(h-y)))*(tanh(10.*y))

          phi1(1,j,k,is)=zero !(half + half * tanh(10. * (h - y))) * cp(is)

        enddo !z
      enddo !y
    enddo !nphi

  endif


  if (iturb==1) then !PURE WHITE NOISE

    call system_clock(count=code)
    call random_seed(size = ii)
    call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

    call random_number(ux1)
    call random_number(uy1)
    call random_number(uz1)

    do is=1,nphi
      do k=1,xsize(3)
        z=(k+xstart(3)-1-1)*dz
        do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)
        enddo !z
      enddo !y
    enddo !nphi
  endif


  if (iturb==2) then !AHLMAN PERTURBATION

    call system_clock(count=code)
    call random_seed(size = ii)
    call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

    call random_number(bxy1)
    call random_number(bxz1)


    do is=1,nphi
      do k=1,xsize(3)
        z=(k+xstart(3)-1-1)*dz
        do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)

          u_in = 1.
          h = 1.

          Av = 0.01*u_in
          v_w = 0.6*h
          y0 = 0.6*h

          vv = +Av*cos(pi*z/v_w)*sin(pi*(y-y0)/v_w)
          wv = -Av*sin(pi*z/v_w)*cos(pi*(y-y0)/v_w)

          Ast = 0.01*u_in
          tau = h/u_in
          ust = (Ast/2.)*(sin(0.5*(t/tau))+sin(0.1*(t/tau)))

          u = (0.5+0.5*tanh(10.*(h-y)))*(tanh(10.*y))

          !randon disturbances + streamwise vortices + periodic disturbances
          bxy1(j,k) = u*( bxy1(j,k)*noise1 + vv + ust)
          bxz1(j,k) = u*( bxz1(j,k)*noise1 + wv + ust)

        enddo !z
      enddo !y
    enddo !nphi
  endif

  !INFLOW top !!!not fully tested
  if (nclyn==2) then
    if (xend(2)==ny) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          byxn(i,k)=ux1(i,xsize(2)-1,k)
          byyn(i,k)=uy1(i,xsize(2)-1,k)
          byzn(i,k)=uz1(i,xsize(2)-1,k)
        enddo
      enddo
    endif
  endif



  !OUTFLOW
  do k=1,xsize(3)
    do j=1,xsize(2)
      if (ux1(1,j,k).gt.uxmax) uxmax=ux1(1,j,k)
      if (ux1(1,j,k).lt.uxmin) uxmin=ux1(1,j,k)
    enddo
  enddo
  call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)
  if (nrank==0) write(*,*) "Inflow x velocity min max=",real(uxmin1,4),real(uxmax1,4)


  do k=1,xsize(3)
    do j=1,xsize(2)
      if (ux1(nx-1,j,k).gt.uxmax) uxmax=ux1(nx-1,j,k)
      if (ux1(nx-1,j,k).lt.uxmin) uxmin=ux1(nx-1,j,k)
    enddo
  enddo

  call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

  !SPONGE ZONE
  do k=1,xsize(3)
    do j=1,xsize(2)
      do i=1,xsize(1)
        if (ux1(i,j,k).lt.0.) ux1(i,j,k)=ux1(i,j,k)*0.5*(1.-tanh(real(i-1,mytype)*dx-xlx+3.0))
      enddo
    enddo
  enddo
  if (nrank==0) print *,'SPONGE FOR VELOCITY'
  !SPONGE ZONE


  if (nclyS1.eq.2) then
    do is=1,nphi
      if (xstart(2)==1) then
        do k=1,xsize(3)
          do i=1,xsize(1)
            phi1(i,1,k,is)=cp(is)
          enddo
        enddo
      endif
    enddo
  endif

  if (u1==0) cx=(0.5*(uxmax1+uxmin1))*gdt(itr)*udx
  if (u1==1) cx=uxmax1*gdt(itr)*udx
  if (u1==2) cx=u2*gdt(itr)*udx    !works better

  do k=1,xsize(3)
    do j=1,xsize(2)
      bxxn(j,k)=ux1(nx,j,k)-cx*(ux1(nx,j,k)-ux1(nx-1,j,k))
      bxyn(j,k)=uy1(nx,j,k)-cx*(uy1(nx,j,k)-uy1(nx-1,j,k))
      bxzn(j,k)=uz1(nx,j,k)-cx*(uz1(nx,j,k)-uz1(nx-1,j,k))
    enddo
  enddo

  if (iscalar==1) then
    do k=1,xsize(3)
      do j=1,xsize(2)
        phi1(nx,j,k,:)=phi1(nx,j,k,:)-cx*(phi1(nx,j,k,:)-phi1(nx-1,j,k,:))
      enddo
    enddo
  endif

  if (nrank==0) write(*,*) "Outflow velocity min max=",uxmin1,uxmax1

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
  integer :: nprobes, ntimes1, ntimes2
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes
  !
  real(mytype),save,allocatable,dimension(:,:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

contains

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: x, xprobes, yprobes, zprobes
    integer :: i,j,k,code
    character :: a

    allocate(usum(zsize(1),zsize(2)),vsum(zsize(1),zsize(2)),wsum(zsize(1),zsize(2)))
    allocate(uusum(zsize(1),zsize(2)),uvsum(zsize(1),zsize(2)),uwsum(zsize(1),zsize(2)))
    allocate(vvsum(zsize(1),zsize(2)),vwsum(zsize(1),zsize(2)),wwsum(zsize(1),zsize(2)))
    usum=zero;vsum=zero;wsum=zero
    uusum=zero;uvsum=zero;uwsum=zero
    vvsum=zero;vwsum=zero;wwsum=zero
    ntimes1 = 0
    ntimes2 = 0
    nprobes  = 0

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

  end subroutine init_post
  !############################################################################
  subroutine postprocessing(ux1,uy1,uz1,phi1,ep1,pre1,diss1) !By Felipe Schuch

    USE decomp_2d_io

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1, pre1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: diss1
    !
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) ::  uprime, vprime, wprime, dissprime, preprime
    !
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u1sum,v1sum,w1sum,ta1 !first
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u2sum,v2sum,w2sum !second
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u3sum,v3sum,w3sum !third
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u4sum,v4sum,w4sum !fourth
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: uvsum,uwsum,vwsum !tensors
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: disssum,presum,tsum !diss, pressure and auxiliary variable

    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u1psum,v1psum,w1psum !first
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u2psum,v2psum,w2psum !second
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u3psum,v3psum,w3psum !third
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: u4psum,v4psum,w4psum !fourth
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: uvpsum,uwpsum,vwpsum !tensors
    real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: dissspum,prepsum,tpsum !diss, pressure and auxiliary variable
    !
    real(mytype),dimension(xszS(1),xszS(2),xszS(3),nphi) :: psum,ppsum
    real(mytype),dimension(xszS(1),xszS(2),xszS(3),nphi) :: upsum,vpsum,wpsum
    !
    integer :: i,j,k,is,tacum
    character(len=30) :: filename

    tacum = itime - initstats1

    if (itime.ge.initstats1) then !start acumulating

      !pressure
      call fine_to_coarseS(1,pre1,tsum)
      presum=presum+tsum

      !u1=ux1
      call fine_to_coarseS(1,ux1,tsum)
      u1sum=u1sum+tsum

      !u2=ux1*ux1
      ta1=ux1*ux1
      call fine_to_coarseS(1,ta1,tsum)
      u2sum=u2sum+tsum

      !u3=ux1*ux1*ux1 (THIRD ORDER MOMENTS - SKEWNESS/assimetria)
      ta1=ux1*ux1*ux1
      call fine_to_coarseS(1,ta1,tsum)
      u3sum=u3sum+tsum

      !u4=ux1*ux1*ux1*ux1 (FOURTH ORDER MOMENTS - FLATNESS/achatamento)
      ta1=ux1*ux1*ux1*ux1
      call fine_to_coarseS(1,ta1,tsum)
      u4sum=u4sum+tsum

      !v1=uy1
      call fine_to_coarseS(1,uy1,tsum)
      v1sum=v1sum+tsum

      !v2=uy1*uy1
      ta1=uy1*uy1
      call fine_to_coarseS(1,ta1,tsum)
      v2sum=v2sum+tsum

      !v3sum=uy1*uy1*uy1 (THIRD ORDER MOMENTS - SKEWNESS/assimetria)
      ta1=uy1*uy1*uy1
      call fine_to_coarseS(1,ta1,tsum)
      v3sum=v3sum+tsum

      !v4sum=uy1*uy1*uy1*uy1 (FOURTH ORDER MOMENTS - FLATNESS/achatamento)
      ta1=uy1*uy1*uy1*uy1
      call fine_to_coarseS(1,ta1,tsum)
      v4sum=v4sum+tsum

      !w1=uz1
      call fine_to_coarseS(1,uz1,tsum)
      w1sum=w1sum+tsum

      !w2=uz1*uz1
      ta1=uz1*uz1
      call fine_to_coarseS(1,ta1,tsum)
      w2sum=w2sum+tsum

      !w3=uz1*uz1*uz1 (THIRD ORDER MOMENTS - SKEWNESS/assimetria)
      ta1=uz1*uz1*uz1
      call fine_to_coarseS(1,ta1,tsum)
      w3sum=w3sum+tsum

      !w4=uz1*uz1*uz1*uz1 (FOURTH ORDER MOMENTS - FLATNESS/achatamento)
      ta1=uz1*uz1*uz1*uz1
      call fine_to_coarseS(1,ta1,tsum)
      w4sum=w4sum+tsum

      !uvsum=ux1*uy1
      ta1=ux1*uy1
      call fine_to_coarseS(1,ta1,tsum)
      uvsum=uvsum+tsum

      !uwsum=ux1*uz1
      ta1=ux1*uz1
      call fine_to_coarseS(1,ta1,tsum)
      uwsum=uwsum+tsum

      !vwsum=uy1*uz1
      ta1=uy1*uz1
      call fine_to_coarseS(1,ta1,tsum)
      vwsum=vwsum+tsum

      if (iscalar==1) then
        do is=1, nphi

          !psum=phi1
          call fine_to_coarseS(1,phi1(:,:,:,is),tsum)
          psum(:,:,:,is)=psum(:,:,:,is)+tsum

          !ppsum=phi1*phi1
          ta1=phi1(:,:,:,is)*phi1(:,:,:,is)
          call fine_to_coarseS(1,ta1,tsum)
          ppsum(:,:,:,is)=ppsum(:,:,:,is)+tsum

          !upsum=phi1*ux1
          ta1=phi1(:,:,:,is)*ux1
          call fine_to_coarseS(1,ta1,tsum)
          upsum(:,:,:,is)=upsum(:,:,:,is)+tsum

          !vpsum=phi1*uy1
          ta1=phi1(:,:,:,is)*uy1
          call fine_to_coarseS(1,ta1,tsum)
          vpsum(:,:,:,is)=vpsum(:,:,:,is)+tsum

          !wpsum=phi1*uz1
          ta1=phi1(:,:,:,is)*uz1
          call fine_to_coarseS(1,ta1,tsum)
          wpsum(:,:,:,is)=wpsum(:,:,:,is)+tsum

        enddo
      endif
      
      endif
      if (itime.ge.initstats2) then

        !AVERAGE OF THE FLUCTUATIONS ON TIME


        preprime = pre1 - presum/tacum
        uprime = ux1 - u1sum/tacum
        vprime = uy1 - u2sum/tacum
        wprime = uz1 - u3sum/tacum

        !pressure
        call fine_to_coarseS(1,preprime,tpsum)
        prepsum=prepsum+tpsum

        !u1=ux1
        call fine_to_coarseS(1,uprime,tpsum)
        u1psum=u1sum+tpsum

        !u2=ux1*ux1
        call fine_to_coarseS(1,uprime*uprime,tpsum)
        u2psum=u2psum+tpsum

        !u3=ux1*ux1*ux1 (THIRD ORDER MOMENTS - SKEWNESS/assimetria)
        call fine_to_coarseS(1,uprime*uprime*uprime,tpsum)
        u3psum=u3psum+tpsum

        !u4=ux1*ux1*ux1*ux1 (FOURTH ORDER MOMENTS - FLATNESS/achatamento)
        call fine_to_coarseS(1,uprime*uprime*uprime*uprime,tpsum)
        u4psum=u4psum+tpsum

        !v1=uy1
        call fine_to_coarseS(1,vprime,tpsum)
        v1psum=v1psum+tpsum

        !v2=uy1*uy1
        call fine_to_coarseS(1,vprime*vprime,tpsum)
        v2psum=v2psum+tpsum

        !v3sum=uy1*uy1*uy1 (THIRD ORDER MOMENTS - SKEWNESS/assimetria)
        call fine_to_coarseS(1,vprime*vprime*vprime,tpsum)
        v3psum=v3psum+tpsum

        !v4sum=uy1*uy1*uy1*uy1 (FOURTH ORDER MOMENTS - FLATNESS/achatamento)
        call fine_to_coarseS(1,vprime*vprime*vprime*vprime,tpsum)
        v4psum=v4psum+tpsum

        !w1=uz1
        call fine_to_coarseS(1,wprime,tpsum)
        w1psum=w1psum+tpsum

        !w2=uz1*uz1
        call fine_to_coarseS(1,wprime*wprime,tpsum)
        w2psum=w2psum+tpsum

        !w3=uz1*uz1*uz1 (THIRD ORDER MOMENTS - SKEWNESS/assimetria)
        call fine_to_coarseS(1,wprime*wprime*wprime,tpsum)
        w3psum=w3psum+tpsum

        !w4=uz1*uz1*uz1*uz1 (FOURTH ORDER MOMENTS - FLATNESS/achatamento)
        call fine_to_coarseS(1,wprime*wprime*wprime*wprime,tpsum)
        w4psum=w4psum+tpsum

        !uvsum=ux1*uy1
        call fine_to_coarseS(1,uprime*vprime,tpsum)
        uvpsum=uvpsum+tpsum

        !uwsum=ux1*uz1
        call fine_to_coarseS(1,uprime*wprime,tpsum)
        uwpsum=uwpsum+tpsum

        !vwsum=uy1*uz1
        call fine_to_coarseS(1,vprime*wprime,tpsum)
        vwpsum=vwpsum+tpsum

     endif
     
     if (mod(itime,imodulo)==0) then

        if (save_pre==1) call decomp_2d_write_one(1,presum,'stats/pre.sum',1)

        if (save_ux==1) call decomp_2d_write_one(1,u1sum,'stats/u1.sum',1)
        if (save_ux==1) call decomp_2d_write_one(1,u2sum,'stats/u2.sum',1)
        if (save_ux==1) call decomp_2d_write_one(1,u3sum,'stats/u3.sum',1)
        if (save_ux==1) call decomp_2d_write_one(1,u4sum,'stats/u4.sum',1)

        if (save_uy==1) call decomp_2d_write_one(1,v1sum,'stats/v1.sum',1)
        if (save_uy==1) call decomp_2d_write_one(1,v2sum,'stats/v2.sum',1)
        if (save_uy==1) call decomp_2d_write_one(1,v3sum,'stats/v3.sum',1)
        if (save_uy==1) call decomp_2d_write_one(1,v4sum,'stats/v4.sum',1)

        if (save_uz==1) call decomp_2d_write_one(1,w1sum,'stats/w1.sum',1)
        if (save_uz==1) call decomp_2d_write_one(1,w2sum,'stats/w2.sum',1)
        if (save_uz==1) call decomp_2d_write_one(1,w3sum,'stats/w3.sum',1)
        if (save_uz==1) call decomp_2d_write_one(1,w4sum,'stats/w4.sum',1)

        call decomp_2d_write_one(1,half*(u2sum+v2sum+w2sum),'stats/k.sum',1)

        call decomp_2d_write_one(1,uvsum,'stats/uv.sum',1)
        call decomp_2d_write_one(1,uwsum,'stats/uw.sum',1)
        call decomp_2d_write_one(1,vwsum,'stats/vw.sum',1)

        !FLUCTUATIONS

        if (save_pre==1) call decomp_2d_write_one(1,prepsum,'stats/prep.sum',1)

        if (save_ux==1) call decomp_2d_write_one(1,u1psum,'stats/u1p.sum',1)
        if (save_ux==1) call decomp_2d_write_one(1,u2psum,'stats/u2p.sum',1)
        if (save_ux==1) call decomp_2d_write_one(1,u3psum,'stats/u3p.sum',1)
        if (save_ux==1) call decomp_2d_write_one(1,u4psum,'stats/u4p.sum',1)

        if (save_uy==1) call decomp_2d_write_one(1,v1psum,'stats/v1p.sum',1)
        if (save_uy==1) call decomp_2d_write_one(1,v2psum,'stats/v2p.sum',1)
        if (save_uy==1) call decomp_2d_write_one(1,v3psum,'stats/v3p.sum',1)
        if (save_uy==1) call decomp_2d_write_one(1,v4psum,'stats/v4p.sum',1)

        if (save_uz==1) call decomp_2d_write_one(1,w1psum,'stats/w1p.sum',1)
        if (save_uz==1) call decomp_2d_write_one(1,w2psum,'stats/w2p.sum',1)
        if (save_uz==1) call decomp_2d_write_one(1,w3psum,'stats/w3p.sum',1)
        if (save_uz==1) call decomp_2d_write_one(1,w4psum,'stats/w4p.sum',1)

        call decomp_2d_write_one(1,half*(u2psum+v2psum+w2psum),'stats/kp.sum',1)

        call decomp_2d_write_one(1,uvpsum,'stats/uvp.sum',1)
        call decomp_2d_write_one(1,uwpsum,'stats/uwp.sum',1)
        call decomp_2d_write_one(1,vwpsum,'stats/vwp.sum',1)

!        if (save_phi==1) then
!          if (iscalar==1) then
!            do is=1, nphi
!              call decomp_2d_write_one(1, psum(:,:,:,is), 'stats/c'//char(48+is)//'.sum',1)
!              call decomp_2d_write_one(1,ppsum(:,:,:,is),'stats/c'//char(48+is)//'c'//char(48+is)//'.sum',1)
!              call decomp_2d_write_one(1,upsum(:,:,:,is),'stats/uc'//char(48+is)//'.sum',1)
!              call decomp_2d_write_one(1,vpsum(:,:,:,is),'stats/vc'//char(48+is)//'.sum',1)
!              call decomp_2d_write_one(1,wpsum(:,:,:,is),'stats/wc'//char(48+is)//'.sum',1)

!              call decomp_2d_write_one(1, pppsum(:,:,:,is), 'stats/c'//char(48+is)//'p.sum',1)
!              call decomp_2d_write_one(1,ppppsum(:,:,:,is),'stats/c'//char(48+is)//'pc'//char(48+is)//'.sum',1)
!              call decomp_2d_write_one(1,uppsum(:,:,:,is),'stats/uc'//char(48+is)//p'.sum',1)
!              call decomp_2d_write_one(1,vppsum(:,:,:,is),'stats/vc'//char(48+is)//p'.sum',1)
!              call decomp_2d_write_one(1,wppsum(:,:,:,is),'stats/wc'//char(48+is)//p'.sum',1)

!            enddo
!          endif
!        endif
      endif

      return

    end subroutine postprocessing
    !############################################################################
    subroutine write_probes(ux1,uy1,uz1,pre1,diss1,phi1) !By Felipe Schuch

      real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1, pre1, diss1
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
  end module post_processing
#endif
