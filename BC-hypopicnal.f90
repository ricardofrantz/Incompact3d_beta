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

  iscalar = 1

  open(10,file='BC-hypopicnal.prm',status='old',form='formatted')
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
  read (10,*) wrotation
  read (10,*) irotation
  read (10,*) initstats1
  read (10,*) initstats2
  read (10,*) a !
  read (10,*) a ! NUMERICAL DISSIPATION
  read (10,*) a !
  read (10,*) jLES
  read (10,*) fpi2

  if (nrank==0) then
     print *,'==================Ekman flow==================='
     write(*,"(' irotation          : ',I15)") irotation
     write(*,"(' wrotation          : ',F15.8)") wrotation
     write(*,"(' initstats1         : ',I15)") initstats1
     write(*,"(' initstats2         : ',I15)") initstats2
     print *,'==========================================================='
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

  real(mytype) :: y,r,um,r3,x,z,h,ct
  real(mytype) :: cx0,cy0,cz0,hg,lg,delta_Ekman
  integer :: k,j,i,ijk,fh,ierror,ii,is,code
  integer (kind=MPI_OFFSET_KIND) :: disp

  integer, dimension (:), allocatable :: seed

  real(mytype), save, allocatable, dimension(:) :: uansorge, vansorge
  allocate(uansorge(ysize(2)),vansorge(ysize(2)))

  delta_Ekman=one

  if (iscalar==1) then

   do is=1,nphi

    !call random_number(phi1(:,:,:,is))

   do k=1,xsize(3)
    !z=(k+xstart(3)-1-1)*dz
    do j=1,xsize(2)
        if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy - 0.8*delta_Ekman
        if (istret.ne.0) y=yp(j+xstart(2)-1) - 0.8*delta_Ekman
      do i=1,xsize(1)
       !x=(i+xstart(1)-1-1)*dx!-xlx/2.
       !phi1(i,j,k,is)=exp(-7.*y*y)*noise*(2.*phi1(i,j,k,is)-1.)
      enddo
     enddo
    enddo

   do k=1,xsize(3)
    do j=1,xsize(2)
        if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
        if (istret.ne.0) y=yp(j+xstart(2)-1)
     do i=1,xsize(1)
      !condição incicial de empuxo segundo Ansorge 2014
       !phi1(i,j,k,is)=phi1(i,j,k,is) + erf(y/1.5_mytype)
       !condição incicial de temperatura, considerou-se um perfil de tarde com maximo junto ao fundo e mínimo em yly 
       !phi1(i,j,k,is)=(1.-y/yly)*cp(is)  
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
     !   do k=1,xsize(3)
     !      do i=1,xsize(1)
     !         ux1(i,2,k)=zero
     !         uy1(i,2,k)=zero
     !         uz1(i,2,k)=zero
     !      enddo
     !   enddo
     !   do k=1,xsize(3)
     !      do i=1,xsize(1)
     !         ux1(i,3,k)=zero
     !         uy1(i,3,k)=zero
     !         uz1(i,3,k)=zero
     !      enddo
     !   enddo
     !   do k=1,xsize(3)
     !      do i=1,xsize(1)
     !         ux1(i,4,k)=zero
     !         uy1(i,4,k)=zero
     !         uz1(i,4,k)=zero
     !      enddo
     !   enddo
     endif
   endif

  !modulation of the random noise + initial velocity profile
  do k=1,xsize(3)
     do j=1,xsize(2)
        if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy- 0.8*delta_Ekman
        if (istret.ne.0) y=yp(j+xstart(2)-1)- 0.8*delta_Ekman
        um=exp(-zptwo*y*y)
        do i=1,xsize(1)
           ux1(i,j,k)=noise*um*(two*ux1(i,j,k)-one)
           uy1(i,j,k)=noise*um*(two*uy1(i,j,k)-one)
           uz1(i,j,k)=noise*um*(two*uz1(i,j,k)-one)
        enddo
     enddo
  enddo



  endif

  if (iin==2) then !import fields as initial condition
   if (nrank==0) print *,'Reading external files"'

    call decomp_2d_read_one(1,ux1,'ux1_init.dat')
    call decomp_2d_read_one(1,uy1,'uy1_init.dat')
    call decomp_2d_read_one(1,uz1,'uz1_init.dat')

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

#ifdef DEBG
  if (nrank .eq. 0) print *,'# init end ok'
#endif

  return
end subroutine init
!********************************************************************
subroutine boundary_conditions (ux,uy,uz,phi,ep1)

  USE param
  USE variables
  USE decomp_2d

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


  if (nclyS1==2) then
     if (xstart(2)==1) then
        do k=1,xsize(3)
           do i=1,xsize(1)
             do is=1,nphi
           phi(i,1,k,is)=zero
             enddo
           enddo
        enddo
     endif
  endif

  if (nclySn==2) then
     if (xend(2)==ny) then
        do k=1,xsize(3)
           do i=1,xsize(1)
             do is=1,nphi
            phi(i,xsize(2),k,is)=cp(is)
             enddo
           enddo
        enddo
     endif
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
  real(mytype),save,allocatable,dimension(:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

contains

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: x, xprobes, yprobes, zprobes
    integer :: i,j,k,code
    character :: a

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post start'
#endif

    allocate(usum(ysize(2)),vsum(ysize(2)),wsum(ysize(2)))
    allocate(uusum(ysize(2)),uvsum(ysize(2)),uwsum(ysize(2)))
    allocate(vvsum(ysize(2)),vwsum(ysize(2)),wwsum(ysize(2)))
    usum=zero;vsum=zero;wsum=zero
    uusum=zero;uvsum=zero;uwsum=zero
    vvsum=zero;vwsum=zero;wwsum=zero
    ntimes1 = 0
    ntimes2 = 0

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
       nyprobes(:) = i+1
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

#ifdef DEBG 
    if (nrank .eq. 0) print *,'# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine postprocessing(ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

    USE MPI

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1
    !
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2, uy2, uz2
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) ::  uprime, vprime, wprime

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi2
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) ::  phiprime
    !
    real(mytype),dimension(ysize(2)) :: um, vm, wm
    real(mytype),dimension(ysize(2)) :: um1,vm1,wm1
    real(mytype),dimension(ysize(2)) :: uum, uvm, uwm, vvm, vwm, wwm
    real(mytype),dimension(ysize(2)) :: uum1,uvm1,uwm1,vvm1,vwm1,wwm1
    !
    integer :: i,j,k,code
    character(len=30) :: filename

     call VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu,pre1)
     call dissmed (ux1,uy1,uz1,phi1,diss1)
     call write_probes(ux1,uy1,uz1,pre1,diss1,phi1)

    if (itime.ge.initstats1) then

       call transpose_x_to_y(ux1,ux2)
       call transpose_x_to_y(uy1,uy2)
       call transpose_x_to_y(uz1,uz2)
       !call transpose_x_to_y(phi1,phi2())

       um = zero; vm = zero; wm = zero
       um1= zero; vm1= zero; wm1= zero
       do j=1,ysize(2)
          um(j) = sum(ux2(:,j,:))
          vm(j) = sum(uy2(:,j,:))
          wm(j) = sum(uz2(:,j,:))
       enddo

       call MPI_ALLREDUCE(um,um1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(vm,vm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(wm,wm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

       usum = usum + um1/real(nx*nz,mytype)
       vsum = vsum + vm1/real(nx*nz,mytype)
       wsum = wsum + wm1/real(nx*nz,mytype)

       ntimes1 = ntimes1 + 1

       if (itime.ge.initstats2) then

          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   uprime(i,j,k) = ux2(i,j,k) - usum(j)/real(ntimes1,mytype)
                   vprime(i,j,k) = uy2(i,j,k) - vsum(j)/real(ntimes1,mytype)
                   wprime(i,j,k) = uz2(i,j,k) - wsum(j)/real(ntimes1,mytype)
                enddo
             enddo
          enddo

          uum=zero ;uvm=zero ;uwm=zero ;vvm=zero ;vwm=zero ;wwm=zero
          uum1=zero;uvm1=zero;uwm1=zero;vvm1=zero;vwm1=zero;wwm1=zero
          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   uum(j) = uum(j) + uprime(i,j,k)*uprime(i,j,k)
                   uvm(j) = uvm(j) + uprime(i,j,k)*vprime(i,j,k)
                   uwm(j) = uwm(j) + uprime(i,j,k)*wprime(i,j,k)
                   vvm(j) = vvm(j) + vprime(i,j,k)*vprime(i,j,k)
                   vwm(j) = vwm(j) + vprime(i,j,k)*wprime(i,j,k)
                   wwm(j) = wwm(j) + wprime(i,j,k)*wprime(i,j,k)
                enddo
             enddo
          enddo

          call MPI_ALLREDUCE(uum,uum1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(uvm,uvm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(uwm,uwm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(vvm,vvm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(vwm,vwm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(wwm,wwm1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

          uusum = uusum + uum1/real(nx*nz,mytype)
          uvsum = uvsum + uvm1/real(nx*nz,mytype)
          uwsum = uwsum + uwm1/real(nx*nz,mytype)
          vvsum = vvsum + vvm1/real(nx*nz,mytype)
          vwsum = vwsum + vwm1/real(nx*nz,mytype)
          wwsum = wwsum + wwm1/real(nx*nz,mytype)

          ntimes2 = ntimes2 + 1
       endif

       if (mod(itime,imodulo).eq.0) then !write results
          if (nrank.eq.0) then
             write(filename,"('./out/stats',I4.4)") itime/imodulo
             open(67,file=trim(filename),status='unknown',form='formatted')
             do j=1,ysize(2)
                write(67,'(10E16.8)') yp(j),&
                     usum(j)/real(ntimes1,mytype),&
                     vsum(j)/real(ntimes1,mytype),&
                     wsum(j)/real(ntimes1,mytype),&
                     uusum(j)/real(ntimes2,mytype),&
                     uvsum(j)/real(ntimes2,mytype),&
                     uwsum(j)/real(ntimes2,mytype),&
                     vvsum(j)/real(ntimes2,mytype),&
                     vwsum(j)/real(ntimes2,mytype),&
                     wwsum(j)/real(ntimes2,mytype)
             enddo
             close(67)
          endif
       endif
    endif

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
          write(filename,"('./out/probe',I4.4)") i
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
end module post_processing


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
                 diss1(i,j,k)=diss1(i,j,k)+2._mytype*xnu*0.5_mytype*0.5_mytype*(A(l,m,i,j,k)+A(m,l,i,j,k))**2
              enddo
           enddo
        enddo
     enddo
  enddo
  return

end subroutine dissipation
!############################################################################
subroutine dissmed (ux1,uy1,uz1,phi1,diss1)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  USE MPI

  implicit none

  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1

  real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3)) :: diss1

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxf1,uyf1,uzf1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,temp1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,temp2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,temp3
  
  real(mytype),dimension(ysize(2)) :: sumxz, sumxz1
  real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu

  real(8) :: tstart,tend
  integer :: i,j,k,is,code
  character(len=30) :: filename

  tstart=MPI_WTIME()
  
  call transpose_x_to_y(ux1,td2)
  call transpose_x_to_y(uy1,te2)
  call transpose_x_to_y(uz1,tf2)
  
  sumxz=0._mytype;sumxz1=0._mytype
  do j=1,ysize(2)
     sumxz(j) = sum(td2(:,j,:))
  enddo
  call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  sumxz1=sumxz1/real(nx*nz,mytype)
  do k=1,ysize(3)
   do j=1,ysize(2)
    do i=1,ysize(1)
         td2(i,j,k) = td2(i,j,k) - sumxz1(j)
    enddo
   enddo
  enddo
  
  sumxz=0._mytype;sumxz1=0._mytype
  do j=1,ysize(2)
     sumxz(j) = sum(te2(:,j,:))
  enddo
  call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  sumxz1=sumxz1/real(nx*nz,mytype)
  do k=1,ysize(3)
   do j=1,ysize(2)
    do i=1,ysize(1)
         te2(i,j,k) = te2(i,j,k) - sumxz1(j)
    enddo
   enddo
  enddo
  
  sumxz=0._mytype;sumxz1=0._mytype
    do j=1,ysize(2)
     sumxz(j) = sum(tf2(:,j,:))
  enddo
  call MPI_ALLREDUCE(sumxz,sumxz1,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  sumxz1=sumxz1/real(nx*nz,mytype)
  do k=1,ysize(3)
   do j=1,ysize(2)
    do i=1,ysize(1)
         tf2(i,j,k) = tf2(i,j,k) - sumxz1(j)
    enddo
   enddo
  enddo

  call transpose_y_to_x(td2,uxf1)
  call transpose_y_to_x(te2,uyf1)
  call transpose_y_to_x(tf2,uzf1)
  !x-derivatives
  call derx (ta1,uxf1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) !du/dx
  call derx (tb1,uyf1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) !dv/dx
  call derx (tc1,uzf1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) !dw/dx

  !y-derivatives
  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)

  call transpose_y_to_x(ta2,td1) !du/dy
  call transpose_y_to_x(tb2,te1) !dv/dy
  call transpose_y_to_x(tc2,tf1) !dw/dy

  !!z-derivatives
  call transpose_y_to_z(td2,td3)
  call transpose_y_to_z(te2,te3)
  call transpose_y_to_z(tf2,tf3)
  call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

  call transpose_z_to_y(ta3,td2)
  call transpose_z_to_y(tb3,te2)
  call transpose_z_to_y(tc3,tf2)

  call transpose_y_to_x(td2,tg1) !du/dy
  call transpose_y_to_x(te2,th1) !dv/dy
  call transpose_y_to_x(tf2,ti1) !dw/dy

  !du/dx=ta1 du/dy=td1 and du/dz=tg1
  !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
  !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

     call dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1) !sij sij

        uvisu=0._mytype
        call fine_to_coarseV(1,diss1,uvisu)
        write(filename,"('./data/dissf',I4.4)") itime/imodulo
        call decomp_2d_write_one(1,uvisu,filename,2)


        call transpose_x_to_y (diss1,temp2)
        call transpose_y_to_z (temp2,temp3)
        call mean_plane_z(temp3,zsize(1),zsize(2),zsize(3),temp3(:,:,1))
        write(filename,"('./data/dissfm',I4.4)") itime/imodulo
        call decomp_2d_write_plane(3,temp3,3,1,filename)

     !di1=0._mytype
     !write(filename,"('./data/utmap',I4.4)") itime/imodulo
     !do ijk=1,nvect1
     !   di1(ijk,1,1)=sqrt(sqrt((td1(ijk,1,1)**2)+(tf1(ijk,1,1)**2))*xnu)
     !enddo
     !call decomp_2d_write_plane(1,di1,2,1,filename)


  tend=MPI_WTIME()
  if (nrank.eq.0) print *,'Time in dissmean (s)', real(tend-tstart,4)

  return

end subroutine dissmed
#endif
