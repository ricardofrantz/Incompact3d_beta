#ifdef IMPLICIT
!#define my_mod_solide

module ludecomp

  use decomp_2d, only : mytype

  implicit none

  private
  public :: ludecomp7

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Interface for septadiag LU decomposition used in implicit mode
  ! ludecomp7_0 is called when ncly=0
  ! ludecomp7_12 is called when ncly=1 or 2
  !
  interface ludecomp7
    module procedure ludecomp7_12
      module procedure ludecomp7_0
      end interface ludecomp7

    contains

      !*******************************************************************
      !Decomposition LU matrice septadiag
      subroutine ludecomp7_12(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ny)

        use decomp_2d, only : mytype

        implicit none

        integer :: i,j,k,ny
        real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm
        real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm
        real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm

        real(mytype), dimension(ny,ny) :: mata, matl, matu, prod

        ggm=zero;hhm=zero;ssm=zero
        vvm=zero;wwm=zero;zzm=zero

        ggm(1)=aam(1)
        hhm(1)=bbm(1)
        ssm(1)=ccm(1)

        zzm(2)=ddm(2)/ggm(1)
        ggm(2)=aam(2)-zzm(2)*hhm(1)
        hhm(2)=bbm(2)-zzm(2)*ssm(1)
        ssm(2)=ccm(2)-zzm(2)*rrm(1)

        wwm(3)=eem(3)/ggm(1)
        zzm(3)=(ddm(3)-wwm(3)*hhm(1))/ggm(2)
        ggm(3)=aam(3)-wwm(3)*ssm(1)-zzm(3)*hhm(2)
        hhm(3)=bbm(3)-wwm(3)*rrm(1)-zzm(3)*ssm(2)
        ssm(3)=ccm(3)-zzm(3)*rrm(2)

        do j=4,ny
          vvm(j)=qqm(j)/ggm(j-3)
          wwm(j)=(eem(j)-vvm(j)*hhm(j-3))/ggm(j-2)
          zzm(j)=(ddm(j)-wwm(j)*hhm(j-2)-vvm(j)*ssm(j-3))/ggm(j-1)
          ggm(j)=aam(j)-vvm(j)*rrm(j-3)-wwm(j)*ssm(j-2)-zzm(j)*hhm(j-1)
          hhm(j)=bbm(j)-wwm(j)*rrm(j-2)-zzm(j)*ssm(j-1)
          ssm(j)=ccm(j)-zzm(j)*rrm(j-1)
        enddo

#ifdef DEBUG
        if (.false.) then
          print *,'aam,bbm,ccm,ddm,eem,qqm,rrm'
          do k=1,ny
            print *,aam(k),bbm(k),ccm(k),ddm(k),eem(k),qqm(k),rrm(k)
          enddo
          print *,'ggm,hhm,ssm,rrm'
          do k=1,ny
            print *,ggm(k),hhm(k),ssm(k),rrm(k)
          enddo
          print *,'vvm,wwm,zzm'
          do k=1,ny
            print *,vvm(k),wwm(k),zzm(k)
          enddo
          mata=zero
          matl=zero
          matu=zero
          j=1
          mata(j,ny-2)=qqm(j)
          mata(j,ny-1)=eem(j)
          mata(j,ny  )=ddm(j)
          mata(j,j   )=aam(j)
          mata(j,j+1 )=bbm(j)
          mata(j,j+2 )=ccm(j)
          mata(j,j+3 )=rrm(j)
          matl(j,j)=one
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          matu(j,j+2)=ssm(j)
          matu(j,j+3)=rrm(j)
          j=2
          mata(j,ny-1)=qqm(j)
          mata(j,ny  )=eem(j)
          mata(j,j-1 )=ddm(j)
          mata(j,j   )=aam(j)
          mata(j,j+1 )=bbm(j)
          mata(j,j+2 )=ccm(j)
          mata(j,j+3 )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          matu(j,j+2)=ssm(j)
          matu(j,j+3)=rrm(j)
          j=3
          mata(j,ny  )=qqm(j)
          mata(j,j-2 )=eem(j)
          mata(j,j-1 )=ddm(j)
          mata(j,j   )=aam(j)
          mata(j,j+1 )=bbm(j)
          mata(j,j+2 )=ccm(j)
          mata(j,j+3 )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matl(j,j-2)=wwm(j)
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          matu(j,j+2)=ssm(j)
          matu(j,j+3)=rrm(j)
          do j=4,ny-3
            mata(j,j-3)=qqm(j)
            mata(j,j-2)=eem(j)
            mata(j,j-1)=ddm(j)
            mata(j,j  )=aam(j)
            mata(j,j+1)=bbm(j)
            mata(j,j+2)=ccm(j)
            mata(j,j+3)=rrm(j)
            matl(j,j)=one
            matl(j,j-1)=zzm(j)
            matl(j,j-2)=wwm(j)
            matl(j,j-3)=vvm(j)
            matu(j,j)=ggm(j)
            matu(j,j+1)=hhm(j)
            matu(j,j+2)=ssm(j)
            matu(j,j+3)=rrm(j)
          enddo
          j=ny-2
          mata(j,j-3)=qqm(j)
          mata(j,j-2)=eem(j)
          mata(j,j-1)=ddm(j)
          mata(j,j  )=aam(j)
          mata(j,j+1)=bbm(j)
          mata(j,j+2)=ccm(j)
          mata(j,1  )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matl(j,j-2)=wwm(j)
          matl(j,j-3)=vvm(j)
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          matu(j,j+2)=ssm(j)
          j=ny-1
          mata(j,j-3)=qqm(j)
          mata(j,j-2)=eem(j)
          mata(j,j-1)=ddm(j)
          mata(j,j  )=aam(j)
          mata(j,j+1)=bbm(j)
          mata(j,1  )=ccm(j)
          mata(j,2  )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matl(j,j-2)=wwm(j)
          matl(j,j-3)=vvm(j)
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          j=ny
          mata(j,j-3)=qqm(j)
          mata(j,j-2)=eem(j)
          mata(j,j-1)=ddm(j)
          mata(j,j  )=aam(j)
          mata(j,1  )=bbm(j)
          mata(j,2  )=ccm(j)
          mata(j,3  )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matl(j,j-2)=wwm(j)
          matl(j,j-3)=vvm(j)
          matu(j,j)=ggm(j)

          print *,'A = '
          do j=1,ny
            print *,mata(j,:)
          enddo
          print *,'L = '
          do j=1,ny
            print *,matl(j,:)
          enddo
          print *,'U = '
          do j=1,ny
            print *,matu(j,:)
          enddo

          do i=1,ny
            do j=1,ny
              prod(i,j)=zero
              do k=1,ny
                prod(i,j)=prod(i,j)+matl(i,k)*matu(k,j)
              enddo
            enddo
          enddo
          print *,'A-L*U', maxval(abs(mata-prod)), maxloc(abs(mata-prod))
          do j=1,ny
            print *,mata(j,:)-prod(j,:)
          enddo
        endif
#endif

      end subroutine ludecomp7_12
      !*******************************************************************
      !Decomposition LU matrice cyclique septadiag
      subroutine ludecomp7_0(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,ny)

        use decomp_2d, only : mytype

        implicit none

        integer :: i,j,k,ny
        real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm
        real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm
        real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm
        real(mytype),dimension(ny), intent(out) :: l1m,l2m,l3m
        real(mytype),dimension(ny), intent(out) :: u1m,u2m,u3m

        real(mytype), dimension(ny,ny) :: mata, matl, matu, prod

        ! a=diag, b=diag+1, c=diag+2, r=diag+3
        !         d=diag-1, e=diag-2, q=diag-3
        ggm=zero;hhm=zero;ssm=zero ! U part (diag, diag+1, diag+2, diag+3=rrm)
        u1m=zero;u2m=zero;u3m=zero ! U cyclic part (last col, last-1, last-2)
        vvm=zero;wwm=zero;zzm=zero ! L part (diag=1, diag-3, diag-2, diag-1)
        l1m=zero;l2m=zero;l3m=zero ! L cyclic part (last row, last-1, last-2)

        ggm(1)=aam(1)
        hhm(1)=bbm(1)
        ssm(1)=ccm(1)
        u1m(1)=ddm(1)
        u2m(1)=eem(1)
        u3m(1)=qqm(1)
        l1m(1)=bbm(ny)/ggm(1)
        l2m(1)=ccm(ny-1)/ggm(1)
        l3m(1)=rrm(ny-2)/ggm(1)

        zzm(2)=ddm(2)/ggm(1)
        ggm(2)=aam(2)-zzm(2)*hhm(1)
        hhm(2)=bbm(2)-zzm(2)*ssm(1)
        ssm(2)=ccm(2)-zzm(2)*rrm(1)
        u1m(2)=eem(2)-zzm(2)*u1m(1)
        u2m(2)=qqm(2)-zzm(2)*u2m(1)
        u3m(2)=      -zzm(2)*u3m(1)
        l1m(2)=(ccm(ny)  -hhm(1)*l1m(1))/ggm(2)
        l2m(2)=(rrm(ny-1)-hhm(1)*l2m(1))/ggm(2)
        l3m(2)=(         -hhm(1)*l3m(1))/ggm(2)

        wwm(3)=eem(3)/ggm(1)
        zzm(3)=(ddm(3)-wwm(3)*hhm(1))/ggm(2)
        ggm(3)=aam(3)-wwm(3)*ssm(1)-zzm(3)*hhm(2)
        hhm(3)=bbm(3)-wwm(3)*rrm(1)-zzm(3)*ssm(2)
        ssm(3)=ccm(3)              -zzm(3)*rrm(2)
        u1m(3)=qqm(3)-wwm(3)*u1m(1)-zzm(3)*u1m(2)
        u2m(3)=      -wwm(3)*u2m(1)-zzm(3)*u2m(2)
        u3m(3)=      -wwm(3)*u3m(1)-zzm(3)*u3m(2)
        l1m(3)=(rrm(ny)-ssm(1)*l1m(1)-hhm(2)*l1m(2))/ggm(3)
        l2m(3)=(       -ssm(1)*l2m(1)-hhm(2)*l2m(2))/ggm(3)
        l3m(3)=(       -ssm(1)*l3m(1)-hhm(2)*l3m(2))/ggm(3)

        do j=4,ny-3
          vvm(j)=qqm(j)/ggm(j-3)
          wwm(j)=(eem(j)-vvm(j)*hhm(j-3))/ggm(j-2)
          zzm(j)=(ddm(j)-wwm(j)*hhm(j-2)-vvm(j)*ssm(j-3))/ggm(j-1)
          ggm(j)=aam(j)-vvm(j)*rrm(j-3)-wwm(j)*ssm(j-2)-zzm(j)*hhm(j-1)
          hhm(j)=bbm(j)-wwm(j)*rrm(j-2)-zzm(j)*ssm(j-1)
          ssm(j)=ccm(j)-zzm(j)*rrm(j-1)
          u1m(j)= - vvm(j)*u1m(j-3) - wwm(j)*u1m(j-2) - zzm(j)*u1m(j-1)
          u2m(j)= - vvm(j)*u2m(j-3) - wwm(j)*u2m(j-2) - zzm(j)*u2m(j-1)
          u3m(j)= - vvm(j)*u3m(j-3) - wwm(j)*u3m(j-2) - zzm(j)*u3m(j-1)
          l1m(j)= - (rrm(j-3)*l1m(j-3)+ssm(j-2)*l1m(j-2)+hhm(j-1)*l1m(j-1))/ggm(j)
          l2m(j)= - (rrm(j-3)*l2m(j-3)+ssm(j-2)*l2m(j-2)+hhm(j-1)*l2m(j-1))/ggm(j)
          l3m(j)= - (rrm(j-3)*l3m(j-3)+ssm(j-2)*l3m(j-2)+hhm(j-1)*l3m(j-1))/ggm(j)
        enddo

        do j=ny-2,ny
          vvm(j)=qqm(j)/ggm(j-3)
          wwm(j)=(eem(j)-vvm(j)*hhm(j-3))/ggm(j-2)
          zzm(j)=(ddm(j)-wwm(j)*hhm(j-2)-vvm(j)*ssm(j-3))/ggm(j-1)
          ggm(j)=aam(j)-vvm(j)*rrm(j-3)-wwm(j)*ssm(j-2)-zzm(j)*hhm(j-1)
          hhm(j)=bbm(j)-wwm(j)*rrm(j-2)-zzm(j)*ssm(j-1)
          ssm(j)=ccm(j)-zzm(j)*rrm(j-1)
        enddo
        j=ny-2
        u1m(j)= - vvm(j)*u1m(j-3)-wwm(j)*u1m(j-2)-zzm(j)*u1m(j-1) &
        - l3m(j-1)*rrm(j-1)
        u2m(j)= - vvm(j)*u2m(j-3)-wwm(j)*u2m(j-2)-zzm(j)*u2m(j-1) &
        - l3m(j-2)*rrm(j-2) - l3m(j-1)*ssm(j-1)
        u3m(j)= - vvm(j)*u3m(j-3)-wwm(j)*u3m(j-2)-zzm(j)*u3m(j-1) &
        - l3m(j-3)*rrm(j-3) - l3m(j-2)*ssm(j-2) - l3m(j-1)*hhm(j-1)
        do k=1,j-1
          u1m(j)=u1m(j) - u1m(k)*l3m(k)
          u2m(j)=u2m(j) - u2m(k)*l3m(k)
          u3m(j)=u3m(j) - u3m(k)*l3m(k)
        enddo
        l1m(j)= - ( rrm(j-3)*l1m(j-3)+ssm(j-2)*l1m(j-2)+hhm(j-1)*l1m(j-1) &
        + vvm(ny)*u3m(j-1)+wwm(ny)*u3m(j) ) / ( ggm(j)+u3m(j) )
        l2m(j)= - (rrm(j-3)*l2m(j-3)+ssm(j-2)*l2m(j-2)+hhm(j-1)*l2m(j-1) &
        + vvm(ny-1)*u3m(j-2)+wwm(ny-1)*u3m(j-1)+zzm(ny-1)*u3m(j) ) / (ggm(j)+u3m(j))
        do k=1,j-1
          l1m(j)=l1m(j) - l1m(k)*u3m(k) / (ggm(j)+u3m(j))
          l2m(j)=l2m(j) - l2m(k)*u3m(k) / (ggm(j)+u3m(j))
        enddo
        j=ny-1
        u1m(j)= - vvm(j)*u1m(j-3)-wwm(j)*u1m(j-2)-zzm(j)*u1m(j-1) &
        - l2m(j-2)*rrm(j-2) - l2m(j-1)*ssm(j-1)
        u2m(j)= - vvm(j)*u2m(j-3)-wwm(j)*u2m(j-2)-zzm(j)*u2m(j-1) &
        - l2m(j-3)*rrm(j-3) - l2m(j-2)*ssm(j-2) - l2m(j-1)*hhm(j-1)
        do k=1,j-1
          u1m(j)=u1m(j) - u1m(k)*l2m(k)
          u2m(j)=u2m(j) - u2m(k)*l2m(k)
        enddo
        l1m(j)= - (rrm(j-3)*l1m(j-3)+ssm(j-2)*l1m(j-2)+hhm(j-1)*l1m(j-1) &
        + vvm(ny)*u2m(j-2)+wwm(ny)*u2m(j-1)+zzm(ny)*u2m(j) ) / (ggm(j)+u2m(j));
        do k=1,j-1
          l1m(j)=l1m(j) - l1m(k)*u2m(k) / (ggm(j)+u2m(j));
        enddo
        j=ny
        u1m(j)= - vvm(j)*u1m(j-3)-wwm(j)*u1m(j-2)-zzm(j)*u1m(j-1) &
        - l1m(j-3)*rrm(j-3)-l1m(j-2)*ssm(j-2)-l1m(j-1)*hhm(j-1);
        do k=1,j-1
          u1m(j)=u1m(j) - u1m(k)*l1m(k)
        enddo

#ifdef DEBUG
        if (.false.) then
          print *,'aam,bbm,ccm,ddm,eem,qqm,rrm'
          do k=1,ny
            print *,aam(k),bbm(k),ccm(k),ddm(k),eem(k),qqm(k),rrm(k)
          enddo
          print *,'ggm,hhm,ssm,rrm,u1m,u2m,u3m'
          do k=1,ny
            print *,ggm(k),hhm(k),ssm(k),rrm(k),u1m(k),u2m(k),u3m(k)
          enddo
          print *,'vvm,wwm,zzm,l1m,l2m,l3m'
          do k=1,ny
            print *,vvm(k),wwm(k),zzm(k),l1m(k),l2m(k),l3m(k)
          enddo
          mata=zero
          matl=zero
          matu=zero
          j=1
          mata(j,ny-2)=qqm(j)
          mata(j,ny-1)=eem(j)
          mata(j,ny  )=ddm(j)
          mata(j,j   )=aam(j)
          mata(j,j+1 )=bbm(j)
          mata(j,j+2 )=ccm(j)
          mata(j,j+3 )=rrm(j)
          matl(j,j)=one
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          matu(j,j+2)=ssm(j)
          matu(j,j+3)=rrm(j)
          j=2
          mata(j,ny-1)=qqm(j)
          mata(j,ny  )=eem(j)
          mata(j,j-1 )=ddm(j)
          mata(j,j   )=aam(j)
          mata(j,j+1 )=bbm(j)
          mata(j,j+2 )=ccm(j)
          mata(j,j+3 )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          matu(j,j+2)=ssm(j)
          matu(j,j+3)=rrm(j)
          j=3
          mata(j,ny  )=qqm(j)
          mata(j,j-2 )=eem(j)
          mata(j,j-1 )=ddm(j)
          mata(j,j   )=aam(j)
          mata(j,j+1 )=bbm(j)
          mata(j,j+2 )=ccm(j)
          mata(j,j+3 )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matl(j,j-2)=wwm(j)
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          matu(j,j+2)=ssm(j)
          matu(j,j+3)=rrm(j)
          do j=4,ny-3
            mata(j,j-3)=qqm(j)
            mata(j,j-2)=eem(j)
            mata(j,j-1)=ddm(j)
            mata(j,j  )=aam(j)
            mata(j,j+1)=bbm(j)
            mata(j,j+2)=ccm(j)
            mata(j,j+3)=rrm(j)
            matl(j,j)=one
            matl(j,j-1)=zzm(j)
            matl(j,j-2)=wwm(j)
            matl(j,j-3)=vvm(j)
            matu(j,j)=ggm(j)
            matu(j,j+1)=hhm(j)
            matu(j,j+2)=ssm(j)
            matu(j,j+3)=rrm(j)
          enddo
          j=ny-2
          mata(j,j-3)=qqm(j)
          mata(j,j-2)=eem(j)
          mata(j,j-1)=ddm(j)
          mata(j,j  )=aam(j)
          mata(j,j+1)=bbm(j)
          mata(j,j+2)=ccm(j)
          mata(j,1  )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matl(j,j-2)=wwm(j)
          matl(j,j-3)=vvm(j)
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          matu(j,j+2)=ssm(j)
          j=ny-1
          mata(j,j-3)=qqm(j)
          mata(j,j-2)=eem(j)
          mata(j,j-1)=ddm(j)
          mata(j,j  )=aam(j)
          mata(j,j+1)=bbm(j)
          mata(j,1  )=ccm(j)
          mata(j,2  )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matl(j,j-2)=wwm(j)
          matl(j,j-3)=vvm(j)
          matu(j,j)=ggm(j)
          matu(j,j+1)=hhm(j)
          j=ny
          mata(j,j-3)=qqm(j)
          mata(j,j-2)=eem(j)
          mata(j,j-1)=ddm(j)
          mata(j,j  )=aam(j)
          mata(j,1  )=bbm(j)
          mata(j,2  )=ccm(j)
          mata(j,3  )=rrm(j)
          matl(j,j)=one
          matl(j,j-1)=zzm(j)
          matl(j,j-2)=wwm(j)
          matl(j,j-3)=vvm(j)
          matu(j,j)=ggm(j)
          do j=1,ny
            matu(j,ny)=matu(j,ny)+u1m(j)
          enddo
          do j=1,ny-1
            matu(j,ny-1)=matu(j,ny-1)+u2m(j)
          enddo
          do j=1,ny-2
            matu(j,ny-2)=matu(j,ny-2)+u3m(j)
          enddo
          do j=1,ny-1
            matl(ny,j)=matl(ny,j)+l1m(j)
          enddo
          do j=1,ny-2
            matl(ny-1,j)=matl(ny-1,j)+l2m(j)
          enddo
          do j=1,ny-3
            matl(ny-2,j)=matl(ny-2,j)+l3m(j)
          enddo
          print *,'A = '
          do j=1,ny
            print *,mata(j,:)
          enddo
          print *,'L = '
          do j=1,ny
            print *,matl(j,:)
          enddo
          print *,'U = '
          do j=1,ny
            print *,matu(j,:)
          enddo

          do i=1,ny
            do j=1,ny
              prod(i,j)=zero
              do k=1,ny
                prod(i,j)=prod(i,j)+matl(i,k)*matu(k,j)
              enddo
            enddo
          enddo
          print *,'A-L*U', maxval(abs(mata-prod)), maxloc(abs(mata-prod))
          do j=1,ny
            print *,mata(j,:)-prod(j,:)
          enddo
        endif
#endif

      end subroutine ludecomp7_0
      !

    end module ludecomp

    module matinv

      use decomp_2d, only : mytype

      implicit none

      private
      public :: septinv

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Interface for septadiag matrix inversion used in implicit mode
      ! septinv_0 is called when ncly=0
      ! septinv_12 is called when ncly=1 or 2
      !
      interface septinv
        module procedure septinv_12
          module procedure septinv_0
          end interface septinv

        contains

          !********************************************************************
          !INVERSE SEPTADIAG
          subroutine septinv_12(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,nx,ny,nz)
            !
            !********************************************************************
            use decomp_2d, only : mytype

            implicit none

            integer :: i,j,k,nx,ny,nz
            real(mytype),dimension(nx,ny,nz),intent(out) :: xsol
            real(mytype),dimension(nx,ny,nz),intent(in) :: bbb
            real(mytype),dimension(ny),intent(in)    :: ggm,hhm,ssm,rrm
            real(mytype),dimension(ny),intent(in)    :: vvm,wwm,zzm

            do k=1,nz
              do i=1,nx

                !Descente
                xsol(i,1,k)=bbb(i,1,k)
                xsol(i,2,k)=bbb(i,2,k)-zzm(2)*xsol(i,1,k)
                xsol(i,3,k)=bbb(i,3,k)-zzm(3)*xsol(i,2,k)-wwm(3)*xsol(i,1,k)

                do j=4,ny
                  xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
                  -zzm(j)*xsol(i,j-1,k)
                enddo
                !

                !Remontee
                xsol(i,ny,k)=xsol(i,ny,k)/ggm(ny)
                xsol(i,ny-1,k)=(xsol(i,ny-1,k)-hhm(ny-1)*xsol(i,ny,k))/ggm(ny-1)
                xsol(i,ny-2,k)=(xsol(i,ny-2,k)-hhm(ny-2)*xsol(i,ny-1,k)- &
                ssm(ny-2)*xsol(i,ny,k))/ggm(ny-2)

                do j=ny-3,1,-1
                  xsol(i,j,k)=(xsol(i,j,k)-hhm(j)*xsol(i,j+1,k)-ssm(j)*xsol(i,j+2,k) &
                  -rrm(j)*xsol(i,j+3,k))/ggm(j)
                enddo

              enddo
            enddo

          end subroutine septinv_12
          !********************************************************************
          !INVERSE SEPTADIAG CYCLIQUE
          subroutine septinv_0(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,nx,ny,nz)

            use decomp_2d, only : mytype

            implicit none

            integer :: i,j,k,kk,nx,ny,nz
            real(mytype),dimension(nx,ny,nz), intent(out) :: xsol
            real(mytype),dimension(nx,ny,nz), intent(in) :: bbb
            real(mytype),dimension(ny), intent(in)    :: ggm,hhm,ssm,rrm
            real(mytype),dimension(ny), intent(in)    :: vvm,wwm,zzm
            real(mytype),dimension(ny), intent(in) :: l1m,l2m,l3m
            real(mytype),dimension(ny), intent(in) :: u1m,u2m,u3m

            do k=1,nz
              do i=1,nx

                !Descente avec matrice triangulaire inf
                xsol(i,1,k)=bbb(i,1,k)
                xsol(i,2,k)=bbb(i,2,k)-zzm(2)*xsol(i,1,k)
                xsol(i,3,k)=bbb(i,3,k)-zzm(3)*xsol(i,2,k)-wwm(3)*xsol(i,1,k)
                do j=4,ny-3
                  xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
                  -zzm(j)*xsol(i,j-1,k)
                enddo
                j=ny-2
                xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
                -zzm(j)*xsol(i,j-1,k)
                do kk=1,j-1
                  xsol(i,j,k)=xsol(i,j,k)-l3m(kk)*xsol(i,kk,k)
                enddo
                j=ny-1
                xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
                -zzm(j)*xsol(i,j-1,k)
                do kk=1,j-1
                  xsol(i,j,k)=xsol(i,j,k)-l2m(kk)*xsol(i,kk,k)
                enddo
                j=ny
                xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
                -zzm(j)*xsol(i,j-1,k)
                do kk=1,j-1
                  xsol(i,j,k)=xsol(i,j,k)-l1m(kk)*xsol(i,kk,k)
                enddo
                !

                !Remontee avec matrice triangulaire sup.
                xsol(i,ny,k)=xsol(i,ny,k)/(ggm(ny)+u1m(ny))
                xsol(i,ny-1,k)=(xsol(i,ny-1,k)-(hhm(ny-1)+u1m(ny-1))*xsol(i,ny,k))/(ggm(ny-1)+u2m(ny-1))
                xsol(i,ny-2,k)=(xsol(i,ny-2,k)-(hhm(ny-2)+u2m(ny-2))*xsol(i,ny-1,k)- &
                (ssm(ny-2)+u1m(ny-2))*xsol(i,ny,k))/(ggm(ny-2)+u3m(ny-2))

                do j=ny-3,1,-1
                  xsol(i,j,k)=(xsol(i,j,k)-hhm(j)*xsol(i,j+1,k)-ssm(j)*xsol(i,j+2,k) &
                  -rrm(j)*xsol(i,j+3,k) &
                  -u3m(j)*xsol(i,ny-2,k) &
                  -u2m(j)*xsol(i,ny-1,k) &
                  -u1m(j)*xsol(i,ny,k) )/ggm(j)
                enddo
              enddo
            enddo

          end subroutine septinv_0
          !********************************************************************
        end module matinv
        !********************************************************************
        subroutine  inttimp (ux1,uy1,uz1,gx,gy,gz,hx,hy,hz,ta1,tb1,tc1,px1,py1,pz1,&
          td1,te1,tf1,ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2)

          USE param
          USE variables
          USE decomp_2d
          use derivY
          use matinv

          implicit none

          integer :: ijk,nxyz,i,j,k,code
          !real(mytype) :: xcst,x,z,r,liss ! pour VERIF inversion matrice, modification module param
          real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
          real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
          real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx,gy,gz
          real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx,hy,hz
          real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1
          real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: td1,te1,tf1
          real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
          real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2
          real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: td2,te2,tf2,di2
          real(mytype),dimension(ysize(1),ysize(3)) :: bc1,bc2,bc3

          td1=zero;te1=zero;tf1=zero
          ta2=zero;tb2=zero;tc2=zero;td2=zero;te2=zero;tf2=zero
          ux2=zero;uy2=zero;uz2=zero

          nxyz=xsize(1)*xsize(2)*xsize(3)

          !!!!!!!!!!!!!
          !! CN2+AB2 !!
          !!!!!!!!!!!!!
          if (nscheme.eq.1) then
            if ((itime.eq.1).and.(ilit.eq.0)) then

              !START WITH EXPLICIT EULER + IMPLICIT CN SCHEMES
              do k=1,xsize(3)
                do j=1,xsize(2)
                  do i=1,xsize(1)
                    !uhat
                    td1(i,j,k)= gdt(itr)*ta1(i,j,k)
                    te1(i,j,k)= gdt(itr)*tb1(i,j,k)
                    tf1(i,j,k)= gdt(itr)*tc1(i,j,k)

                    !save (N+Lxz)(n-1)
                    gx(i,j,k)=ta1(i,j,k)
                    gy(i,j,k)=tb1(i,j,k)
                    gz(i,j,k)=tc1(i,j,k)

                  enddo
                enddo
              enddo

            else !CONTINUE WITH EXPLICIT AB2 + IMPLICIT CN SCHEMES

              do k=1,xsize(3)
                do j=1,xsize(2)
                  do i=1,xsize(1)
                    !uhat
                    td1(i,j,k)= adt(itr)*ta1(i,j,k)+bdt(itr)*gx(i,j,k)-px1(i,j,k)
                    te1(i,j,k)= adt(itr)*tb1(i,j,k)+bdt(itr)*gy(i,j,k)-py1(i,j,k)
                    tf1(i,j,k)= adt(itr)*tc1(i,j,k)+bdt(itr)*gz(i,j,k)-pz1(i,j,k)

                    !save (N+Lxz)(n-1)
                    gx(i,j,k)=ta1(i,j,k)
                    gy(i,j,k)=tb1(i,j,k)
                    gz(i,j,k)=tc1(i,j,k)
                  enddo
                enddo
              enddo

            endif
          endif

          !!!!!!!!!!!!!!!
          !!! CN2+AB3 !!!
          !!!!!!!!!!!!!!!
          if (nscheme==2) then
            if ((itime.eq.1).and.(ilit.eq.0)) then

              if (nrank==0) print *,'start AB3 with Euler',itime
              !START WITH EXPLICIT EULER + IMPLICIT CN SCHEMES

              do k=1,xsize(3)
                do j=1,xsize(2)
                  do i=1,xsize(1)
                    !uhat
                    td1(i,j,k)= dt*ta1(i,j,k)
                    te1(i,j,k)= dt*tb1(i,j,k)
                    tf1(i,j,k)= dt*tc1(i,j,k)

                    !save (N+Lxz)(n-1)
                    gx(i,j,k)=ta1(i,j,k)
                    gy(i,j,k)=tb1(i,j,k)
                    gz(i,j,k)=tc1(i,j,k)

                  enddo
                enddo
              enddo

            else if ((itime.eq.2).and.(ilit.eq.0)) then

              if (nrank==0) print *,'then with AB2',itime
                !THEN AB2 + IMPLICIT CN SCHEMES
                do k=1,xsize(3)
                  do j=1,xsize(2)
                    do i=1,xsize(1)
                      !uhat
                      td1(i,j,k)= 1.5*dt*ta1(i,j,k)-0.5*dt*gx(i,j,k)-px1(i,j,k)
                      te1(i,j,k)= 1.5*dt*tb1(i,j,k)-0.5*dt*gy(i,j,k)-py1(i,j,k)
                      tf1(i,j,k)= 1.5*dt*tc1(i,j,k)-0.5*dt*gz(i,j,k)-pz1(i,j,k)

                      !save (N+Lxz)(n-1)&(n-2)
                      hx(i,j,k)=gx(i,j,k)
                      hy(i,j,k)=gy(i,j,k)
                      hz(i,j,k)=gz(i,j,k)
                      gx(i,j,k)=ta1(i,j,k)
                      gy(i,j,k)=tb1(i,j,k)
                      gz(i,j,k)=tc1(i,j,k)
                    enddo
                  enddo
                enddo

              else !FINALLY EXPLICIT AB3 + IMPLICIT CN SCHEMES

                do k=1,xsize(3)
                  do j=1,xsize(2)
                    do i=1,xsize(1)
                      !uhat
                      td1(i,j,k)= adt(itr)*ta1(i,j,k)+bdt(itr)*gx(i,j,k) &
                      +cdt(itr)*hx(i,j,k)-px1(i,j,k)
                      te1(i,j,k)= adt(itr)*tb1(i,j,k)+bdt(itr)*gy(i,j,k) &
                      +cdt(itr)*hy(i,j,k)-py1(i,j,k)
                      tf1(i,j,k)= adt(itr)*tc1(i,j,k)+bdt(itr)*gz(i,j,k) &
                      +cdt(itr)*hz(i,j,k)-pz1(i,j,k)

                      !save (N+Lxz)(n-1)
                      hx(i,j,k)=gx(i,j,k)
                      hy(i,j,k)=gy(i,j,k)
                      hz(i,j,k)=gz(i,j,k)
                      gx(i,j,k)=ta1(i,j,k)
                      gy(i,j,k)=tb1(i,j,k)
                      gz(i,j,k)=tc1(i,j,k)
                    enddo
                  enddo
                enddo

              endif
            endif

            !Y-PENCIL FOR MATRIX INVERSION
            call transpose_x_to_y(ux1,ux2)
            call transpose_x_to_y(uy1,uy2)
            call transpose_x_to_y(uz1,uz2)
            call transpose_x_to_y(td1,ta2)
            call transpose_x_to_y(te1,tb2)
            call transpose_x_to_y(tf1,tc2)

            bc1(:,:)=ux2(:,ny-1,:)
            bc2(:,:)=uy2(:,ny-1,:)
            bc3(:,:)=uz2(:,ny-1,:)


            !ta2: A.uhat
            !td2:(A+xcstB).un
            if (ncly.eq.0) then
              call multmatrix7(td2,ta2,ux2)
              call multmatrix7(te2,tb2,uy2)
              call multmatrix7(tf2,tc2,uz2)
            elseif (ncly.eq.1) then
            elseif (ncly.eq.2) then
              call multmatrix7(td2,ta2,ux2)
              call multmatrix7(te2,tb2,uy2)
              call multmatrix7(tf2,tc2,uz2)
            endif

            !SECOND MEMBRE COMPLET BBB=A uhat+(A+xcst.B)u^n
            !in  ta2,tb2,tc2
            ta2=td2+ta2
            tb2=te2+tb2
            tc2=tf2+tc2
            !
            ! CONDITIONS AUX LIMITES
            ta2(:,ny,:)=bc1(:,:)
            tb2(:,ny,:)=bc2(:,:)
            tc2(:,ny,:)=bc3(:,:)

            !Inversion systeme lineaire Mx=b: (A-xcst.B)u^n+1=uhat+(A+xcst.B)u^n
            ux2=zero; uy2=zero; uz2=zero;
            if (ncly.eq.0) then
              call septinv(ux2,ta2,ggm0,hhm0,ssm0,rrm0,vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ysize(1),ysize(2),ysize(3))
              call septinv(uy2,tb2,ggm0,hhm0,ssm0,rrm0,vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ysize(1),ysize(2),ysize(3))
              call septinv(uz2,tc2,ggm0,hhm0,ssm0,rrm0,vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ysize(1),ysize(2),ysize(3))
            elseif (ncly.eq.1) then
              call septinv(ux2,ta2,ggm11,hhm11,ssm11,rrm11,vvm11,wwm11,zzm11,ysize(1),ysize(2),ysize(3))
              call septinv(uy2,tb2,ggm10,hhm10,ssm10,rrm10,vvm10,wwm10,zzm10,ysize(1),ysize(2),ysize(3))
              call septinv(uz2,tc2,ggm11,hhm11,ssm11,rrm11,vvm11,wwm11,zzm11,ysize(1),ysize(2),ysize(3))
            elseif (ncly.eq.2) then
              call septinv(ux2,ta2,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ysize(1),ysize(2),ysize(3))
              call septinv(uy2,tb2,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ysize(1),ysize(2),ysize(3))
              call septinv(uz2,tc2,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ysize(1),ysize(2),ysize(3))
            endif

            call transpose_y_to_x(ux2,ux1)
            call transpose_y_to_x(uy2,uy1)
            call transpose_y_to_x(uz2,uz1)

#ifdef DEBUG
            if (.false.) then
              !VERIFICATION INVERSION CORRECTE, CAS NCLY=1
              !Mx-b=0 ?
              ! ux
              if (nrank.eq.0) print *,'Check matrix inversion for ux - NCLY=1'
              call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
              td2(:,1,:)=ux2(:,1,:)+alsajy*(2.*ux2(:,2,:)) - xcst*td2(:,1,:)*pp2y(1)
              do j=2,ysize(2)-1
                td2(:,j,:)=alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) - xcst*td2(:,j,:)*pp2y(j)
              enddo
              td2(:,ysize(2),:)=ux2(:,ysize(2),:)+alsajy*(2.*ux2(:,ysize(2)-1,:)) - xcst*td2(:,ysize(2),:)*pp2y(ysize(2))
              td2=td2-ta2
              call transpose_y_to_x(td2,td1)
              call  test_scalar_min_max(td1)
              ! uy
              if (nrank.eq.0) print *,'Check matrix inversion for uy - NCLY=1'
              call deryy(td2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
              td2(:,1,:)=uy2(:,1,:)+0.*alsajy*(uy2(:,2,:)+uy2(:,ysize(2),:)) - xcst*td2(:,j,:)*pp2y(1)
              do j=2,ysize(2)-1
                td2(:,j,:)=alsajy*uy2(:,j-1,:) + uy2(:,j,:) + alsajy*uy2(:,j+1,:) - xcst*td2(:,j,:)*pp2y(j)
              enddo
              td2(:,ysize(2),:)=uy2(:,ysize(2),:)+0.*alsajy*(uy2(:,ysize(2)-1,:)+uy2(:,1,:)) - xcst*td2(:,ysize(2),:)*pp2y(ysize(2))
              td2=td2-tb2
              call transpose_y_to_x(td2,td1)
              call  test_scalar_min_max(td1)
              ! uz
              if (nrank.eq.0) print *,'Check matrix inversion for uz - NCLY=1'
              call deryy(td2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
              td2(:,1,:)=uz2(:,1,:)+alsajy*(2.*uz2(:,2,:)) - xcst*td2(:,1,:)*pp2y(1)
              do j=2,ysize(2)-1
                td2(:,j,:)=alsajy*uz2(:,j-1,:) + uz2(:,j,:) + alsajy*uz2(:,j+1,:) - xcst*td2(:,j,:)*pp2y(j)
              enddo
              td2(:,ysize(2),:)=uz2(:,ysize(2),:)+alsajy*(2.*uz2(:,ysize(2)-1,:)) - xcst*td2(:,ysize(2),:)*pp2y(ysize(2))
              td2=td2-tc2
              call transpose_y_to_x(td2,td1)
              call  test_scalar_min_max(td1)
            endif
#endif

            ! IF ab2 OR ab3
            if ((nscheme.eq.1).or.(nscheme.eq.2)) then
              ! IF NOT (first time step with no restart)
              if ( (ilit.eq.1).or.(itime.gt.1) ) then
                ux1=ux1+px1
                uy1=uy1+py1
                uz1=uz1+pz1
              endif
            endif

            return
          end subroutine inttimp
          !********************************************************************
          !PREMULT FOR INTTIMP WITH SEPTADIAG
          subroutine multmatrix7(td2,ta2,ux2)

            USE param
            USE variables
            USE derivY
            USE decomp_2d

            implicit none

            integer :: i,j,k,code
            !real(mytype) :: xcst ! modification module param
            real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: ux2
            real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: td2,ta2
            real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: di2

            !xcst= xnu*gdt(itr)*0.5

            !A.uhat
            if (.true.) then

              if (istret.ne.0) then
                do j=1,ysize(2)
                  ta2(:,j,:)=ta2(:,j,:)/pp2y(j)
                enddo
              endif

              if (ncly.eq.0) then

                td2(:,1,:) = alsajy*ta2(:,2,:) + ta2(:,1,:) + alsajy*ta2(:,ysize(2),:)
                do j=2,ysize(2)-1
                  td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
                enddo
                td2(:,ysize(2),:) = alsajy*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:) + alsajy*ta2(:,1,:)
                ta2=td2

              elseif (ncly.eq.1) then


              elseif (ncly.eq.2) then

                td2(:,1,:) = 0.
                td2(:,2,:) = alsa2y*ta2(:,1,:) + ta2(:,2,:) + alsa2y*ta2(:,3,:)
                td2(:,3,:) = alsa3y*ta2(:,2,:) + ta2(:,3,:) + alsa3y*ta2(:,4,:)
                do j=4,ysize(2)-3
                  td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
                enddo
                td2(:,ysize(2)-2,:) = alsaty*ta2(:,ysize(2)-3,:) + ta2(:,ysize(2)-2,:) + alsaty*ta2(:,ysize(2)-1,:)
                td2(:,ysize(2)-1,:) = alsamy*ta2(:,ysize(2)-2,:) + ta2(:,ysize(2)-1,:) + alsamy*ta2(:,ysize(2),:)
                td2(:,ysize(2),:) = 0.
                ta2=td2

              endif

            else

              do k=1,ysize(3)
                do i=1,ysize(1)
                  td2(i,1,k)=zero!
                  td2(i,2,k)=alsa2y*ta2(i,1,k)+ta2(i,2,k)+alsa2y*ta2(i,3,k)
                  td2(i,3,k)=alsa3y*ta2(i,2,k)+ta2(i,3,k)+alsa3y*ta2(i,4,k)
                  do j=4,ysize(2)-3
                    td2(i,j,k)=alsajy*ta2(i,j-1,k)+ta2(i,j,k)+alsajy*ta2(i,j+1,k)
                  enddo
                  td2(i,ysize(2)-2,k)=alsaty*ta2(i,ysize(2)-3,k)+ta2(i,ysize(2)-2,k)+alsaty*ta2(i,ysize(2)-1,k)
                  td2(i,ysize(2)-1,k)=alsamy*ta2(i,ysize(2)-2,k)+ta2(i,ysize(2)-1,k)+alsamy*ta2(i,ysize(2),k)
                  td2(i,ysize(2),k)=zero!
                enddo
              enddo

              do k=1,ysize(3)
                do j=1,ysize(2)
                  do i=1,ysize(1)
                    ta2(i,j,k)=td2(i,j,k)
                  enddo
                enddo
              enddo

            endif

            !(A+nu*dt.B).un

            if (.true.) then

              if (ncly.eq.0) then

                call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)

              elseif (ncly.eq.1) then



              elseif (ncly.eq.2) then

                call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)

              endif

              if (istret.ne.0) then
                do j=1,ysize(2)
                  ux2(:,j,:)=ux2(:,j,:)/pp2y(j)
                enddo
              endif

              if (ncly.eq.0) then

                td2(:,1,:) = alsajy*ux2(:,2,:) + ux2(:,1,:) + alsajy*ux2(:,ysize(2),:) &
                + xcst*td2(:,1,:)
                do j=2,ysize(2)-1
                  td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                  + xcst*td2(:,j,:)
                enddo
                td2(:,ysize(2),:) = alsajy*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) + alsajy*ux2(:,1,:) &
                + xcst*td2(:,ysize(2),:)

              elseif (ncly.eq.1) then



              elseif (ncly.eq.2) then

                td2(:,1,:) = 0.
                td2(:,2,:) = alsa2y*ux2(:,1,:) + ux2(:,2,:) + alsa2y*ux2(:,3,:) &
                + xcst*td2(:,2,:)
                td2(:,3,:) = alsa3y*ux2(:,2,:) + ux2(:,3,:) + alsa3y*ux2(:,4,:) &
                + xcst*td2(:,3,:)
                do j=4,ysize(2)-3
                  td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                  + xcst*td2(:,j,:)
                enddo
                td2(:,ysize(2)-2,:) = alsaty*ux2(:,ysize(2)-3,:) + ux2(:,ysize(2)-2,:) + alsaty*ux2(:,ysize(2)-1,:) &
                + xcst*td2(:,ysize(2)-2,:)
                td2(:,ysize(2)-1,:) = alsamy*ux2(:,ysize(2)-2,:) + ux2(:,ysize(2)-1,:) + alsamy*ux2(:,ysize(2),:) &
                + xcst*td2(:,ysize(2)-1,:)
                td2(:,ysize(2),:) = 0.

              endif

            else

              if(istret==0) then

                do k=1,ysize(3)
                  do i=1,ysize(1)

                    !for ux
                    td2(i,1,k)= 0.!
                    td2(i,2,k)=(alsa2y+xcst*as2y)*ux2(i,1,k)+(1.-xcst*2.*as2y)*ux2(i,2,k)+(alsa2y+xcst*as2y)*ux2(i,3,k)
                    td2(i,3,k)=xcst*bs3y*ux2(i,1,k)+(alsa3y+xcst*as3y)*ux2(i,2,k) &
                    + (1.-xcst*2.*(as3y+bs3y))*ux2(i,3,k) + (alsa3y+xcst*as3y)*ux2(i,4,k) &
                    +xcst*bs3y*ux2(i,5,k)
                    do j=4,ysize(2)-3
                      td2(i,j,k)=xcst*csjy*ux2(i,j-3,k)+xcst*bsjy*ux2(i,j-2,k)+(alsajy+xcst*asjy)*ux2(i,j-1,k) &
                      + (1.-xcst*2.*(asjy+bsjy+csjy))*ux2(i,j,k) + (alsajy+xcst*asjy)*ux2(i,j+1,k) &
                      +xcst*bsjy*ux2(i,j+2,k)+xcst*csjy*ux2(i,j+3,k)
                    enddo
                    td2(i,ysize(2)-2,k)=xcst*bsty*ux2(i,ysize(2)-4,k)+(alsaty+xcst*asty)*ux2(i,ysize(2)-3,k) &
                    + (1.-xcst*2.*(asty+bsty))*ux2(i,ysize(2)-2,k) + (alsaty+xcst*asty)*ux2(i,ysize(2)-1,k) &
                    +xcst*bsty*ux2(i,ysize(2),k)
                    td2(i,ysize(2)-1,k)=(alsamy+xcst*asmy)*ux2(i,ysize(2)-2,k)+(1.-xcst*2.*asmy)*ux2(i,ysize(2)-1,k) &
                    +(alsamy+xcst*asmy)*ux2(i,ysize(2),k)
                    td2(i,ysize(2),k)=zero!xcst*csny*ux2(i,ysize(2)-2,k)+(alsany+xcst*bsny)*ux2(i,ysize(2)-1,k) &
                    !+(1.+xcst*asny)*ux2(i,ysize(2),k)

                  enddo
                enddo

              else

                do k=1,ysize(3)
                  do i=1,ysize(1)

                    !for ux
                    td2(i,1,k)= 0.!(1.+xcst*as1y)*ux2(i,1,k)+(alsa1y+xcst*bs1y)*ux2(i,2,k) &
                    !+ xcst*cs1y*ux2(i,3,k)
                    td2(i,2,k)=(alsa2y+xcst*as2y*pp2y(2))*ux2(i,1,k)+(1.-xcst*2.*as2y*pp2y(2))*ux2(i,2,k)&
                    +(alsa2y+xcst*as2y*pp2y(2))*ux2(i,3,k)
                    td2(i,3,k)=xcst*bs3y*ux2(i,1,k)*pp2y(3)+(alsa3y+xcst*as3y*pp2y(3))*ux2(i,2,k) &
                    + (1.-xcst*2.*(as3y+bs3y)*pp2y(3))*ux2(i,3,k) + (alsa3y+xcst*as3y*pp2y(3))*ux2(i,4,k) &
                    +xcst*bs3y*ux2(i,5,k)*pp2y(3)
                    do j=4,ysize(2)-3
                      td2(i,j,k)=xcst*csjy*ux2(i,j-3,k)*pp2y(j)+xcst*bsjy*ux2(i,j-2,k)*pp2y(j) &
                      +(alsajy+xcst*asjy*pp2y(j))*ux2(i,j-1,k) &
                      + (1.-xcst*2.*(asjy+bsjy+csjy)*pp2y(j))*ux2(i,j,k) + (alsajy+xcst*asjy*pp2y(j))*ux2(i,j+1,k) &
                      +xcst*bsjy*ux2(i,j+2,k)*pp2y(j)+xcst*csjy*ux2(i,j+3,k)*pp2y(j)
                    enddo
                    td2(i,ysize(2)-2,k)=xcst*bsty*ux2(i,ysize(2)-4,k)*pp2y(ysize(2)-2)+ &
                    (alsaty+xcst*asty*pp2y(ysize(2)-2))*ux2(i,ysize(2)-3,k) &
                    + (1.-xcst*2.*(asty+bsty)*pp2y(ysize(2)-2))*ux2(i,ysize(2)-2,k) &
                    + (alsaty+xcst*asty*pp2y(ysize(2)-2))*ux2(i,ysize(2)-1,k) &
                    +xcst*bsty*ux2(i,ysize(2),k)*pp2y(ysize(2)-2)
                    td2(i,ysize(2)-1,k)=(alsamy+xcst*asmy*pp2y(ysize(2)-1))*ux2(i,ysize(2)-2,k)+ &
                    (1.-xcst*2.*asmy*pp2y(ysize(2)-1))*ux2(i,ysize(2)-1,k) &
                    +(alsamy+xcst*asmy*pp2y(ysize(2)-1))*ux2(i,ysize(2),k)
                    td2(i,ysize(2),k)=zero!xcst*csny*ux2(i,ysize(2)-2,k)+(alsany+xcst*bsny)*ux2(i,ysize(2)-1,k) &
                    !+(1.+xcst*asny)*ux2(i,ysize(2),k)

                  enddo
                enddo

              endif

            endif


          end subroutine multmatrix7
          !*******************************************************************
          subroutine implicit_schemes()

            USE param
            USE derivY
            USE variables
            USE var
            use ludecomp

            implicit none

            integer  :: i,j,k

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! MATRICE M2 TIME IMPLICIT !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!!
            !FOR VELOCITY  !!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!
            xcst = xnu*dt*0.5
            !!!!!!!!!!!!!!!!!!!!!!
            !CAS SEPTADIAGONALE !!
            !!!!!!!!!!!!!!!!!!!!!!

            !!! NCL = 2, dirichlet imposé, fonction nulle à la paroi
            !
            !DIAG
            aam(1     )=as1y
            aam(ny    )=asny
            aam(2     )=-two*as2y
            aam(ny-1  )=-two*asmy
            aam(3     )=-two*(as3y+bs3y)
            aam(ny-2  )=-two*(asty+bsty)
            aam(4:ny-3)=-two*(asjy+bsjy+csjy)
            if (istret==0) then
              aam = one-xcst*aam
            else
              aam = one/pp2y -xcst*aam
            endif
            !CL sur aam
            aam(1 )=one
            aam(ny)=one
            !
            !DIAG SUP 1
            bbm(1     )=bs1y
            bbm(ny    )=bsny
            bbm(2     )=as2y
            bbm(ny-1  )=asmy
            bbm(3     )=as3y
            bbm(ny-2  )=asty
            bbm(4:ny-3)=asjy
            bbm = -xcst*bbm
            if (istret==0) then
              bbm(2     )=bbm(2     )+alsa2y
              bbm(ny-1  )=bbm(ny-1  )+alsamy
              bbm(3     )=bbm(3     )+alsa3y
              bbm(ny-2  )=bbm(ny-2  )+alsaty
              bbm(4:ny-3)=bbm(4:ny-3)+alsajy
            else
              bbm(2     )=bbm(2     )+alsa2y/pp2y(3)
              bbm(ny-1  )=bbm(ny-1  )+alsamy/pp2y(ny)
              bbm(3     )=bbm(3     )+alsa3y/pp2y(4)
              bbm(ny-2  )=bbm(ny-2  )+alsaty/pp2y(ny-1)
              bbm(4:ny-3)=bbm(4:ny-3)+alsajy/pp2y(5:ny-2)
            endif
            !CL sur bbm
            bbm(1 )=zero
            bbm(ny)=zero
            !
            !DIAG SUP 2
            ccm(1     )=cs1y
            ccm(ny    )=csny
            ccm(2     )=zero
            ccm(ny-1  )=zero
            ccm(3     )=bs3y
            ccm(ny-2  )=bsty
            ccm(4:ny-3)=bsjy
            ccm = -xcst*ccm
            !CL sur ccm
            ccm(1 )=zero
            ccm(ny)=zero
            !
            !DIAG SUP 3
            rrm(1     )=ds1y
            rrm(ny    )=dsny
            rrm(2     )=zero
            rrm(ny-1  )=zero
            rrm(3     )=zero
            rrm(ny-2  )=zero
            rrm(4:ny-3)=csjy
            rrm = -xcst*rrm
            !CL sur rrm
            rrm(1 )=zero
            rrm(ny)=zero
            !
            !DIAG INF 1
            if (istret==0) then
              ddm=bbm
            else
              ddm(1     )=bs1y
              ddm(ny    )=bsny
              ddm(2     )=as2y
              ddm(ny-1  )=asmy
              ddm(3     )=as3y
              ddm(ny-2  )=asty
              ddm(4:ny-3)=asjy
              ddm = -xcst*ddm
              ddm(2     )=ddm(2     )+alsa2y/pp2y(1)
              ddm(ny-1  )=ddm(ny-1  )+alsamy/pp2y(ny-2)
              ddm(3     )=ddm(3     )+alsa3y/pp2y(2)
              ddm(ny-2  )=ddm(ny-2  )+alsaty/pp2y(ny-3)
              ddm(4:ny-3)=ddm(4:ny-3)+alsajy/pp2y(3:ny-4)
              !CL sur ddm
              ddm(1 )=zero
              ddm(ny)=zero
            endif
            !
            !DIAG INF 2
            eem=ccm
            !
            !DIAG INF 3
            qqm=rrm

            !!! NCL = 1, npaire=0, dirichlet imposé, fonction impaire
            !
            ! DIAG
            aam10(1     )=zero
            aam10(ny    )=zero
            aam10(2     )=-two*asjy-3.*bsjy-two*csjy
            aam10(ny-1  )=aam10(2)
            aam10(3:ny-2)=-two*(asjy+bsjy+csjy)
            if (istret==0) then
              aam10 = 1. - xcst*aam10
            else
              aam10 = 1./pp2y - xcst*aam10
            endif
            !
            !DIAG SUP 1
            bbm10(1     )=zero
            bbm10(ny    )=zero
            bbm10(2     )=asjy-csjy
            bbm10(ny-1  )=asjy
            bbm10(3     )=asjy
            bbm10(ny-2  )=asjy-csjy
            bbm10(4:ny-3)=asjy
            if (istret==0) then
              bbm10 = alsajy - xcst*bbm10
            else
              bbm10(2:ny-1) = alsajy/pp2y(3:ny) - xcst*bbm10(2:ny-1)
            endif
            !CL sur bbm10
            bbm10(1 )=zero
            bbm10(ny)=zero
            !
            !DIAG SUP 2
            ccm10(1     )=zero
            ccm10(ny    )=zero
            ccm10(2     )=bsjy
            ccm10(ny-1  )=zero
            ccm10(3     )=bsjy
            ccm10(ny-2  )=bsjy
            ccm10(4:ny-3)=bsjy
            ccm10 = -xcst*ccm10
            !
            !DIAG SUP 3
            rrm10(1     )=zero
            rrm10(ny    )=zero
            rrm10(2     )=csjy
            rrm10(ny-1  )=zero
            rrm10(3     )=csjy
            rrm10(ny-2  )=zero
            rrm10(4:ny-3)=csjy
            rrm10 = -xcst*rrm10
            !
            !DIAG INF 1
            ddm10(1     )=zero
            ddm10(ny    )=zero
            ddm10(2     )=asjy
            ddm10(ny-1  )=asjy-csjy
            ddm10(3     )=asjy-csjy
            ddm10(ny-2  )=asjy
            ddm10(4:ny-3)=asjy
            if (istret==0) then
              ddm10 = alsajy - xcst*ddm10
            else
              ddm10(2:ny-1) = alsajy/pp2y(1:ny-2) - xcst*ddm10(2:ny-1)
            endif
            !CL sur ddm10
            ddm10(1 )=zero
            ddm10(ny)=zero
            !
            !DIAG INF 2
            eem10(1     )=zero
            eem10(ny    )=zero
            eem10(2     )=zero
            eem10(ny-1  )=bsjy
            eem10(3     )=bsjy
            eem10(ny-2  )=bsjy
            eem10(4:ny-3)=bsjy
            eem10 = -xcst*eem10
            !
            !DIAG INF 3
            qqm10(1     )=zero
            qqm10(ny    )=zero
            qqm10(2     )=zero
            qqm10(ny-1  )=csjy
            qqm10(3     )=zero
            qqm10(ny-2  )=csjy
            qqm10(4:ny-3)=csjy
            qqm10 = -xcst*qqm10

            !!! NCL = 1, npaire=1, neumann imposé, fonction paire
            !
            ! DIAG
            aam11(1     )=-two*(asjy+bsjy+csjy)
            aam11(ny    )=aam11(1)
            aam11(2     )=-two*asjy-bsjy-two*csjy
            aam11(ny-1  )=aam11(2)
            aam11(3:ny-2)=-two*(asjy+bsjy+csjy)
            if (istret==0) then
              aam11 = 1. - xcst*aam11
            else
              aam11 = 1./pp2y - xcst*aam11
            endif
            !
            !DIAG SUP 1
            bbm11(1     )=two*asjy
            bbm11(ny    )=zero
            bbm11(2     )=asjy+csjy
            bbm11(ny-1  )=asjy
            bbm11(3     )=asjy
            bbm11(ny-2  )=asjy+csjy
            bbm11(4:ny-3)=asjy
            if (istret==0) then
              bbm11 = alsajy - xcst*bbm11
            else
              bbm11(1:ny-1) = alsajy/pp2y(2:ny) - xcst*bbm11(1:ny-1)
            endif
            !CL sur bbm11
            if (istret==0) then
              bbm11(1 )=bbm11(1)+alsajy
            else
              bbm11(1 )=bbm11(1)+alsajy/pp2y(2)
            endif
            bbm11(ny)=zero
            !
            !DIAG SUP 2
            ccm11(1     )=two*bsjy
            ccm11(ny    )=zero
            ccm11(2     )=bsjy
            ccm11(ny-1  )=zero
            ccm11(3     )=bsjy
            ccm11(ny-2  )=bsjy
            ccm11(4:ny-3)=bsjy
            ccm11 = -xcst*ccm11
            !
            !DIAG SUP 3
            rrm11(1     )=two*csjy
            rrm11(ny    )=zero
            rrm11(2     )=csjy
            rrm11(ny-1  )=zero
            rrm11(3     )=csjy
            rrm11(ny-2  )=zero
            rrm11(4:ny-3)=csjy
            rrm11 = -xcst*rrm11
            !
            !DIAG INF 1
            ddm11(1     )=zero
            ddm11(ny    )=two*asjy
            ddm11(2     )=asjy
            ddm11(ny-1  )=asjy+csjy
            ddm11(3     )=asjy+csjy
            ddm11(ny-2  )=asjy
            ddm11(4:ny-3)=asjy
            if (istret==0) then
              ddm11 = alsajy - xcst*ddm11
            else
              ddm11(2:ny) = alsajy/pp2y(1:ny-1) - xcst*ddm11(2:ny)
            endif
            !CL sur ddm11
            ddm11(1 )=zero
            if (istret==0) then
              ddm11(ny)=ddm11(ny)+alsajy!a1
            else
              ddm11(ny)=ddm11(ny)+alsajy/pp2y(ny-1)!a1
            endif
            !
            !DIAG INF 2
            eem11(1     )=zero
            eem11(ny    )=two*bsjy
            eem11(2     )=zero
            eem11(ny-1  )=bsjy
            eem11(3     )=bsjy
            eem11(ny-2  )=bsjy
            eem11(4:ny-3)=bsjy
            eem11 = -xcst*eem11
            !
            !DIAG INF 3
            qqm11(1     )=zero
            qqm11(ny    )=two*csjy
            qqm11(2     )=zero
            qqm11(ny-1  )=csjy
            qqm11(3     )=zero
            qqm11(ny-2  )=csjy
            qqm11(4:ny-3)=csjy
            qqm11 = -xcst*qqm11

            !!! NXL = 0
            !DIAG
            if (istret==0) then
              aam0 = 1.-xcst*(-two*(asjy+bsjy+csjy))
            else
              aam0 = 1./pp2y-xcst*(-two*(asjy+bsjy+csjy))
            endif
            !
            !DIAG SUP 1
            if (istret==0) then
              bbm0 = alsajy-xcst*asjy
            else
              bbm0(1:ny-1) = alsajy/pp2y(2:ny) -xcst*asjy
              bbm0(ny) = alsajy/pp2y(1) -xcst*asjy
            endif
            !
            !DIAG SUP 2
            ccm0 = -xcst*bsjy
            !
            !DIAG SUP 3
            rrm0 = -xcst*csjy
            !
            !DIAG INF 1
            if (istret==0) then
              ddm0=bbm0
            else
              ddm0(1) = alsajy/pp2y(ny) -xcst*asjy
              ddm0(2:ny) = alsajy/pp2y(1:ny-1) -xcst*asjy
            endif
            !
            !DIAG INF 2
            eem0=ccm0
            !
            !DIAG INF 3
            qqm0=rrm0

            call ludecomp7(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,&
            vvm,wwm,zzm,ny)
            call ludecomp7(aam10,bbm10,ccm10,ddm10,eem10,qqm10,ggm10,hhm10,ssm10,rrm10,&
            vvm10,wwm10,zzm10,ny)
            call ludecomp7(aam11,bbm11,ccm11,ddm11,eem11,qqm11,ggm11,hhm11,ssm11,rrm11,&
            vvm11,wwm11,zzm11,ny)
            call ludecomp7(aam0,bbm0,ccm0,ddm0,eem0,qqm0,ggm0,hhm0,ssm0,rrm0,&
            vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ny)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! FIN MATRICE M2 TIME IMPLICIT    !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          end subroutine implicit_schemes
          !************************************************************
          subroutine scalarimp(ux1,uy1,uz1,phi1,phis1,phiss1,di1,tg1,th1,ti1,td1,&
            uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3)

            !IMPLICIT TIME INTEGRATION FOR D2/DY2
            USE param
            USE variables
            USE decomp_2d
            USE derivY
            USE MPI
            use matinv
#ifdef my_mod_solide
            use conjugate_ht, only : temp_bot,temp_top,ny_sol_bot,cl_bot,cl_top,update_temp_solide
#endif

            implicit none

            real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
            real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi1,phis1,phiss1
            real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: di1,tg1,th1,ti1,td1
            !real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phisauv

            real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,phi2,ta2,tb2,tc2,td2,di2

            real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3,di3,ta3,tb3
            integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz,code

            !parametres CL
            real(mytype) :: x,y,z,r,lambda,phislbda,adiab,tjet,liss

#ifdef my_mod_solide
            real(mytype),dimension(ysize(1),2,ysize(3)) :: mytmptemp
#endif

            tg1=zero;th1=zero;ti1=zero;td1=zero
            ta2=zero;tb2=zero;tc2=zero;td2=zero
            ta3=zero;tb3=zero

            nvect1=xsize(1)*xsize(2)*xsize(3)
            nvect2=ysize(1)*ysize(2)*ysize(3)
            nvect3=zsize(1)*zsize(2)*zsize(3)

            !X PENCILS
            tg1=ux1*phi1
            call derx (th1,tg1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! npaire=0 pour phi
            !call derx (th1,phi1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! npaire=1
            !call derx (th1,phi1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! npaire=0
            !call derxxt (tg1,phi1,di1,sx,sfxpt,ssxpt,swxpt,xsize(1),xsize(2),xsize(3),1) ! npaire=1
            call derxxt (tg1,phi1,di1,sx,sfxt,ssxt,swxt,xsize(1),xsize(2),xsize(3),0) ! npaire=0
            call transpose_x_to_y(phi1,phi2)
            call transpose_x_to_y(uy1,uy2)
            call transpose_x_to_y(uz1,uz2)

            !Y PENCILS
            tc2=phi2*uy2
            call dery (tb2,tc2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! npaire=0 pour phi
            !call dery (tb2,phi2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! npaire=1
            !call dery (tb2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! npaire=0

            if (istret.ne.0) then

              !   call dery (tc2,phi2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! npaire=1
              call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! npaire=0

              do k=1,ysize(3)
                do j=1,ysize(2)
                  do i=1,ysize(1)
                    ta2(i,j,k)=-pp4y(j)*tc2(i,j,k)
                  enddo
                enddo
              enddo

            endif

            call transpose_y_to_z(phi2,phi3)
            call transpose_y_to_z(uz2,uz3)

            !Z PENCILS
            ta3=uz3*phi3
            call derz (tb3,ta3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) ! npaire=0 pour phi
            !call derz (tb3,phi3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) ! npaire=1
            !call derz (tb3,phi3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! npaire=0
            !call derzzt (ta3,phi3,di3,sz,sfzpt,sszpt,swzpt,zsize(1),zsize(2),zsize(3),1) ! npaire=1
            call derzzt (ta3,phi3,di3,sz,sfzt,sszt,swzt,zsize(1),zsize(2),zsize(3),0) ! npaire=0

            call transpose_z_to_y(ta3,tc2)
            call transpose_z_to_y(tb3,td2)

            !Y PENCILS ADD TERMS
            do k=1,ysize(3)
              do j=1,ysize(2)
                do i=1,ysize(1)
                  tc2(i,j,k)=tc2(i,j,k)+ta2(i,j,k) !SECOND DERIVATIVE
                  !         td2(i,j,k)=uz2(i,j,k)*td2(i,j,k)+uy2(i,j,k)*tb2(i,j,k) !FIRST DERIVATIVE
                  td2(i,j,k)=td2(i,j,k)+tb2(i,j,k) !FIRST DERIVATIVE
                enddo
              enddo
            enddo

            call transpose_y_to_x(tc2,ti1)
            call transpose_y_to_x(td2,td1)


            !X PENCILS ADD TERMS
            do k=1,xsize(3)
              do j=1,xsize(2)
                do i=1,xsize(1)
                  tg1(i,j,k)= tg1(i,j,k)+ti1(i,j,k) !SECOND DERIVATIVE
                  !         th1(i,j,k)= ux1(i,j,k)*th1(i,j,k)+td1(i,j,k) !FIRST DERIVATIVE
                  th1(i,j,k)= th1(i,j,k)+td1(i,j,k) !FIRST DERIVATIVE
                enddo
              enddo
            enddo

            do k=1,xsize(3)
              do j=1,xsize(2)
                do i=1,xsize(1)
                  tg1(i,j,k)= (xnu/nsc(1))*tg1(i,j,k)-th1(i,j,k)
                enddo
              enddo
            enddo

            !!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! Kasagi Source term !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            tg1=tg1+ux1*xnu/nsc(1)

            ! TIME ADVANCEMENT EXPLICIT AB + IMPLICIT CN2 (d2/dy2)
            nxyz=xsize(1)*xsize(2)*xsize(3)

            !!!!!!!!!!!!!
            !! CN2+AB2 !!
            !!!!!!!!!!!!!
            if (nscheme.eq.1) then
              if (itime.eq.1.and.ilit.eq.0) then !start with EXPLICIT EULER + CN2

                do k=1,xsize(3)
                  do j=1,xsize(2)
                    do i=1,xsize(1)
                      td1(i,j,k)= gdt(itr)*tg1(i,j,k)
                      phis1(i,j,k)= tg1(i,j,k)
                    enddo
                  enddo
                enddo

              else !CONTINUE WITH EXPLICIT AB2 + IMPLICIT CN SCHEMES

                do k=1,xsize(3)
                  do j=1,xsize(2)
                    do i=1,xsize(1)
                      td1(i,j,k)= adt(itr)*tg1(i,j,k)+bdt(itr)*phis1(i,j,k)
                      phis1(i,j,k)= tg1(i,j,k)
                    enddo
                  enddo
                enddo

              endif

            endif


            !!!!!!!!!!!!!!!!!!!!
            !CN2+AB3         !!!
            !!!!!!!!!!!!!!!!!!!!
            if (nscheme==2) then
              if ((itime.eq.1).and.(ilit.eq.0)) then
                if (nrank==0) print *,'Temperature start with Euler',itime
                !START WITH EXPLICIT EULER + IMPLICIT CN SCHEMES

                do k=1,xsize(3)
                  do j=1,xsize(2)
                    do i=1,xsize(1)
                      td1(i,j,k)= dt*tg1(i,j,k)
                      phis1(i,j,k)= tg1(i,j,k)
                    enddo
                  enddo
                enddo

              else !CONTINUE WITH EXPLICIT AB2 + IMPLICIT CN SCHEMES
                if  ((itime.eq.2).and.(ilit.eq.0)) then
                  if (nrank==0) print *,'then with AB2',itime

                    do k=1,xsize(3)
                      do j=1,xsize(2)
                        do i=1,xsize(1)
                          td1(i,j,k)= 1.5*dt*tg1(i,j,k)-0.5*dt*phis1(i,j,k)
                          phiss1(i,j,k)=phis1(i,j,k)
                          phis1(i,j,k)= tg1(i,j,k)
                        enddo
                      enddo
                    enddo

                  else !FINALLY EXPLICIT AB3 + IMPLICIT CN SCHEMES

                    do k=1,xsize(3)
                      do j=1,xsize(2)
                        do i=1,xsize(1)
                          td1(i,j,k)= adt(itr)*tg1(i,j,k)+bdt(itr)*phis1(i,j,k)&
                          +cdt(itr)*phiss1(i,j,k)
                          phiss1(i,j,k)=phis1(i,j,k)
                          phis1(i,j,k)= tg1(i,j,k)
                        enddo
                      enddo
                    enddo

                  endif
                endif
              endif

              !Y-PENCIL FOR MATRIX INVERSION
              call transpose_x_to_y(phi1,phi2)
              call transpose_x_to_y(td1,ta2)

#ifdef my_mod_solide
              mytmptemp(:,1,:)=phi2(:,1,:)
              mytmptemp(:,2,:)=phi2(:,ysize(2),:)
#endif

              !ta2: A.T_hat
              !td2:(A+xcstB).Tn
              call multmatrix7T(td2,ta2,phi2)

              !right hand side
              ta2 = ta2 + td2
              if (ncly==2) then
#ifdef my_mod_solide
                ta2(:,1,:) = 0.5*(mytmptemp(:,1,:)+temp_bot(:,ny_sol_bot+3,:))
                ta2(:,ysize(2),:) = 0.5*(mytmptemp(:,2,:)+temp_top(:,1,:))
                !if (itime.le.50000) then
                !   ta2(:,1,:)=g_0!(itime-1)*ta2(:,1,:)/50000.+g_0*(50001-itime)/50000.
                !   ta2(:,ysize(2),:)=g_n!(itime-1)*ta2(:,ysize(2),:)/50000.+g_n*(50001-itime)/50000.
                !endif
#else
                ta2(:,1       ,:)=g_0
                ta2(:,ysize(2),:)=g_n
#endif
              endif

              !Inversion systeme lineaire Mx=b: (A-xcst.B)u^n+1=uhat+(A+xcst.B)u^n
              phi2=zero
              if (ncly==0) then
                call septinv(phi2,ta2,ggmt0,hhmt0,ssmt0,rrmt0,vvmt0,wwmt0,zzmt0,l1mt,l2mt,l3mt,u1mt,u2mt,u3mt,ysize(1),ysize(2),ysize(3))
              elseif (ncly==1) then
                call septinv(phi2,ta2,ggmt1,hhmt1,ssmt1,rrmt1,vvmt1,wwmt1,zzmt1,ysize(1),ysize(2),ysize(3))
              elseif (ncly==2) then
                call septinv(phi2,ta2,ggmt,hhmt,ssmt,rrmt,vvmt,wwmt,zzmt,ysize(1),ysize(2),ysize(3))
              endif

#ifdef my_mod_solide
              call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! npaire=0
              cl_bot(:,1,:)=tc2(:,1       ,:)
              cl_top(:,1,:)=tc2(:,ysize(2),:)
              !if (itime.le.50000) then
              !   cl_bot(:,1,:)=(itime-1)*cl_bot(:,1,:)/50000.+1.*(50001-itime)/50000.
              !   cl_top(:,1,:)=(itime-1)*cl_top(:,1,:)/50000.-1.*(50001-itime)/50000.
              !endif
              call update_temp_solide()
#endif

              call transpose_y_to_x(phi2,phi1)

#ifdef DEBUG
              if (.false.) then
                !VERIFICATION INVERSION CORRECTE, CAS NCLY=0
                !Mx-b=0 ?
                ! T
                if (nrank.eq.0) print *,'Check matrix inversion for T - NCLY=0'
                call deryy(td2,phi2,di2,sy,sfyt,ssyt,swyt,ysize(1),ysize(2),ysize(3),0)
                td2(:,1,:)=phi2(:,1,:)+alsajy*(phi2(:,2,:)+phi2(:,ysize(2),:)) - xcst*td2(:,1,:)
                do j=2,ysize(2)-1
                  td2(:,j,:)=alsajy*phi2(:,j-1,:) + phi2(:,j,:) + alsajy*phi2(:,j+1,:) - xcst*td2(:,j,:)
                enddo
                td2(:,ysize(2),:)=phi2(:,ysize(2),:)+alsajy*(phi2(:,ysize(2)-1,:)+phi2(:,1,:)) - xcst*td2(:,ysize(2),:)
                td2=td2-ta2
                call transpose_y_to_x(td2,td1)
                call  test_scalar_min_max(td1)
              endif
#endif

            end subroutine scalarimp

            !********************************************************************
            !PREMULT FOR INTTIMP WITH SEPTADIAG
            subroutine multmatrix7T(td2,ta2,ux2)
              !
              !********************************************************************
              USE param
              USE variables
              USE derivY
              USE decomp_2d

              implicit none

              integer :: i,j,k,code
              !real(mytype) :: xcst modification module param
              real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: ux2
              real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: td2,ta2
              real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: di2

              !xcst= xnu*gdt(itr)*0.5

              !A.uhat
              if (.true.) then

                if (istret.ne.0) then
                  do j=1,ysize(2)
                    ta2(:,j,:)=ta2(:,j,:)/pp2y(j)
                  enddo
                endif

                if (ncly.eq.0) then

                  td2(:,1,:) = alsajyt*ta2(:,2,:) + ta2(:,1,:) + alsajyt*ta2(:,ysize(2),:)
                  do j=2,ysize(2)-1
                    td2(:,j,:) = alsajyt*ta2(:,j-1,:) + ta2(:,j,:) + alsajyt*ta2(:,j+1,:)
                  enddo
                  td2(:,ysize(2),:) = alsajyt*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:) + alsajyt*ta2(:,1,:)
                  ta2=td2

                elseif (ncly.eq.1) then

                  !      td2(:,1,:) = 2.*alsajyt*ta2(:,2,:) + ta2(:,1,:) ! npaire=1
                  td2(:,1,:) = ta2(:,1,:) ! npaire=0
                  do j=2,ysize(2)-1
                    td2(:,j,:) = alsajyt*ta2(:,j-1,:) + ta2(:,j,:) + alsajyt*ta2(:,j+1,:)
                  enddo
                  !      td2(:,ysize(2),:) = 2.*alsajyt*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:) ! npaire=1
                  td2(:,ysize(2),:) = ta2(:,ysize(2),:) ! npaire=0
                  ta2=td2

                elseif (ncly.eq.2) then

                  td2(:,1,:) = 0.
                  td2(:,2,:) = alsa2y*ta2(:,1,:) + ta2(:,2,:) + alsa2y*ta2(:,3,:)
                  td2(:,3,:) = alsa3y*ta2(:,2,:) + ta2(:,3,:) + alsa3y*ta2(:,4,:)
                  do j=4,ysize(2)-3
                    td2(:,j,:) = alsajyt*ta2(:,j-1,:) + ta2(:,j,:) + alsajyt*ta2(:,j+1,:)
                  enddo
                  td2(:,ysize(2)-2,:) = alsaty*ta2(:,ysize(2)-3,:) + ta2(:,ysize(2)-2,:) + alsaty*ta2(:,ysize(2)-1,:)
                  td2(:,ysize(2)-1,:) = alsamy*ta2(:,ysize(2)-2,:) + ta2(:,ysize(2)-1,:) + alsamy*ta2(:,ysize(2),:)
                  td2(:,ysize(2),:) = 0.
                  ta2=td2

                endif

              else
                do k=1,ysize(3)
                  do i=1,ysize(1)

                    td2(i,1,k)= 0.
                    td2(i,2,k)=alsa2y*ta2(i,1,k)+ta2(i,2,k)+alsa2y*ta2(i,3,k)
                    td2(i,3,k)=alsa3y*ta2(i,2,k)+ta2(i,3,k)+alsa3y*ta2(i,4,k)
                    do j=4,ysize(2)-3
                      td2(i,j,k)=alsajyt*ta2(i,j-1,k)+ta2(i,j,k)+alsajyt*ta2(i,j+1,k)
                    enddo
                    td2(i,ysize(2)-2,k)=alsaty*ta2(i,ysize(2)-3,k)+ta2(i,ysize(2)-2,k)+alsaty*ta2(i,ysize(2)-1,k)
                    td2(i,ysize(2)-1,k)=alsamy*ta2(i,ysize(2)-2,k)+ta2(i,ysize(2)-1,k)+alsamy*ta2(i,ysize(2),k)
                    td2(i,ysize(2),k)=zero!

                  enddo
                enddo
                do k=1,ysize(3)
                  do j=1,ysize(2)
                    do i=1,ysize(1)
                      ta2(i,j,k)=td2(i,j,k)
                    enddo
                  enddo
                enddo
              endif

              !(A+nu*dt.B).un
              if (.true.) then

                if (ncly.eq.0) then

                  call deryyt(td2,ux2,di2,sy,sfyt,ssyt,swyt,ysize(1),ysize(2),ysize(3),0)

                elseif (ncly.eq.1) then

                  !      call deryyt(td2,ux2,di2,sy,sfypt,ssypt,swypt,ysize(1),ysize(2),ysize(3),1) ! npaire=1
                  call deryyt(td2,ux2,di2,sy,sfyt,ssyt,swyt,ysize(1),ysize(2),ysize(3),0) ! npaire=0

                elseif (ncly.eq.2) then

                  call deryyt(td2,ux2,di2,sy,sfyt,ssyt,swyt,ysize(1),ysize(2),ysize(3),0)

                endif

                if (istret.ne.0) then
                  do j=1,ysize(2)
                    ux2(:,j,:)=ux2(:,j,:)/pp2y(j)
                  enddo
                endif

                if (ncly.eq.0) then

                  td2(:,1,:) = alsajyt*ux2(:,2,:) + ux2(:,1,:) + alsajyt*ux2(:,ysize(2),:) &
                  + xcst_pr*td2(:,1,:)
                  do j=2,ysize(2)-1
                    td2(:,j,:) = alsajyt*ux2(:,j-1,:) + ux2(:,j,:) + alsajyt*ux2(:,j+1,:) &
                    + xcst_pr*td2(:,j,:)
                  enddo
                  td2(:,ysize(2),:) = alsajyt*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) + alsajyt*ux2(:,1,:) &
                  + xcst_pr*td2(:,ysize(2),:)

                elseif (ncly.eq.1) then

                  !      td2(:,1,:) = 2.*alsajyt*ux2(:,2,:) + ux2(:,1,:) & ! npaire=1
                  td2(:,1,:) = ux2(:,1,:) & ! npaire=0
                  + xcst_pr*td2(:,1,:)
                  do j=2,ysize(2)-1
                    td2(:,j,:) = alsajyt*ux2(:,j-1,:) + ux2(:,j,:) + alsajyt*ux2(:,j+1,:) &
                    + xcst_pr*td2(:,j,:)
                  enddo
                  !      td2(:,ysize(2),:) = 2.*alsajyt*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) & ! npaire=1
                  td2(:,ysize(2),:) = ux2(:,ysize(2),:) & ! npaire=0
                  + xcst_pr*td2(:,ysize(2),:)

                elseif (ncly.eq.2) then

                  td2(:,1,:) = 0.
                  td2(:,2,:) = alsa2y*ux2(:,1,:) + ux2(:,2,:) + alsa2y*ux2(:,3,:) &
                  + xcst_pr*td2(:,2,:)
                  td2(:,3,:) = alsa3y*ux2(:,2,:) + ux2(:,3,:) + alsa3y*ux2(:,4,:) &
                  + xcst_pr*td2(:,3,:)
                  do j=4,ysize(2)-3
                    td2(:,j,:) = alsajyt*ux2(:,j-1,:) + ux2(:,j,:) + alsajyt*ux2(:,j+1,:) &
                    + xcst_pr*td2(:,j,:)
                  enddo
                  td2(:,ysize(2)-2,:) = alsaty*ux2(:,ysize(2)-3,:) + ux2(:,ysize(2)-2,:) + alsaty*ux2(:,ysize(2)-1,:) &
                  + xcst_pr*td2(:,ysize(2)-2,:)
                  td2(:,ysize(2)-1,:) = alsamy*ux2(:,ysize(2)-2,:) + ux2(:,ysize(2)-1,:) + alsamy*ux2(:,ysize(2),:) &
                  + xcst_pr*td2(:,ysize(2)-1,:)
                  td2(:,ysize(2),:) = 0.

                endif

              else

                if(istret==0) then

                  do k=1,ysize(3)
                    do i=1,ysize(1)

                      !for ux
                      td2(i,1,k)= 0.

                      td2(i,2,k)= (alsa2y+xcst_pr*as2y)*ux2(i,1,k)+(1.-xcst_pr*2.*as2y)*ux2(i,2,k)+(alsa2y+xcst_pr*as2y)*ux2(i,3,k)

                      td2(i,3,k)=xcst_pr*bs3y*ux2(i,1,k)+(alsa3y+xcst_pr*as3y)*ux2(i,2,k) &
                      + (1.-xcst_pr*2.*(as3y+bs3y))*ux2(i,3,k) + (alsa3y+xcst_pr*as3y)*ux2(i,4,k) &
                      +xcst_pr*bs3y*ux2(i,5,k)
                      do j=4,ysize(2)-3
                        td2(i,j,k)=xcst_pr*csjyt*ux2(i,j-3,k)+xcst_pr*bsjyt*ux2(i,j-2,k)+(alsajyt+xcst_pr*asjyt)*ux2(i,j-1,k) &
                        + (1.-xcst_pr*2.*(asjyt+bsjyt+csjyt))*ux2(i,j,k) + (alsajyt+xcst_pr*asjyt)*ux2(i,j+1,k) &
                        +xcst_pr*bsjyt*ux2(i,j+2,k)+xcst_pr*csjyt*ux2(i,j+3,k)
                      enddo
                      td2(i,ysize(2)-2,k)=xcst_pr*bsty*ux2(i,ysize(2)-4,k)+(alsaty+xcst_pr*asty)*ux2(i,ysize(2)-3,k) &
                      + (1.-xcst_pr*2.*(asty+bsty))*ux2(i,ysize(2)-2,k) + (alsaty+xcst_pr*asty)*ux2(i,ysize(2)-1,k) &
                      +xcst_pr*bsty*ux2(i,ysize(2),k)
                      td2(i,ysize(2)-1,k)=(alsamy+xcst_pr*asmy)*ux2(i,ysize(2)-2,k)+(1.-xcst_pr*2.*asmy)*ux2(i,ysize(2)-1,k) &
                      +(alsamy+xcst_pr*asmy)*ux2(i,ysize(2),k)
                      td2(i,ysize(2),k)=zero

                    enddo
                  enddo

                else

                  do k=1,ysize(3)
                    do i=1,ysize(1)

                      !for ux
                      td2(i,1,k)= 0.
                      td2(i,2,k)=(alsa2y+xcst_pr*as2y*pp2y(2))*ux2(i,1,k)+(1.-xcst_pr*2.*as2y*pp2y(2))*ux2(i,2,k)&
                      +(alsa2y+xcst_pr*as2y*pp2y(2))*ux2(i,3,k)
                      td2(i,3,k)=xcst_pr*bs3y*ux2(i,1,k)*pp2y(3)+(alsa3y+xcst_pr*as3y*pp2y(3))*ux2(i,2,k) &
                      + (1.-xcst_pr*2.*(as3y+bs3y)*pp2y(3))*ux2(i,3,k) + (alsa3y+xcst_pr*as3y*pp2y(3))*ux2(i,4,k) &
                      +xcst_pr*bs3y*ux2(i,5,k)*pp2y(3)
                      do j=4,ysize(2)-3
                        td2(i,j,k)=xcst_pr*csjyt*ux2(i,j-3,k)*pp2y(j)+xcst_pr*bsjyt*ux2(i,j-2,k)*pp2y(j) &
                        +(alsajyt+xcst_pr*asjyt*pp2y(j))*ux2(i,j-1,k) &
                        + (1.-xcst_pr*2.*(asjyt+bsjyt+csjyt)*pp2y(j))*ux2(i,j,k) + (alsajyt+xcst_pr*asjyt*pp2y(j))*ux2(i,j+1,k) &
                        +xcst_pr*bsjyt*ux2(i,j+2,k)*pp2y(j)+xcst_pr*csjyt*ux2(i,j+3,k)*pp2y(j)
                      enddo
                      td2(i,ysize(2)-2,k)=xcst_pr*bsty*ux2(i,ysize(2)-4,k)*pp2y(ysize(2)-2)+ &
                      (alsaty+xcst_pr*asty*pp2y(ysize(2)-2))*ux2(i,ysize(2)-3,k) &
                      + (1.-xcst_pr*2.*(asty+bsty)*pp2y(ysize(2)-2))*ux2(i,ysize(2)-2,k) &
                      + (alsaty+xcst_pr*asty*pp2y(ysize(2)-2))*ux2(i,ysize(2)-1,k) &
                      +xcst_pr*bsty*ux2(i,ysize(2),k)*pp2y(ysize(2)-2)
                      td2(i,ysize(2)-1,k)=(alsamy+xcst_pr*asmy*pp2y(ysize(2)-1))*ux2(i,ysize(2)-2,k)+ &
                      (1.-xcst_pr*2.*asmy*pp2y(ysize(2)-1))*ux2(i,ysize(2)-1,k) &
                      +(alsamy+xcst_pr*asmy*pp2y(ysize(2)-1))*ux2(i,ysize(2),k)
                      td2(i,ysize(2),k)=zero

                    enddo
                  enddo

                endif

              endif


            end subroutine multmatrix7T

            !!********************************************************************
            !!
            !subroutine derxxt(tx,ux,rx,sx,sfx,ssx,swx,nx,ny,nz,npaire)
            !!
            !!********************************************************************

            !USE param
            !USE derivX

            !implicit none

            !integer :: nx,ny,nz,npaire,i,j,k
            !real(mytype), dimension(nx,ny,nz) :: tx,ux,rx
            !real(mytype), dimension(ny,nz) :: sx
            !real(mytype),  dimension(nx):: sfx,ssx,swx

            !if (nclx==0) then
            !   do k=1,nz
            !   do j=1,ny
            !      tx(1,j,k)=asixt*(ux(2,j,k)-ux(1   ,j,k)&
            !           -ux(1,j,k)+ux(nx  ,j,k))&
            !           +bsixt*(ux(3,j,k)-ux(1   ,j,k)&
            !           -ux(1,j,k)+ux(nx-1,j,k))&
            !           +csixt*(ux(4,j,k)-ux(1   ,j,k)&
            !           -ux(1,j,k)+ux(nx-2,j,k))
            !      rx(1,j,k)=-1.
            !      tx(2,j,k)=asixt*(ux(3,j,k)-ux(2   ,j,k)&
            !           -ux(2,j,k)+ux(1   ,j,k))&
            !           +bsixt*(ux(4,j,k)-ux(2   ,j,k)&
            !           -ux(2,j,k)+ux(nx  ,j,k))&
            !           +csixt*(ux(5,j,k)-ux(2   ,j,k)&
            !           -ux(2,j,k)+ux(nx-1,j,k))
            !      rx(2,j,k)=zero
            !      tx(3,j,k)=asixt*(ux(4,j,k)-ux(3 ,j,k)&
            !           -ux(3,j,k)+ux(2 ,j,k))&
            !           +bsixt*(ux(5,j,k)-ux(3 ,j,k)&
            !           -ux(3,j,k)+ux(1 ,j,k))&
            !           +csixt*(ux(6,j,k)-ux(3 ,j,k)&
            !           -ux(3,j,k)+ux(nx,j,k))
            !      rx(3,j,k)=zero
            !      do i=4,nx-3
            !         tx(i,j,k)=asixt*(ux(i+1,j,k)-ux(i  ,j,k)&
            !              -ux(i  ,j,k)+ux(i-1,j,k))&
            !              +bsixt*(ux(i+2,j,k)-ux(i  ,j,k)&
            !              -ux(i  ,j,k)+ux(i-2,j,k))&
            !              +csixt*(ux(i+3,j,k)-ux(i  ,j,k)&
            !              -ux(i  ,j,k)+ux(i-3,j,k))
            !         rx(i,j,k)=zero
            !      enddo
            !      tx(nx-2,j,k)=asixt*(ux(nx-1,j,k)-ux(nx-2,j,k)&
            !           -ux(nx-2,j,k)+ux(nx-3,j,k))&
            !           +bsixt*(ux(nx  ,j,k)-ux(nx-2,j,k)&
            !           -ux(nx-2,j,k)+ux(nx-4,j,k))&
            !           +csixt*(ux(1   ,j,k)-ux(nx-2,j,k)&
            !           -ux(nx-2,j,k)+ux(nx-5,j,k))
            !      rx(nx-2,j,k)=zero
            !      tx(nx-1,j,k)=asixt*(ux(nx  ,j,k)-ux(nx-1,j,k)&
            !           -ux(nx-1,j,k)+ux(nx-2,j,k))&
            !           +bsixt*(ux(1   ,j,k)-ux(nx-1,j,k)&
            !           -ux(nx-1,j,k)+ux(nx-3,j,k))&
            !           +csixt*(ux(2   ,j,k)-ux(nx-1,j,k)&
            !           -ux(nx-1,j,k)+ux(nx-4,j,k))
            !      rx(nx-1,j,k)=zero
            !      tx(nx  ,j,k)=asixt*(ux(1 ,j,k)-ux(nx  ,j,k)&
            !           -ux(nx,j,k)+ux(nx-1,j,k))&
            !           +bsixt*(ux(2 ,j,k)-ux(nx  ,j,k)&
            !           -ux(nx,j,k)+ux(nx-2,j,k))&
            !           +csixt*(ux(3 ,j,k)-ux(nx  ,j,k)&
            !           -ux(nx,j,k)+ux(nx-3,j,k))
            !      rx(nx  ,j,k)=alsaixt
            !      do i=2,nx
            !         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
            !         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*ssx(i)
            !      enddo
            !         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
            !         rx(nx,j,k)=rx(nx,j,k)*swx(nx)
            !      do i=nx-1,1,-1
            !         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
            !         rx(i,j,k)=(rx(i,j,k)-sfx(i)*rx(i+1,j,k))*swx(i)
            !      enddo
            !      sx(j,k)=(   tx(1,j,k)-alsaixt*tx(nx,j,k))/&
            !           (1.+rx(1,j,k)-alsaixt*rx(nx,j,k))
            !      do i=1,nx
            !         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
            !      enddo
            !   enddo
            !   enddo
            !endif

            !if (nclx==1) then
            !   if (npaire==1) then
            !      do k=1,nz
            !      do j=1,ny
            !         tx(1,j,k)=asixt*(ux(2,j,k)-ux(1,j,k)&
            !              -ux(1,j,k)+ux(2,j,k))&
            !              +bsixt*(ux(3,j,k)-ux(1,j,k)&
            !              -ux(1,j,k)+ux(3,j,k))&
            !              +csixt*(ux(4,j,k)-ux(1,j,k)&
            !              -ux(1,j,k)+ux(4,j,k))
            !         tx(2,j,k)=asixt*(ux(3,j,k)-ux(2,j,k)&
            !              -ux(2,j,k)+ux(1,j,k))&
            !              +bsixt*(ux(4,j,k)-ux(2,j,k)&
            !              -ux(2,j,k)+ux(2,j,k))&
            !              +csixt*(ux(5,j,k)-ux(2,j,k)&
            !              -ux(2,j,k)+ux(3,j,k))
            !         tx(3,j,k)=asixt*(ux(4,j,k)-ux(3,j,k)&
            !              -ux(3,j,k)+ux(2,j,k))&
            !              +bsixt*(ux(5,j,k)-ux(3,j,k)&
            !              -ux(3,j,k)+ux(1,j,k))&
            !              +csixt*(ux(6,j,k)-ux(3,j,k)&
            !              -ux(3,j,k)+ux(2,j,k))
            !         do i=4,nx-3
            !            tx(i,j,k)=asixt*(ux(i+1,j,k)-ux(i  ,j,k)&
            !                 -ux(i  ,j,k)+ux(i-1,j,k))&
            !                 +bsixt*(ux(i+2,j,k)-ux(i  ,j,k)&
            !                 -ux(i  ,j,k)+ux(i-2,j,k))&
            !                 +csixt*(ux(i+3,j,k)-ux(i  ,j,k)&
            !                 -ux(i  ,j,k)+ux(i-3,j,k))
            !         enddo
            !         tx(nx-2,j,k)=asixt*(ux(nx-1,j,k)-ux(nx-2,j,k)&
            !              -ux(nx-2,j,k)+ux(nx-3,j,k))&
            !              +bsixt*(ux(nx  ,j,k)-ux(nx-2,j,k)&
            !              -ux(nx-2,j,k)+ux(nx-4,j,k))&
            !              +csixt*(ux(nx-1,j,k)-ux(nx-2,j,k)&
            !              -ux(nx-2,j,k)+ux(nx-5,j,k))
            !         tx(nx-1,j,k)=asixt*(ux(nx  ,j,k)-ux(nx-1,j,k)&
            !              -ux(nx-1,j,k)+ux(nx-2,j,k))&
            !              +bsixt*(ux(nx-1,j,k)-ux(nx-1,j,k)&
            !              -ux(nx-1,j,k)+ux(nx-3,j,k))&
            !              +csixt*(ux(nx-2,j,k)-ux(nx-1,j,k)&
            !              -ux(nx-1,j,k)+ux(nx-4,j,k))
            !         tx(nx  ,j,k)=asixt*(ux(nx-1,j,k)-ux(nx  ,j,k)&
            !              -ux(nx  ,j,k)+ux(nx-1,j,k))&
            !              +bsixt*(ux(nx-2,j,k)-ux(nx  ,j,k)&
            !              -ux(nx  ,j,k)+ux(nx-2,j,k))&
            !              +csixt*(ux(nx-3,j,k)-ux(nx  ,j,k)&
            !              -ux(nx  ,j,k)+ux(nx-3,j,k))
            !         do i=2,nx
            !            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
            !         enddo
            !         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
            !         do i=nx-1,1,-1
            !            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
            !         enddo
            !      enddo
            !      enddo
            !   endif
            !   if (npaire==0) then
            !      do k=1,nz
            !      do j=1,ny
            !         tx(1,j,k)=zero
            !         tx(2,j,k)=asixt*(ux(3,j,k)-ux(2,j,k)&
            !              -ux(2,j,k)+ux(1,j,k))&
            !              +bsixt*(ux(4,j,k)-ux(2,j,k)&
            !              -ux(2,j,k)-ux(2,j,k))&
            !              +csixt*(ux(5,j,k)-ux(2,j,k)&
            !              -ux(2,j,k)-ux(3,j,k))
            !         tx(3,j,k)=asixt*(ux(4,j,k)-ux(3,j,k)&
            !              -ux(3,j,k)+ux(2,j,k))&
            !              +bsixt*(ux(5,j,k)-ux(3,j,k)&
            !              -ux(3,j,k)+ux(1,j,k))&
            !              +csixt*(ux(6,j,k)-ux(3,j,k)&
            !              -ux(3,j,k)-ux(2,j,k))
            !         do i=4,nx-3
            !            tx(i,j,k)=asixt*(ux(i+1,j,k)-ux(i  ,j,k)&
            !                 -ux(i  ,j,k)+ux(i-1,j,k))&
            !                 +bsixt*(ux(i+2,j,k)-ux(i  ,j,k)&
            !                 -ux(i  ,j,k)+ux(i-2,j,k))&
            !                 +csixt*(ux(i+3,j,k)-ux(i  ,j,k)&
            !                 -ux(i  ,j,k)+ux(i-3,j,k))
            !         enddo
            !         tx(nx-2,j,k)=asixt*( ux(nx-1,j,k)-ux(nx-2,j,k)&
            !              -ux(nx-2,j,k)+ux(nx-3,j,k))&
            !              +bsixt*( ux(nx  ,j,k)-ux(nx-2,j,k)&
            !              -ux(nx-2,j,k)+ux(nx-4,j,k))&
            !              +csixt*(-ux(nx-1,j,k)-ux(nx-2,j,k)&
            !              -ux(nx-2,j,k)+ux(nx-5,j,k))
            !         tx(nx-1,j,k)=asixt*( ux(nx  ,j,k)-ux(nx-1,j,k)&
            !              -ux(nx-1,j,k)+ux(nx-2,j,k))&
            !              +bsixt*(-ux(nx-1,j,k)-ux(nx-1,j,k)&
            !              -ux(nx-1,j,k)+ux(nx-3,j,k))&
            !              +csixt*(-ux(nx-2,j,k)-ux(nx-1,j,k)&
            !              -ux(nx-1,j,k)+ux(nx-4,j,k))
            !         tx(nx  ,j,k)=zero
            !         do i=2,nx
            !            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
            !         enddo
            !         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
            !         do i=nx-1,1,-1
            !            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
            !         enddo
            !      enddo
            !      enddo
            !   endif
            !endif
            !
            !if (nclx==2) then
            !   do k=1,nz
            !   do j=1,ny
            !      tx(1,j,k)=as1x*ux(1,j,k)+bs1x*ux(2,j,k)&
            !           +cs1x*ux(3,j,k)+ds1x*ux(4,j,k)
            !      tx(2,j,k)=as2x*(ux(3,j,k)-ux(2,j,k)&
            !           -ux(2,j,k)+ux(1,j,k))
            !      tx(3,j,k)=as3x*(ux(4,j,k)-ux(3,j,k)&
            !           -ux(3,j,k)+ux(2,j,k))&
            !           +bs3x*(ux(5,j,k)-ux(3,j,k)&
            !           -ux(3,j,k)+ux(1,j,k))
            !      do i=4,nx-3
            !         tx(i,j,k)=asixt*(ux(i+1,j,k)-ux(i  ,j,k)&
            !              -ux(i  ,j,k)+ux(i-1,j,k))&
            !              +bsixt*(ux(i+2,j,k)-ux(i  ,j,k)&
            !              -ux(i  ,j,k)+ux(i-2,j,k))&
            !              +csixt*(ux(i+3,j,k)-ux(i  ,j,k)&
            !              -ux(i  ,j,k)+ux(i-3,j,k))
            !      enddo
            !      tx(nx-2,j,k)=astx*(ux(nx-1,j,k)-ux(nx-2,j,k)&
            !           -ux(nx-2,j,k)+ux(nx-3,j,k))&
            !           +bstx*(ux(nx  ,j,k)-ux(nx-2,j,k)&
            !           -ux(nx-2,j,k)+ux(nx-4,j,k))
            !      tx(nx-1,j,k)=asmx*(ux(nx  ,j,k)-ux(nx-1,j,k)&
            !           -ux(nx-1,j,k)+ux(nx-2,j,k))
            !      tx(nx  ,j,k)=asnx*ux(nx  ,j,k)+bsnx*ux(nx-1,j,k)&
            !           +csnx*ux(nx-2,j,k)+dsnx*ux(nx-3,j,k)
            !      do i=2,nx
            !         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
            !      enddo
            !      tx(nx,j,k)=tx(nx,j,k)*swx(nx)
            !      do i=nx-1,1,-1
            !         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
            !      enddo
            !   enddo
            !   enddo
            !endif

            !return
            !end subroutine derxxt

            !!********************************************************************
            !!
            !subroutine deryyt(ty,uy,ry,sy,sfy,ssy,swy,nx,ny,nz,npaire)
            !!
            !!********************************************************************

            !USE param
            !USE derivY

            !implicit none

            !integer :: nx,ny,nz,npaire,i,j,k
            !real(mytype), dimension(nx,ny,nz) :: ty,uy,ry
            !real(mytype), dimension(nx,nz) :: sy
            !real(mytype), dimension(ny) :: sfy,ssy,swy

            !if (ncly==0) then
            !   do k=1,nz
            !   do i=1,nx
            !      ty(i,1,k)=asjyt*(uy(i,2,k)-uy(i,1,k)&
            !           -uy(i,1,k)+uy(i,ny,k))&
            !           +bsjyt*(uy(i,3,k)-uy(i,1,k)&
            !           -uy(i,1,k)+uy(i,ny-1,k))&
            !           +csjyt*(uy(i,4,k)-uy(i,1,k)&
            !           -uy(i,1,k)+uy(i,ny-2,k))
            !      ry(i,1,k)=-1.
            !      ty(i,2,k)=asjyt*(uy(i,3,k)-uy(i,2,k)&
            !           -uy(i,2,k)+uy(i,1,k))&
            !           +bsjyt*(uy(i,4,k)-uy(i,2,k)&
            !           -uy(i,2,k)+uy(i,ny,k))&
            !           +csjyt*(uy(i,5,k)-uy(i,2,k)&
            !           -uy(i,2,k)+uy(i,ny-1,k))
            !      ry(i,2,k)=zero
            !      ty(i,3,k)=asjyt*(uy(i,4,k)-uy(i,3,k)&
            !           -uy(i,3,k)+uy(i,2,k))&
            !           +bsjyt*(uy(i,5,k)-uy(i,3,k)&
            !           -uy(i,3,k)+uy(i,1,k))&
            !           +csjyt*(uy(i,6,k)-uy(i,3,k)&
            !           -uy(i,3,k)+uy(i,ny,k))
            !      ry(i,3,k)=zero
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do j=4,ny-3
            !   do i=1,nx
            !      ty(i,j,k)=asjyt*(uy(i,j+1,k)-uy(i,j,k)&
            !           -uy(i,j,k)+uy(i,j-1,k))&
            !           +bsjyt*(uy(i,j+2,k)-uy(i,j,k)&
            !           -uy(i,j,k)+uy(i,j-2,k))&
            !           +csjyt*(uy(i,j+3,k)-uy(i,j,k)&
            !           -uy(i,j,k)+uy(i,j-3,k))
            !      ry(i,j,k)=zero
            !   enddo
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do i=1,nx
            !      ty(i,ny-2,k)=asjyt*(uy(i,ny-1,k)-uy(i,ny-2,k)&
            !           -uy(i,ny-2,k)+uy(i,ny-3,k))&
            !           +bsjyt*(uy(i,ny  ,k)-uy(i,ny-2,k)&
            !           -uy(i,ny-2,k)+uy(i,ny-4,k))&
            !           +csjyt*(uy(i,1   ,k)-uy(i,ny-2,k)&
            !           -uy(i,ny-2,k)+uy(i,ny-5,k))
            !      ry(i,ny-2,k)=zero
            !      ty(i,ny-1,k)=asjyt*(uy(i,ny  ,k)-uy(i,ny-1,k)&
            !           -uy(i,ny-1,k)+uy(i,ny-2,k))&
            !           +bsjyt*(uy(i,1   ,k)-uy(i,ny-1,k)&
            !           -uy(i,ny-1,k)+uy(i,ny-3,k))&
            !           +csjyt*(uy(i,2   ,k)-uy(i,ny-1,k)&
            !           -uy(i,ny-1,k)+uy(i,ny-4,k))
            !      ry(i,ny-1,k)=zero
            !      ty(i,ny  ,k)=asjyt*(uy(i,1 ,k)-uy(i,ny  ,k)&
            !           -uy(i,ny,k)+uy(i,ny-1,k))&
            !           +bsjyt*(uy(i,2 ,k)-uy(i,ny  ,k)&
            !           -uy(i,ny,k)+uy(i,ny-2,k))&
            !           +csjyt*(uy(i,3 ,k)-uy(i,ny  ,k)&
            !           -uy(i,ny,k)+uy(i,ny-3,k))
            !      ry(i,ny  ,k)=alsajyt
            !   enddo
            !   enddo
            !if (iimplicit.eq.1) return
            !   do k=1,nz
            !   do j=2,ny
            !   do i=1,nx
            !      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
            !      ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*ssy(j)
            !   enddo
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do i=1,nx
            !      ty(i,ny,k)=ty(i,ny,k)*swy(ny)
            !      ry(i,ny,k)=ry(i,ny,k)*swy(ny)
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do j=ny-1,1,-1
            !   do i=1,nx
            !      ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
            !      ry(i,j,k)=(ry(i,j,k)-sfy(j)*ry(i,j+1,k))*swy(j)
            !   enddo
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do i=1,nx
            !      sy(i,k)=(   ty(i,1,k)-alsajyt*ty(i,ny,k))/&
            !           (1.+ry(i,1,k)-alsajyt*ry(i,ny,k))
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do j=1,ny
            !   do i=1,nx
            !      ty(i,j,k)=ty(i,j,k)-sy(i,k)*ry(i,j,k)
            !   enddo
            !   enddo
            !   enddo
            !endif

            !if (ncly==1) then
            !   if (npaire==1) then
            !      do k=1,nz
            !      do i=1,nx
            !         ty(i,1,k)=asjyt*(uy(i,2,k)-uy(i,1,k)&
            !              -uy(i,1,k)+uy(i,2,k))&
            !              +bsjyt*(uy(i,3,k)-uy(i,1,k)&
            !              -uy(i,1,k)+uy(i,3,k))&
            !              +csjyt*(uy(i,4,k)-uy(i,1,k)&
            !              -uy(i,1,k)+uy(i,4,k))
            !         ty(i,2,k)=asjyt*(uy(i,3,k)-uy(i,2,k)&
            !              -uy(i,2,k)+uy(i,1,k))&
            !              +bsjyt*(uy(i,4,k)-uy(i,2,k)&
            !              -uy(i,2,k)+uy(i,2,k))&
            !              +csjyt*(uy(i,5,k)-uy(i,2,k)&
            !              -uy(i,2,k)+uy(i,3,k))
            !         ty(i,3,k)=asjyt*(uy(i,4,k)-uy(i,3,k)&
            !              -uy(i,3,k)+uy(i,2,k))&
            !              +bsjyt*(uy(i,5,k)-uy(i,3,k)&
            !              -uy(i,3,k)+uy(i,1,k))&
            !              +csjyt*(uy(i,6,k)-uy(i,3,k)&
            !              -uy(i,3,k)+uy(i,2,k))
            !      enddo
            !      enddo
            !      do k=1,nz
            !      do j=4,ny-3
            !      do i=1,nx
            !         ty(i,j,k)=asjyt*(uy(i,j+1,k)-uy(i,j  ,k)&
            !              -uy(i,j  ,k)+uy(i,j-1,k))&
            !              +bsjyt*(uy(i,j+2,k)-uy(i,j  ,k)&
            !              -uy(i,j  ,k)+uy(i,j-2,k))&
            !              +csjyt*(uy(i,j+3,k)-uy(i,j  ,k)&
            !              -uy(i,j  ,k)+uy(i,j-3,k))
            !      enddo
            !      enddo
            !      enddo
            !      do k=1,nz
            !      do i=1,nx
            !         ty(i,ny-2,k)=asjyt*(uy(i,ny-1,k)-uy(i,ny-2,k)&
            !              -uy(i,ny-2,k)+uy(i,ny-3,k))&
            !              +bsjyt*(uy(i,ny  ,k)-uy(i,ny-2,k)&
            !              -uy(i,ny-2,k)+uy(i,ny-4,k))&
            !              +csjyt*(uy(i,ny-1,k)-uy(i,ny-2,k)&
            !              -uy(i,ny-2,k)+uy(i,ny-5,k))
            !         ty(i,ny-1,k)=asjyt*(uy(i,ny  ,k)-uy(i,ny-1,k)&
            !              -uy(i,ny-1,k)+uy(i,ny-2,k))&
            !              +bsjyt*(uy(i,ny-1,k)-uy(i,ny-1,k)&
            !              -uy(i,ny-1,k)+uy(i,ny-3,k))&
            !              +csjyt*(uy(i,ny-2,k)-uy(i,ny-1,k)&
            !              -uy(i,ny-1,k)+uy(i,ny-4,k))
            !         ty(i,ny  ,k)=asjyt*(uy(i,ny-1,k)-uy(i,ny  ,k)&
            !              -uy(i,ny  ,k)+uy(i,ny-1,k))&
            !              +bsjyt*(uy(i,ny-2,k)-uy(i,ny  ,k)&
            !              -uy(i,ny  ,k)+uy(i,ny-2,k))&
            !              +csjyt*(uy(i,ny-3,k)-uy(i,ny  ,k)&
            !              -uy(i,ny  ,k)+uy(i,ny-3,k))
            !      enddo
            !      enddo
            !if (iimplicit.eq.1) return
            !      do k=1,nz
            !      do j=2,ny
            !      do i=1,nx
            !         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
            !      enddo
            !      enddo
            !      enddo
            !      do k=1,nz
            !      do i=1,nx
            !         ty(i,ny,k)=ty(i,ny,k)*swy(ny)
            !      enddo
            !      enddo
            !      do k=1,nz
            !      do j=ny-1,1,-1
            !      do i=1,nx
            !         ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
            !      enddo
            !      enddo
            !      enddo
            !   endif
            !   if (npaire==0) then
            !      do k=1,nz
            !      do i=1,nx
            !         ty(i,1,k)=zero
            !         ty(i,2,k)=asjyt*(uy(i,3,k)-uy(i,2,k)&
            !              -uy(i,2,k)+uy(i,1,k))&
            !              +bsjyt*(uy(i,4,k)-uy(i,2,k)&
            !              -uy(i,2,k)-uy(i,2,k))&
            !              +csjyt*(uy(i,5,k)-uy(i,2,k)&
            !              -uy(i,2,k)-uy(i,3,k))
            !         ty(i,3,k)=asjyt*(uy(i,4,k)-uy(i,3,k)&
            !              -uy(i,3,k)+uy(i,2,k))&
            !              +bsjyt*(uy(i,5,k)-uy(i,3,k)&
            !              -uy(i,3,k)+uy(i,1,k))&
            !              +csjyt*(uy(i,6,k)-uy(i,3,k)&
            !              -uy(i,3,k)-uy(i,2,k))
            !      enddo
            !      enddo
            !      do k=1,nz
            !      do j=4,ny-3
            !      do i=1,nx
            !         ty(i,j,k)=asjyt*(uy(i,j+1,k)-uy(i,j  ,k)&
            !              -uy(i,j  ,k)+uy(i,j-1,k))&
            !              +bsjyt*(uy(i,j+2,k)-uy(i,j  ,k)&
            !              -uy(i,j  ,k)+uy(i,j-2,k))&
            !              +csjyt*(uy(i,j+3,k)-uy(i,j  ,k)&
            !              -uy(i,j  ,k)+uy(i,j-3,k))
            !      enddo
            !      enddo
            !      enddo
            !      do k=1,nz
            !      do i=1,nx
            !         ty(i,ny-2,k)=asjyt*( uy(i,ny-1,k)-uy(i,ny-2,k)&
            !              -uy(i,ny-2,k)+uy(i,ny-3,k))&
            !              +bsjyt*( uy(i,ny  ,k)-uy(i,ny-2,k)&
            !              -uy(i,ny-2,k)+uy(i,ny-4,k))&
            !              +csjyt*(-uy(i,ny-1,k)-uy(i,ny-2,k)&
            !              -uy(i,ny-2,k)+uy(i,ny-5,k))
            !         ty(i,ny-1,k)=asjyt*( uy(i,ny  ,k)-uy(i,ny-1,k)&
            !              -uy(i,ny-1,k)+uy(i,ny-2,k))&
            !              +bsjyt*(-uy(i,ny-1,k)-uy(i,ny-1,k)&
            !              -uy(i,ny-1,k)+uy(i,ny-3,k))&
            !              +csjyt*(-uy(i,ny-2,k)-uy(i,ny-1,k)&
            !              -uy(i,ny-1,k)+uy(i,ny-4,k))
            !         ty(i,ny  ,k)=zero
            !      enddo
            !      enddo
            !if (iimplicit.eq.1) return
            !      do k=1,nz
            !      do j=2,ny
            !      do i=1,nx
            !         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
            !      enddo
            !      enddo
            !      enddo
            !      do k=1,nz
            !      do i=1,nx
            !         ty(i,ny,k)=ty(i,ny,k)*swy(ny)
            !      enddo
            !      enddo
            !      do k=1,nz
            !      do j=ny-1,1,-1
            !      do i=1,nx
            !         ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
            !      enddo
            !      enddo
            !      enddo
            !   endif
            !endif

            !if (ncly==2) then
            !   do k=1,nz
            !   do i=1,nx
            !      ty(i,1,k)=as1y*uy(i,1,k)+bs1y*uy(i,2,k)&
            !           +cs1y*uy(i,3,k)+ds1y*uy(i,4,k)
            !      ty(i,2,k)=as2y*(uy(i,3,k)-uy(i,2,k)&
            !           -uy(i,2,k)+uy(i,1,k))
            !      ty(i,3,k)=as3y*(uy(i,4,k)-uy(i,3,k)&
            !           -uy(i,3,k)+uy(i,2,k))&
            !           +bs3y*(uy(i,5,k)-uy(i,3,k)&
            !           -uy(i,3,k)+uy(i,1,k))
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do j=4,ny-3
            !   do i=1,nx
            !      ty(i,j,k)=asjyt*(uy(i,j+1,k)-uy(i,j  ,k)&
            !           -uy(i,j  ,k)+uy(i,j-1,k))&
            !           +bsjyt*(uy(i,j+2,k)-uy(i,j  ,k)&
            !           -uy(i,j  ,k)+uy(i,j-2,k))&
            !           +csjyt*(uy(i,j+3,k)-uy(i,j  ,k)&
            !           -uy(i,j  ,k)+uy(i,j-3,k))
            !   enddo
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do i=1,nx
            !      ty(i,ny-2,k)=asty*(uy(i,ny-1,k)-uy(i,ny-2,k)&
            !           -uy(i,ny-2,k)+uy(i,ny-3,k))&
            !           +bsty*(uy(i,ny  ,k)-uy(i,ny-2,k)&
            !           -uy(i,ny-2,k)+uy(i,ny-4,k))
            !      ty(i,ny-1,k)=asmy*(uy(i,ny  ,k)-uy(i,ny-1,k)&
            !           -uy(i,ny-1,k)+uy(i,ny-2,k))
            !      ty(i,ny  ,k)=asny*uy(i,ny  ,k)+bsny*uy(i,ny-1,k)&
            !           +csny*uy(i,ny-2,k)+dsny*uy(i,ny-3,k)
            !   enddo
            !   enddo
            !if (iimplicit.eq.1) return
            !   do k=1,nz
            !   do j=2,ny
            !   do i=1,nx
            !      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
            !   enddo
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do i=1,nx
            !      ty(i,ny,k)=ty(i,ny,k)*swy(ny)
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do j=ny-1,1,-1
            !   do i=1,nx
            !      ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
            !   enddo
            !   enddo
            !   enddo
            !endif

            !return
            !end subroutine deryyt

            !!********************************************************************
            !!
            !subroutine derzzt(tz,uz,rz,sz,sfz,ssz,swz,nx,ny,nz,npaire)
            !!
            !!********************************************************************

            !USE param
            !USE derivZ

            !implicit none

            !integer :: nx,ny,nz,npaire,i,j,k
            !real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
            !real(mytype), dimension(nx,ny) :: sz
            !real(mytype), dimension(nz) :: sfz,ssz,swz

            !if (nclz==0) then
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,1)=askzt*(uz(i,j,2)-uz(i,j,1   )&
            !           -uz(i,j,1)+uz(i,j,nz  ))&
            !           +bskzt*(uz(i,j,3)-uz(i,j,1   )&
            !           -uz(i,j,1)+uz(i,j,nz-1))&
            !           +cskzt*(uz(i,j,4)-uz(i,j,1   )&
            !           -uz(i,j,1)+uz(i,j,nz-2))
            !      rz(i,j,1)=-1.
            !      tz(i,j,2)=askzt*(uz(i,j,3)-uz(i,j,2 )&
            !           -uz(i,j,2)+uz(i,j,1 ))&
            !           +bskzt*(uz(i,j,4)-uz(i,j,2 )&
            !           -uz(i,j,2)+uz(i,j,nz))&
            !           +cskzt*(uz(i,j,5)-uz(i,j,2 )&
            !           -uz(i,j,2)+uz(i,j,nz-1))
            !      rz(i,j,2)=zero
            !      tz(i,j,3)=askzt*(uz(i,j,4)-uz(i,j,3 )&
            !           -uz(i,j,3)+uz(i,j,2 ))&
            !           +bskzt*(uz(i,j,5)-uz(i,j,3 )&
            !           -uz(i,j,3)+uz(i,j,1 ))&
            !           +cskzt*(uz(i,j,6)-uz(i,j,3 )&
            !           -uz(i,j,3)+uz(i,j,nz))
            !      rz(i,j,3)=zero
            !   enddo
            !   enddo
            !   do k=4,nz-3
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,k)=askzt*(uz(i,j,k+1)-uz(i,j,k  )&
            !           -uz(i,j,k  )+uz(i,j,k-1))&
            !           +bskzt*(uz(i,j,k+2)-uz(i,j,k  )&
            !           -uz(i,j,k  )+uz(i,j,k-2))&
            !           +cskzt*(uz(i,j,k+3)-uz(i,j,k  )&
            !           -uz(i,j,k  )+uz(i,j,k-3))
            !      rz(i,j,k)=zero
            !   enddo
            !   enddo
            !   enddo
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,nz-2)=askzt*(uz(i,j,nz-1)-uz(i,j,nz-2)&
            !           -uz(i,j,nz-2)+uz(i,j,nz-3))&
            !           +bskzt*(uz(i,j,nz  )-uz(i,j,nz-2)&
            !           -uz(i,j,nz-2)+uz(i,j,nz-4))&
            !           +cskzt*(uz(i,j,1   )-uz(i,j,nz-2)&
            !           -uz(i,j,nz-2)+uz(i,j,nz-5))
            !      rz(i,j,nz-2)=zero
            !      tz(i,j,nz-1)=askzt*(uz(i,j,nz  )-uz(i,j,nz-1)&
            !           -uz(i,j,nz-1)+uz(i,j,nz-2))&
            !           +bskzt*(uz(i,j,1   )-uz(i,j,nz-1)&
            !           -uz(i,j,nz-1)+uz(i,j,nz-3))&
            !           +cskzt*(uz(i,j,2   )-uz(i,j,nz-1)&
            !           -uz(i,j,nz-1)+uz(i,j,nz-4))
            !      rz(i,j,nz-1)=zero
            !      tz(i,j,nz  )=askzt*(uz(i,j,1 )-uz(i,j,nz  )&
            !           -uz(i,j,nz)+uz(i,j,nz-1))&
            !           +bskzt*(uz(i,j,2 )-uz(i,j,nz  )&
            !           -uz(i,j,nz)+uz(i,j,nz-2))&
            !           +cskzt*(uz(i,j,3 )-uz(i,j,nz  )&
            !           -uz(i,j,nz)+uz(i,j,nz-3))
            !      rz(i,j,nz  )=alsakz
            !   enddo
            !   enddo
            !   do k=2,nz
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
            !      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*ssz(k)
            !   enddo
            !   enddo
            !   enddo
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,nz)=tz(i,j,nz)*swz(nz)
            !      rz(i,j,nz)=rz(i,j,nz)*swz(nz)
            !   enddo
            !   enddo
            !   do k=nz-1,1,-1
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
            !      rz(i,j,k)=(rz(i,j,k)-sfz(k)*rz(i,j,k+1))*swz(k)
            !   enddo
            !   enddo
            !   enddo
            !   do j=1,ny
            !   do i=1,nx
            !      sz(i,j)=(   tz(i,j,1)-alsakz*tz(i,j,nz))/&
            !           (1.+rz(i,j,1)-alsakz*rz(i,j,nz))
            !   enddo
            !   enddo
            !   do k=1,nz
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
            !   enddo
            !   enddo
            !   enddo
            !endif

            !if (nclz==1) then
            !   if (npaire==1) then
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,1)=askzt*(uz(i,j,2)-uz(i,j,1)&
            !              -uz(i,j,1)+uz(i,j,2))&
            !              +bskzt*(uz(i,j,3)-uz(i,j,1)&
            !              -uz(i,j,1)+uz(i,j,3))&
            !              +cskzt*(uz(i,j,4)-uz(i,j,1)&
            !              -uz(i,j,1)+uz(i,j,4))
            !         tz(i,j,2)=askzt*(uz(i,j,3)-uz(i,j,2)&
            !              -uz(i,j,2)+uz(i,j,1))&
            !              +bskzt*(uz(i,j,4)-uz(i,j,2)&
            !              -uz(i,j,2)+uz(i,j,2))&
            !              +cskzt*(uz(i,j,5)-uz(i,j,2)&
            !              -uz(i,j,2)+uz(i,j,3))
            !         tz(i,j,3)=askzt*(uz(i,j,4)-uz(i,j,3)&
            !              -uz(i,j,3)+uz(i,j,2))&
            !              +bskzt*(uz(i,j,5)-uz(i,j,3)&
            !              -uz(i,j,3)+uz(i,j,1))&
            !              +cskzt*(uz(i,j,6)-uz(i,j,3)&
            !              -uz(i,j,3)+uz(i,j,2))
            !      enddo
            !      enddo
            !      do k=4,nz-3
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,k)=askzt*(uz(i,j,k+1)-uz(i,j,k  )&
            !              -uz(i,j,k  )+uz(i,j,k-1))&
            !              +bskzt*(uz(i,j,k+2)-uz(i,j,k  )&
            !              -uz(i,j,k  )+uz(i,j,k-2))&
            !              +cskzt*(uz(i,j,k+3)-uz(i,j,k  )&
            !              -uz(i,j,k  )+uz(i,j,k-3))
            !      enddo
            !      enddo
            !      enddo
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,nz-2)=askzt*(uz(i,j,nz-1)-uz(i,j,nz-2)&
            !              -uz(i,j,nz-2)+uz(i,j,nz-3))&
            !              +bskzt*(uz(i,j,nz  )-uz(i,j,nz-2)&
            !              -uz(i,j,nz-2)+uz(i,j,nz-4))&
            !              +cskzt*(uz(i,j,nz-1)-uz(i,j,nz-2)&
            !              -uz(i,j,nz-2)+uz(i,j,nz-5))
            !         tz(i,j,nz-1)=askzt*(uz(i,j,nz  )-uz(i,j,nz-1)&
            !              -uz(i,j,nz-1)+uz(i,j,nz-2))&
            !              +bskzt*(uz(i,j,nz-1)-uz(i,j,nz-1)&
            !              -uz(i,j,nz-1)+uz(i,j,nz-3))&
            !              +cskzt*(uz(i,j,nz-2)-uz(i,j,nz-1)&
            !              -uz(i,j,nz-1)+uz(i,j,nz-4))
            !         tz(i,j,nz  )=askzt*(uz(i,j,nz-1)-uz(i,j,nz  )&
            !              -uz(i,j,nz  )+uz(i,j,nz-1))&
            !              +bskzt*(uz(i,j,nz-2)-uz(i,j,nz  )&
            !              -uz(i,j,nz  )+uz(i,j,nz-2))&
            !              +cskzt*(uz(i,j,nz-3)-uz(i,j,nz  )&
            !              -uz(i,j,nz  )+uz(i,j,nz-3))
            !      enddo
            !      enddo
            !      do k=2,nz
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
            !      enddo
            !      enddo
            !      enddo
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
            !      enddo
            !      enddo
            !      do k=nz-1,1,-1
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
            !      enddo
            !      enddo
            !      enddo
            !   endif
            !   if (npaire==0) then
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,1)=zero
            !         tz(i,j,2)=askzt*(uz(i,j,3)-uz(i,j,2)&
            !              -uz(i,j,2)+uz(i,j,1))&
            !              +bskzt*(uz(i,j,4)-uz(i,j,2)&
            !              -uz(i,j,2)-uz(i,j,2))&
            !              +cskzt*(uz(i,j,5)-uz(i,j,2)&
            !              -uz(i,j,2)-uz(i,j,3))
            !         tz(i,j,3)=askzt*(uz(i,j,4)-uz(i,j,3)&
            !              -uz(i,j,3)+uz(i,j,2))&
            !              +bskzt*(uz(i,j,5)-uz(i,j,3)&
            !              -uz(i,j,3)+uz(i,j,1))&
            !              +cskzt*(uz(i,j,6)-uz(i,j,3)&
            !              -uz(i,j,3)-uz(i,j,2))
            !      enddo
            !      enddo
            !      do k=4,nz-3
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,k)=askzt*(uz(i,j,k+1)-uz(i,j,k  )&
            !              -uz(i,j,k  )+uz(i,j,k-1))&
            !              +bskzt*(uz(i,j,k+2)-uz(i,j,k  )&
            !              -uz(i,j,k  )+uz(i,j,k-2))&
            !              +cskzt*(uz(i,j,k+3)-uz(i,j,k  )&
            !              -uz(i,j,k  )+uz(i,j,k-3))
            !      enddo
            !      enddo
            !      enddo
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,nz-2)=askzt*( uz(i,j,nz-1)-uz(i,j,nz-2)&
            !              -uz(i,j,nz-2)+uz(i,j,nz-3))&
            !              +bskzt*( uz(i,j,nz  )-uz(i,j,nz-2)&
            !              -uz(i,j,nz-2)+uz(i,j,nz-4))&
            !              +cskzt*(-uz(i,j,nz-1)-uz(i,j,nz-2)&
            !              -uz(i,j,nz-2)+uz(i,j,nz-5))
            !         tz(i,j,nz-1)=askzt*( uz(i,j,nz  )-uz(i,j,nz-1)&
            !              -uz(i,j,nz-1)+uz(i,j,nz-2))&
            !              +bskzt*(-uz(i,j,nz-1)-uz(i,j,nz-1)&
            !              -uz(i,j,nz-1)+uz(i,j,nz-3))&
            !              +cskzt*(-uz(i,j,nz-2)-uz(i,j,nz-1)&
            !              -uz(i,j,nz-1)+uz(i,j,nz-4))
            !         tz(i,j,nz  )=zero
            !      enddo
            !      enddo
            !      do k=2,nz
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
            !      enddo
            !      enddo
            !      enddo
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
            !      enddo
            !      enddo
            !      do k=nz-1,1,-1
            !      do j=1,ny
            !      do i=1,nx
            !         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
            !      enddo
            !      enddo
            !      enddo
            !   endif
            !endif

            !if (nclz==2) then
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,1)=as1z*uz(i,j,1)+bs1z*uz(i,j,2)&
            !           +cs1z*uz(i,j,3)+ds1z*uz(i,j,4)
            !      tz(i,j,2)=as2z*(uz(i,j,3)-uz(i,j,2)&
            !           -uz(i,j,2)+uz(i,j,1))
            !      tz(i,j,3)=as3z*(uz(i,j,4)-uz(i,j,3)&
            !           -uz(i,j,3)+uz(i,j,2))&
            !           +bs3z*(uz(i,j,5)-uz(i,j,3)&
            !           -uz(i,j,3)+uz(i,j,1))
            !   enddo
            !   enddo
            !   do k=4,nz-3
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,k)=askzt*(uz(i,j,k+1)-uz(i,j,k  )&
            !           -uz(i,j,k  )+uz(i,j,k-1))&
            !           +bskzt*(uz(i,j,k+2)-uz(i,j,k  )&
            !           -uz(i,j,k  )+uz(i,j,k-2))&
            !           +cskzt*(uz(i,j,k+3)-uz(i,j,k  )&
            !           -uz(i,j,k  )+uz(i,j,k-3))
            !   enddo
            !   enddo
            !   enddo
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,nz-2)=astz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
            !           -uz(i,j,nz-2)+uz(i,j,nz-3))&
            !           +bstz*(uz(i,j,nz  )-uz(i,j,nz-2)&
            !           -uz(i,j,nz-2)+uz(i,j,nz-4))
            !      tz(i,j,nz-1)=asmz*(uz(i,j,nz  )-uz(i,j,nz-1)&
            !           -uz(i,j,nz-1)+uz(i,j,nz-2))
            !      tz(i,j,nz  )=asnz*uz(i,j,nz  )+bsnz*uz(i,j,nz-1)&
            !           +csnz*uz(i,j,nz-2)+dsnz*uz(i,j,nz-3)
            !   enddo
            !   enddo
            !   do k=2,nz
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
            !   enddo
            !   enddo
            !   enddo
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,nz)=tz(i,j,nz)*swz(nz)
            !   enddo
            !   enddo
            !   do k=nz-1,1,-1
            !   do j=1,ny
            !   do i=1,nx
            !      tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
            !   enddo
            !   enddo
            !   enddo
            !endif

            !return
            !end subroutine derzzt

            !********************************************************************
            !
            subroutine scalar_schemes(fpi2t)
              !
              !********************************************************************
              use decomp_2d, only : mytype
              use param
              use derivX
              use derivY
              use derivZ
              use variables
              use ludecomp

              implicit none

              real(mytype), intent(in) :: fpi2t
              integer :: i,j,k
              real(mytype), dimension(ny) :: aamt10,bbmt10,ccmt10,ddmt10,eemt10,ggmt10,hhmt10,wwmt10,zzmt10
              real(mytype), dimension(ny) :: rrmt10,qqmt10,vvmt10,ssmt10
              real(mytype), dimension(ny) :: aamt11,bbmt11,ccmt11,ddmt11,eemt11,ggmt11,hhmt11,wwmt11,zzmt11
              real(mytype), dimension(ny) :: rrmt11,qqmt11,vvmt11,ssmt11

              xcst_pr=xnu*dt*0.5/nsc(1)

              alsaixt=(45.*fpi2t*pi*pi-272.)/(2.*(45.*fpi2t*pi*pi-208.))
              asixt  =((6.-9.*alsaixt)/4.)/dx2
              bsixt  =((-3.+24.*alsaixt)/5.)/(4.*dx2)
              csixt  =((2.-11.*alsaixt)/20.)/(9.*dx2)

              alsajyt=(45.*fpi2t*pi*pi-272.)/(2.*(45.*fpi2t*pi*pi-208.))
              asjyt  =((6.-9.*alsajyt)/4.)/dy2
              bsjyt  =((-3.+24.*alsajyt)/5.)/(4.*dy2)
              csjyt  =((2.-11.*alsajyt)/20.)/(9.*dy2)

              alsakzt=(45.*fpi2t*pi*pi-272.)/(2.*(45.*fpi2t*pi*pi-208.))
              askzt  =((6.-9.*alsakzt)/4.)/dz2
              bskzt  =((-3.+24.*alsakzt)/5.)/(4.*dz2)
              cskzt  =((2.-11.*alsakzt)/20.)/(9.*dz2)

              if (nclx.eq.0) then
                sfxt(1)   =alsaixt
                sfxt(2)   =alsaixt
                sfxt(nx-2)=alsaixt
                sfxt(nx-1)=alsaixt
                sfxt(nx)  =zero
                scxt(1)   =two
                scxt(2)   =one
                scxt(nx-2)=one
                scxt(nx-1)=one
                scxt(nx  )=one+alsaixt*alsaixt
                sbxt(1)   =alsaixt
                sbxt(2)   =alsaixt
                sbxt(nx-2)=alsaixt
                sbxt(nx-1)=alsaixt
                sbxt(nx  )=zero
                do i=3,nx-3
                  sfxt(i)=alsaixt
                  scxt(i)=one
                  sbxt(i)=alsaixt
                enddo
              endif

              if (nclx.eq.1) then
                sfxt(1)   =alsaixt+alsaixt
                sfxt(2)   =alsaixt
                sfxt(nx-2)=alsaixt
                sfxt(nx-1)=alsaixt
                sfxt(nx)  =zero
                scxt(1)   =one
                scxt(2)   =one
                scxt(nx-2)=one
                scxt(nx-1)=one
                scxt(nx  )=one
                sbxt(1)   =alsaixt
                sbxt(2)   =alsaixt
                sbxt(nx-2)=alsaixt
                sbxt(nx-1)=alsaixt+alsaixt
                sbxt(nx  )=zero
                do i=3,nx-3
                  sfxt(i)=alsaixt
                  scxt(i)=one
                  sbxt(i)=alsaixt
                enddo
              endif

              if (nclx.eq.2) then
                sfxt(1)   =alsa1x
                sfxt(2)   =alsa2x
                sfxt(3)   =alsa3x
                sfxt(nx-3)=alsaixt
                sfxt(nx-2)=alsatx
                sfxt(nx-1)=alsamx
                sfxt(nx)  =zero
                scxt(1)   =one
                scxt(2)   =one
                scxt(3)   =one
                scxt(nx-3)=one
                scxt(nx-2)=one
                scxt(nx-1)=one
                scxt(nx  )=one
                sbxt(1)   =alsa2x
                sbxt(2)   =alsa3x
                sbxt(3)   =alsaixt
                sbxt(nx-3)=alsatx
                sbxt(nx-2)=alsamx
                sbxt(nx-1)=alsanx
                sbxt(nx  )=zero
                do i=4,nx-4
                  sfxt(i)=alsaixt
                  scxt(i)=one
                  sbxt(i)=alsaixt
                enddo
              endif

              if (ncly.eq.0) then
                sfyt(1)   =alsajyt
                sfyt(2)   =alsajyt
                sfyt(ny-2)=alsajyt
                sfyt(ny-1)=alsajyt
                sfyt(ny)  =zero
                scyt(1)   =two
                scyt(2)   =one
                scyt(ny-2)=one
                scyt(ny-1)=one
                scyt(ny  )=one+alsajyt*alsajyt
                sbyt(1)   =alsajyt
                sbyt(2)   =alsajyt
                sbyt(ny-2)=alsajyt
                sbyt(ny-1)=alsajyt
                sbyt(ny  )=zero
                do j=3,ny-3
                  sfyt(j)=alsajyt
                  scyt(j)=one
                  sbyt(j)=alsajyt
                enddo
              endif

              if (ncly.eq.1) then
                sfyt(1)   =alsajyt+alsajyt
                sfyt(2)   =alsajyt
                sfyt(ny-2)=alsajyt
                sfyt(ny-1)=alsajyt
                sfyt(ny)  =zero
                scyt(1)   =one
                scyt(2)   =one
                scyt(ny-2)=one
                scyt(ny-1)=one
                scyt(ny  )=one
                sbyt(1)   =alsajyt
                sbyt(2)   =alsajyt
                sbyt(ny-2)=alsajyt
                sbyt(ny-1)=alsajyt+alsajyt
                sbyt(ny  )=zero
                do j=3,ny-3
                  sfyt(j)=alsajyt
                  scyt(j)=one
                  sbyt(j)=alsajyt
                enddo
              endif

              if (ncly.eq.2) then
                sfyt(1)   =alsa1y
                sfyt(2)   =alsa2y
                sfyt(3)   =alsa3y
                sfyt(ny-3)=alsajyt
                sfyt(ny-2)=alsaty
                sfyt(ny-1)=alsamy
                sfyt(ny)  =zero
                scyt(1)   =one
                scyt(2)   =one
                scyt(3)   =one
                scyt(ny-3)=one
                scyt(ny-2)=one
                scyt(ny-1)=one
                scyt(ny  )=one
                sbyt(1)   =alsa2y
                sbyt(2)   =alsa3y
                sbyt(3)   =alsajyt
                sbyt(ny-3)=alsaty
                sbyt(ny-2)=alsamy
                sbyt(ny-1)=alsany
                sbyt(ny  )=zero
                do j=4,ny-4
                  sfyt(j)=alsajyt
                  scyt(j)=one
                  sbyt(j)=alsajyt
                enddo
              endif

              if (nclz.eq.0) then
                sfzt(1)   =alsakzt
                sfzt(2)   =alsakzt
                sfzt(nz-2)=alsakzt
                sfzt(nz-1)=alsakzt
                sfzt(nz)  =zero
                sczt(1)   =two
                sczt(2)   =one
                sczt(nz-2)=one
                sczt(nz-1)=one
                sczt(nz  )=one+alsakzt*alsakzt
                sbzt(1)   =alsakzt
                sbzt(2)   =alsakzt
                sbzt(nz-2)=alsakzt
                sbzt(nz-1)=alsakzt
                sbzt(nz  )=zero
                do k=3,nz-3
                  sfzt(k)=alsakzt
                  sczt(k)=one
                  sbzt(k)=alsakzt
                enddo
              endif

              if (nclz.eq.1) then
                sfzt(1)   =alsakzt+alsakzt
                sfzt(2)   =alsakzt
                sfzt(nz-2)=alsakzt
                sfzt(nz-1)=alsakzt
                sfzt(nz)  =zero
                sczt(1)   =one
                sczt(2)   =one
                sczt(nz-2)=one
                sczt(nz-1)=one
                sczt(nz  )=one
                sbzt(1)   =alsakzt
                sbzt(2)   =alsakzt
                sbzt(nz-2)=alsakzt
                sbzt(nz-1)=alsakzt+alsakzt
                sbzt(nz  )=zero
                do k=3,nz-3
                  sfzt(k)=alsakzt
                  sczt(k)=one
                  sbzt(k)=alsakzt
                enddo
              endif

              if (nclz.eq.2) then
                sfzt(1)   =alsa1z
                sfzt(2)   =alsa2z
                sfzt(3)   =alsa3z
                sfzt(nz-3)=alsakzt
                sfzt(nz-2)=alsatz
                sfzt(nz-1)=alsamz
                sfzt(nz)  =zero
                sczt(1)   =one
                sczt(2)   =one
                sczt(3)   =one
                sczt(nz-3)=one
                sczt(nz-2)=one
                sczt(nz-1)=one
                sczt(nz  )=one
                sbzt(1)   =alsa2z
                sbzt(2)   =alsa3z
                sbzt(3)   =alsakzt
                sbzt(nz-3)=alsatz
                sbzt(nz-2)=alsamz
                sbzt(nz-1)=alsanz
                sbzt(nz  )=zero
                do k=4,nz-4
                  sfzt(k)=alsakzt
                  sczt(k)=one
                  sbzt(k)=alsakzt
                enddo
              endif

              if (nclx.eq.1) then
                sfxt(1)=zero
              endif
              if (ncly.eq.1) then
                sfyt(1)=zero
              endif
              if (nclz.eq.1) then
                sfzt(1)=zero
              endif

              call prepare (sbxt,scxt,sfxt,ssxt,swxt,nx)
              call prepare (sbyt,scyt,sfyt,ssyt,swyt,ny)
              call prepare (sbzt,sczt,sfzt,sszt,swzt,nz)

              call prepare (sbxt,scxt,sfxpt,ssxpt,swxpt,nx)
              call prepare (sbyt,scyt,sfypt,ssypt,swypt,ny)
              call prepare (sbzt,sczt,sfzpt,sszpt,swzpt,nz)

              if (nclx.eq.1) then
                sbxt(nx-1)=zero
                call prepare (sbxt,scxt,sfxt,ssxt,swxt,nx)
              endif
              if (ncly.eq.1) then
                sbyt(ny-1)=zero
                call prepare (sbyt,scyt,sfyt,ssyt,swyt,ny)
              endif
              if (nclz.eq.1) then
                sbzt(nz-1)=zero
                call prepare (sbzt,sczt,sfzt,sszt,swzt,nz)
              endif

              !!!!!!!!!!!!!!!!!!!!!!
              !CAS SEPTADIAGONALE !!
              !!!!!!!!!!!!!!!!!!!!!!

              !!! NCL = 2, dirichlet imposé, fonction nulle à la paroi
              !
              !DIAG
              aamt(1     )=as1y
              aamt(ny    )=asny
              aamt(2     )=-two*as2y
              aamt(ny-1  )=-two*asmy
              aamt(3     )=-two*(as3y+bs3y)
              aamt(ny-2  )=-two*(asty+bsty)
              aamt(4:ny-3)=-two*(asjyt+bsjyt+csjyt)
              if (istret==0) then
                aamt = one-xcst_pr*aamt
              else
                aamt = one/pp2y-xcst_pr*aamt
              endif
              !CL sur aamt
              if (istret==0) then
                aamt(1 )=alpha_0+beta_0*(11./6./dy)
                aamt(ny)=alpha_n+beta_n*(11./6./dy)
              else
                aamt(1 )=alpha_0+beta_0*ppy(1 )*(11./6./dy)
                aamt(ny)=alpha_n+beta_n*ppy(ny)*(11./6./dy)
              endif
              !
              !DIAG SUP 1
              bbmt(1     )=bs1y
              bbmt(ny    )=bsny
              bbmt(2     )=as2y
              bbmt(ny-1  )=asmy
              bbmt(3     )=as3y
              bbmt(ny-2  )=asty
              bbmt(4:ny-3)=asjyt
              bbmt = -xcst_pr*bbmt
              if (istret==0) then
                bbmt(2     )=bbmt(2     )+alsa2y
                bbmt(ny-1  )=bbmt(ny-1  )+alsamy
                bbmt(3     )=bbmt(3     )+alsa3y
                bbmt(ny-2  )=bbmt(ny-2  )+alsaty
                bbmt(4:ny-3)=bbmt(4:ny-3)+alsajyt
              else
                bbmt(2     )=bbmt(2     )+alsa2y/pp2y(3)
                bbmt(ny-1  )=bbmt(ny-1  )+alsamy/pp2y(ny)
                bbmt(3     )=bbmt(3     )+alsa3y/pp2y(4)
                bbmt(ny-2  )=bbmt(ny-2  )+alsaty/pp2y(ny-1)
                bbmt(4:ny-3)=bbmt(4:ny-3)+alsajyt/pp2y(5:ny-2)
              endif
              !CL sur bbmt
              if (istret==0) then
                bbmt(1 )=beta_0*(-18./6./dy)
                bbmt(ny)=zero
              else
                bbmt(1 )=beta_0*ppy(1)*(-18./6./dy)
                bbmt(ny)=zero
              endif
              !
              !DIAG SUP 2
              ccmt(1     )=cs1y
              ccmt(ny    )=csny
              ccmt(2     )=zero
              ccmt(ny-1  )=zero
              ccmt(3     )=bs3y
              ccmt(ny-2  )=bsty
              ccmt(4:ny-3)=bsjyt
              ccmt = -xcst_pr*ccmt
              !CL sur ccmt
              if (istret==0) then
                ccmt(1 )=beta_0*(9./6./dy)
                ccmt(ny)=zero
              else
                ccmt(1 )=beta_0*ppy(1)*(9./6./dy)
                ccmt(ny)=zero
              endif
              !
              !DIAG SUP 3
              rrmt(1     )=ds1y
              rrmt(ny    )=dsny
              rrmt(2     )=zero
              rrmt(ny-1  )=zero
              rrmt(3     )=zero
              rrmt(ny-2  )=zero
              rrmt(4:ny-3)=csjyt
              rrmt = -xcst_pr*rrmt
              !CL sur rrmt
              if (istret==0) then
                rrmt(1 )=beta_0*(-2./6./dy)
                rrmt(ny)=zero
              else
                rrmt(1 )=beta_0*ppy(1)*(-2./6./dy)
                rrmt(ny)=zero
              endif
              !
              !DIAG INF 1
              if (istret==0) then
                ddmt=bbmt
              else
                ddmt(1     )=bs1y
                ddmt(ny    )=bsny
                ddmt(2     )=as2y
                ddmt(ny-1  )=asmy
                ddmt(3     )=as3y
                ddmt(ny-2  )=asty
                ddmt(4:ny-3)=asjyt
                ddmt = -xcst_pr*ddmt
                ddmt(2     )=ddmt(2     )+alsa2y/pp2y(1)
                ddmt(ny-1  )=ddmt(ny-1  )+alsamy/pp2y(ny-2)
                ddmt(3     )=ddmt(3     )+alsa3y/pp2y(2)
                ddmt(ny-2  )=ddmt(ny-2  )+alsaty/pp2y(ny-3)
                ddmt(4:ny-3)=ddmt(4:ny-3)+alsajyt/pp2y(3:ny-4)
              endif
              !CL sur ddmt
              if (istret==0) then
                ddmt(1 )=zero
                ddmt(ny)=beta_n*(-18./6./dy)
              else
                ddmt(1 )=zero
                ddmt(ny)=beta_n*ppy(ny)*(-18./6./dy)
              endif
              !
              !DIAG INF 2
              eemt=ccmt
              !CL sur eemt
              if (istret==0) then
                eemt(1)=zero
                eemt(ny)=beta_n*(9./6./dy)
              else
                eemt(1)=zero
                eemt(ny)=beta_n*ppy(ny)*(9./6./dy)
              endif
              !
              !DIAG INF 3
              qqmt=rrmt
              !CL sur qqmt
              if (istret==0) then
                qqmt(1)=zero
                qqmt(ny)=beta_n*(-two/6./dy)
              else
                qqmt(1)=zero
                qqmt(ny)=beta_n*ppy(ny)*(-two/6./dy)
              endif

              !!! NCL = 1, npaire=0, dirichlet imposé, fonction impaire
              !
              ! DIAG
              aamt10(1     )=zero
              aamt10(ny    )=zero
              aamt10(2     )=-two*asjyt-3.*bsjyt-two*csjyt
              aamt10(ny-1  )=aamt10(2)
              aamt10(3:ny-2)=-two*(asjyt+bsjyt+csjyt)
              if (istret==0) then
                aamt10 = one - xcst_pr*aamt10
              else
                aamt10 = one/pp2y - xcst_pr*aamt10
              endif
              !
              !DIAG SUP 1
              bbmt10(1     )=zero
              bbmt10(ny    )=zero
              bbmt10(2     )=asjyt-csjyt
              bbmt10(ny-1  )=asjyt
              bbmt10(3     )=asjyt
              bbmt10(ny-2  )=asjyt-csjyt
              bbmt10(4:ny-3)=asjyt
              if (istret==0) then
                bbmt10 = alsajyt - xcst_pr*bbmt10
              else
                bbmt10(2:ny-1) = alsajyt/pp2y(3:ny) - xcst_pr*bbmt10(2:ny-1)
              endif
              !CL sur bbm10t
              bbmt10(1 )=zero
              bbmt10(ny)=zero
              !
              !DIAG SUP 2
              ccmt10(1     )=zero
              ccmt10(ny    )=zero
              ccmt10(2     )=bsjyt
              ccmt10(ny-1  )=zero
              ccmt10(3     )=bsjyt
              ccmt10(ny-2  )=bsjyt
              ccmt10(4:ny-3)=bsjyt
              ccmt10 = -xcst_pr*ccmt10
              !
              !DIAG SUP 3
              rrmt10(1     )=zero
              rrmt10(ny    )=zero
              rrmt10(2     )=csjyt
              rrmt10(ny-1  )=zero
              rrmt10(3     )=csjyt
              rrmt10(ny-2  )=zero
              rrmt10(4:ny-3)=csjyt
              rrmt10 = -xcst_pr*rrmt10
              !
              !DIAG INF 1
              ddmt10(1     )=zero
              ddmt10(ny    )=zero
              ddmt10(2     )=asjyt
              ddmt10(ny-1  )=asjyt-csjyt
              ddmt10(3     )=asjyt-csjyt
              ddmt10(ny-2  )=asjyt
              ddmt10(4:ny-3)=asjyt
              if (istret==0) then
                ddmt10 = alsajyt - xcst_pr*ddmt10
              else
                ddmt10(2:ny-1) = alsajyt/pp2y(1:ny-2) - xcst_pr*ddmt10(2:ny-1)
              endif
              !CL sur ddmt10
              ddmt10(1 )=zero
              ddmt10(ny)=zero
              !
              !DIAG INF 2
              eemt10(1     )=zero
              eemt10(ny    )=zero
              eemt10(2     )=zero
              eemt10(ny-1  )=bsjyt
              eemt10(3     )=bsjyt
              eemt10(ny-2  )=bsjyt
              eemt10(4:ny-3)=bsjyt
              eemt10 = -xcst_pr*eemt10
              !
              !DIAG INF 3
              qqmt10(1     )=zero
              qqmt10(ny    )=zero
              qqmt10(2     )=zero
              qqmt10(ny-1  )=csjyt
              qqmt10(3     )=zero
              qqmt10(ny-2  )=csjyt
              qqmt10(4:ny-3)=csjyt
              qqmt10 = -xcst_pr*qqmt10

              !!! NCL = 1, npaire=1, neumann imposé, fonction paire
              !
              ! DIAG
              aamt11(1     )=-two*(asjyt+bsjyt+csjyt)
              aamt11(ny    )=aamt11(1)
              aamt11(2     )=-two*asjyt-bsjyt-two*csjyt
              aamt11(ny-1  )=aamt11(2)
              aamt11(3:ny-2)=-two*(asjyt+bsjyt+csjyt)
              if (istret==0) then
                aamt11 = one - xcst_pr*aamt11
              else
                aamt11 = one/pp2y - xcst_pr*aamt11
              endif
              !
              !DIAG SUP 1
              bbmt11(1     )=two*asjyt
              bbmt11(ny    )=zero
              bbmt11(2     )=asjyt+csjyt
              bbmt11(ny-1  )=asjyt
              bbmt11(3     )=asjyt
              bbmt11(ny-2  )=asjyt+csjyt
              bbmt11(4:ny-3)=asjyt
              if (istret==0) then
                bbmt11 = alsajyt - xcst_pr*bbmt11
              else
                bbmt11(1:ny-1) = alsajyt/pp2y(2:ny) - xcst_pr*bbmt11(1:ny-1)
              endif
              !CL sur bbm11t
              if (istret==0) then
                bbmt11(1 )=bbmt11(1)+alsajyt
              else
                bbmt11(1 )=bbmt11(1)+alsajyt/pp2y(2)
              endif
              bbmt11(ny)=zero
              !
              !DIAG SUP 2
              ccmt11(1     )=two*bsjyt
              ccmt11(ny    )=zero
              ccmt11(2     )=bsjyt
              ccmt11(ny-1  )=zero
              ccmt11(3     )=bsjyt
              ccmt11(ny-2  )=bsjyt
              ccmt11(4:ny-3)=bsjyt
              ccmt11 = -xcst_pr*ccmt11
              !
              !DIAG SUP 3
              rrmt11(1     )=two*csjyt
              rrmt11(ny    )=zero
              rrmt11(2     )=csjyt
              rrmt11(ny-1  )=zero
              rrmt11(3     )=csjyt
              rrmt11(ny-2  )=zero
              rrmt11(4:ny-3)=csjyt
              rrmt11 = -xcst_pr*rrmt11
              !
              !DIAG INF 1
              ddmt11(1     )=zero
              ddmt11(ny    )=two*asjyt
              ddmt11(2     )=asjyt
              ddmt11(ny-1  )=asjyt+csjyt
              ddmt11(3     )=asjyt+csjyt
              ddmt11(ny-2  )=asjyt
              ddmt11(4:ny-3)=asjyt
              if (istret==0) then
                ddmt11 = alsajyt - xcst_pr*ddmt11
              else
                ddmt11(2:ny) = alsajyt/pp2y(1:ny-1) - xcst_pr*ddmt11(2:ny)
              endif
              !CL sur ddmt11
              ddmt11(1 )=zero
              if (istret==0) then
                ddmt11(ny)=ddmt11(ny)+alsajyt!a1
              else
                ddmt11(ny)=ddmt11(ny)+alsajyt/pp2y(ny-1)!a1
              endif
              !
              !DIAG INF 2
              eemt11(1     )=zero
              eemt11(ny    )=two*bsjyt
              eemt11(2     )=zero
              eemt11(ny-1  )=bsjyt
              eemt11(3     )=bsjyt
              eemt11(ny-2  )=bsjyt
              eemt11(4:ny-3)=bsjyt
              eemt11 = -xcst_pr*eemt11
              !
              !DIAG INF 3
              qqmt11(1     )=zero
              qqmt11(ny    )=two*csjyt
              qqmt11(2     )=zero
              qqmt11(ny-1  )=csjyt
              qqmt11(3     )=zero
              qqmt11(ny-2  )=csjyt
              qqmt11(4:ny-3)=csjyt
              qqmt11 = -xcst_pr*qqmt11

              !!! NXL = 0
              !DIAG
              if (istret==0) then
                aamt0 = one-xcst_pr*(-two*(asjyt+bsjyt+csjyt))
              else
                aamt0 = one/pp2y-xcst_pr*(-two*(asjyt+bsjyt+csjyt))
              endif
              !
              !DIAG SUP 1
              if (istret==0) then
                bbmt0 = alsajyt-xcst_pr*asjyt
              else
                bbmt0(1:ny-1) = alsajyt/pp2y(2:ny)-xcst_pr*asjyt
                bbmt0(ny) = alsajyt/pp2y(1)-xcst_pr*asjyt
              endif
              !
              !DIAG SUP 2
              ccmt0 = -xcst_pr*bsjyt
              !
              !DIAG SUP 3
              rrmt0 = -xcst_pr*csjyt
              !
              !DIAG INF 1
              if (istret==0) then
                ddmt0=bbmt0
              else
                ddmt0(1)=alsajyt/pp2y(ny) -xcst_pr*asjyt
                ddmt0(2:ny)=alsajyt/pp2y(1:ny-1) -xcst_pr*asjyt
              endif
              !
              !DIAG INF 2
              eemt0=ccmt0
              !
              !DIAG INF 3
              qqmt0=rrmt0

              call ludecomp7(aamt,bbmt,ccmt,ddmt,eemt,qqmt,ggmt,hhmt,ssmt,rrmt,&
              vvmt,wwmt,zzmt,ny)
              call ludecomp7(aamt10,bbmt10,ccmt10,ddmt10,eemt10,qqmt10,ggmt10,hhmt10,ssmt10,rrmt10,&
              vvmt10,wwmt10,zzmt10,ny)
              call ludecomp7(aamt11,bbmt11,ccmt11,ddmt11,eemt11,qqmt11,ggmt11,hhmt11,ssmt11,rrmt11,&
              vvmt11,wwmt11,zzmt11,ny)
              call ludecomp7(aamt0,bbmt0,ccmt0,ddmt0,eemt0,qqmt0,ggmt0,hhmt0,ssmt0,rrmt0,&
              vvmt0,wwmt0,zzmt0,l1mt,l2mt,l3mt,u1mt,u2mt,u3mt,ny)

              ! npaire=0
              aamt1=aamt10; bbmt1=bbmt10; ccmt1=ccmt10; ddmt1=ddmt10
              eemt1=eemt10; ggmt1=ggmt10; hhmt1=hhmt10; wwmt1=wwmt10; zzmt1=zzmt10
              rrmt1=rrmt10; qqmt1=qqmt10; vvmt1=vvmt10; ssmt1=ssmt10
              !! npaire=1
              !   aamt1=aamt11; bbmt1=bbmt11; ccmt1=ccmt11; ddmt1=ddmt11
              !   eemt1=eemt11; ggmt1=ggmt11; hhmt1=hhmt11; wwmt1=wwmt11; zzmt1=zzmt11
              !   rrmt1=rrmt11; qqmt1=qqmt11; vvmt1=vvmt11; ssmt1=ssmt11

            end subroutine scalar_schemes
#endif
