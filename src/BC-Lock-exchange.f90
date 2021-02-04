!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Bartholomew P., Deskos G., Frantz R.A.S., Schuch F.N., Lamballais E. &
!    Laizet S., 2020, Xcompact3D: An open-source framework for solving
!    turbulence problems on a Cartesian mesh, SoftwareX, vol 12, pp 100550
!
!    2-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    3-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

module lockexch

  use decomp_2d, only : mytype, real_type, real2_type
  use decomp_2d, only : xsize, ysize, zsize
  use decomp_2d, only : xstart, ystart, zstart
  use decomp_2d, only : xend, yend, zend
  use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x

  use variables, only : numscalar

  use var, only : xnu, ri, uset, sc, Fr, prandtl
  use var, only : gravy
  use var, only : nrank
  use var, only : dens1, dens2

  use var, only : zero, half, one, two, five, twelve, thirteen
  use var, only : dx, dy, dz, nx, ny, nz

  use var, only : nrhotime

  use param, only : ilmn, ibirman_eos
  use param, only : itime, ioutput, iprocessing
  use param, only : t

  implicit none

  real(mytype), save :: pfront
  real(mytype), save, allocatable :: area2(:,:), vol1(:,:,:)

  integer :: FS
  character(len=100) :: fileformat
  integer, parameter :: filenum = 67
  character(len=1),parameter :: NL=char(10) !new line character

  logical, save :: init = .FALSE.

  private
  public :: init_lockexch, boundary_conditions_lockexch, postprocess_lockexch, &
       pfront, set_fluid_properties_lockexch

contains

  subroutine boundary_conditions_lockexch (rho1, phi1)

    USE param
    USE variables
    USE decomp_2d
    USE MPI
    use sandbox, only : deposit

    implicit none

    integer  :: i,j,k,is
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    if (xstart(2).eq.1 .and. ilmn) then
       j = 1
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             rho1(i, j, k, 1) = rho1(i, j + 1, k, 1) !! drho/dy=0
          enddo
       enddo
    endif

    if (iscalar .eq. 1) call deposit(phi1)

    return
  end subroutine boundary_conditions_lockexch

  subroutine init_lockexch (rho1,ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    USE var, only : mu1

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: um,x,y
    integer :: k,j,i,ii,is,it,code

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             x=real(i+xstart(1)-1-1,mytype)*dx-pfront
             do is=1,numscalar
                phi1(i,j,k,is)=half * (one - tanh((sc(is) / xnu)**(half) * x)) * cp(is)
             enddo
          enddo
       enddo
    enddo
    do k = 1,xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             do is=1,numscalar
                if (phi1(i,j,k,is).gt.cp(is)) then
                   phi1(i,j,k,is) = cp(is)
                elseif (phi1(i,j,k,is).lt.zero) then
                   phi1(i,j,k,is) = zero
                endif
             enddo
          enddo
       enddo
    enddo

    if (ilmn) then
      call set_fluid_properties_lockexch(rho1, mu1)
      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               x=real(i+xstart(1)-1-1,mytype)*dx-pfront
               rho1(i,j,k,1) = half * (one - tanh((prandtl / xnu)**half * x)) &
                    * (dens1 - dens2) + dens2
            enddo
         enddo
      enddo
      do k = 1,xsize(3)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               if (rho1(i,j,k,1).gt.max(dens1, dens2)) then
                  rho1(i,j,k,1) = max(dens1, dens2)
               elseif (rho1(i,j,k,1).lt.min(dens1, dens2)) then
                  rho1(i,j,k,1) = min(dens1, dens2)
               endif
            enddo
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
               um=exp(-twentyfive*x*x)*init_noise
               ux1(i,j,k)=um*(two*ux1(i,j,k)-one)
               uy1(i,j,k)=um*(two*uy1(i,j,k)-one)
               uz1(i,j,k)=um*(two*uz1(i,j,k)-one)
            enddo
         enddo
      enddo
    endif

    return
  end subroutine init_lockexch


  subroutine postprocess_lockexch(rho1,ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

    use decomp_2d, only : alloc_x
    use param, only : iscalar
    use visu, only : filenamedigits, ifilenameformat, &
        write_xdmf_header, write_field, write_xdmf_footer
    use var, only : ux2,uy2,uz2,phi2,rho2
    use var, only : ux3,uy3,uz3,phi3,rho3
    use tools, only : mean_plane_z
    use sandbox, only : postprocessing_aux, budget

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phisum1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: diss1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: diss2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: diss3
    real(mytype),dimension(ysize(1),ysize(3),numscalar) :: dep2
    real(mytype),dimension(zsize(1),zsize(2),numscalar) :: phim3
    real(mytype),dimension(zsize(1),zsize(2)) :: rhom3
    integer :: i,j,k,is
    character(len=30) :: num

    real(mytype) :: mp(numscalar),dms(numscalar),xf(1:2,1:3),xf2d(1:2,1:2)

    if (.not.init) then
       call alloc_x(vol1, opt_global=.true.)
       do k=xstart(3),xend(3)
          do j=xstart(2),xend(2)
             do i=xstart(1),xend(1)
                vol1(i,j,k)=dx*dy*dz
                if (i .eq. 1 .or. i .eq. nx) vol1(i,j,k) = vol1(i,j,k) * half !five/twelve
                if (j .eq. 1 .or. j .eq. ny) vol1(i,j,k) = vol1(i,j,k) * half !five/twelve
                if (k .eq. 1 .or. k .eq. nz) vol1(i,j,k) = vol1(i,j,k) * half !five/twelve
                ! if (i .eq. 2 .or. i .eq. nx-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
                ! if (j .eq. 2 .or. j .eq. ny-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
                ! if (k .eq. 2 .or. k .eq. nz-1) vol1(i,j,k) = vol1(i,j,k) * thirteen/twelve
             end do
          end do
       end do

       allocate(area2(ystart(1):yend(1),ystart(3):yend(3)))
       area2=dx*dz
       do k=ystart(3),yend(3)
          do i=ystart(1),yend(1)
             if (i .eq. 1 .or. i .eq. nx) area2(i,k) = area2(i,k)/two
             if (k .eq. 1 .or. k .eq. nz)  area2(i,k) = area2(i,k)/two
          end do
       end do

       init = .TRUE.
    endif

    if (mod(itime,iprocessing).ne.0) return

    mp=zero; dms=zero; xf=zero; xf2d=zero

    phisum1(:,:,:)=zero
    do is=1, numscalar
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

    do is = 2, numscalar
       phim3(:,:,1) = phim3(:,:,1) + phim3(:,:,is)
    enddo

    if (ilmn) then
       call transpose_x_to_y(rho1(:,:,:,1), rho2)
       call transpose_y_to_z(rho2, rho3)
       call mean_plane_z(rho3, zsize(1), zsize(2), zsize(3), rhom3)

       if (ibirman_eos) then
          phisum1(:,:,:) = phisum1(:,:,:) + (rho1(:,:,:,1) - dens2) / (dens1 - dens2)
          if (numscalar.gt.0) then
             do j = 1, zsize(2)
                do i = 1, zsize(1)
                   phim3(i,j,1) = phim3(i,j,1) + (rhom3(i,j) - dens2) / (dens1 - dens2)
                enddo
             enddo
          endif
       endif
    endif

    call budget(rho1,ux1,uy1,uz1,phi1,vol1, diss1)
    call dep(phi1,dep2)
    call suspended(phi1,vol1,mp)
    call depositrate (dep2,dms)
    call front(phisum1,xf)
    if (numscalar.gt.0) then
       call front2d(phim3(:,:,1),xf2d)
    elseif (ilmn.and.ibirman_eos) then
       call front2d(rhom3(:,:),xf2d)
    endif

    if (nrank .eq. 0) then
       FS = 1+numscalar+numscalar+3+2 !Number of columns
       write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
       FS = FS*14+1  !Line width
       open(67,file='./data/statistics.csv',status='unknown',form='formatted',&
            access='direct',recl=FS)
       write(67,fileformat,rec=itime/iprocessing+1) t,& !1
            mp,&                                    !numscalar
            dms,&                                   !numscalar
            xf(1,:),&                               !3
            xf2d(1,:),&                             !2
            NL                                      !+1
       close(67)
    end if

    if (filenamedigits .eq. 0) then
      WRITE(num, ifilenameformat) itime
    else
      WRITE(num, ifilenameformat) itime/iprocessing
    endif

    call write_xdmf_header(filenum+1, num, './data/xdmf/xy_planes', nx, ny, 1)
    call write_xdmf_header(filenum+2, num, './data/xdmf/xz_planes', nx, 1, nz)


    call postprocessing_aux(ux1,ux2,ux3,'ux',num)
    call postprocessing_aux(uy1,uy2,uy3,'uy',num)
    call postprocessing_aux(uz1,uz2,uz3,'uz',num)
    call postprocessing_aux(diss1,diss2,diss3,'diss',num)

    if (iscalar.eq.1) then
      do is=1, numscalar
        call postprocessing_aux(phi1(:,:,:,is),phi2(:,:,:,is),phi3(:,:,:,is),'phi'//char(is+48),num)
      enddo
    endif

    call write_xdmf_footer(filenum+1)
    call write_xdmf_footer(filenum+2)

    return

  end subroutine postprocess_lockexch

  subroutine dep(phi1,dep2)

    USE decomp_2d_io
    USE MPI

    use var, only : phi2

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),numscalar) :: phi1
    real(mytype),intent(out),dimension(ystart(1):yend(1),ystart(3):yend(3),numscalar) :: dep2

    real(mytype),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),numscalar) :: tempdep2

    integer :: i, k, is
    character(len=30) :: filename

    dep2 = zero
    do is=1, numscalar
       if (uset(is) .eq. zero) cycle
       call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))

       tempdep2=zero

       do k=ystart(3),yend(3)
          do i=ystart(1),yend(1)
             tempdep2(i,1,k,is) = phi2(i,1,k,is)*uset(is)
             dep2(i,k,is) = tempdep2(i,1,k,is)
             !dep2(i,k,is) = dep2(i,k,is) +  phi2(i,1,k,is)*uset(is)
          end do
       end do

       ! write(filename,"('./out/dep',I1.1,I4.4)") is,itime/iprocessing
       ! call decomp_2d_write_plane(2,tempdep2(:,:,:,is),2,1,filename)
    enddo

  end subroutine dep

  subroutine suspended(phi1,vol1,mp1)

    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: vol1
    real(mytype),intent(out) :: mp1(1:numscalar)

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1
    real(mytype) :: mp(1:numscalar)
    integer :: is,code

    mp=zero; mp1=zero

    do is=1, numscalar
       temp1 = phi1(:,:,:,is)*vol1(:,:,:)
       mp(is)= sum(temp1)
    end do

    call MPI_REDUCE(mp,mp1,numscalar,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    return
  end subroutine suspended

  subroutine depositrate ( dep2, dms1)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(ystart(1):yend(1),ystart(3):yend(3),numscalar) :: dep2
    real(mytype),intent(out) :: dms1(numscalar)

    real(mytype) :: dms(numscalar)
    integer :: i,k,is,code

    dms=zero; dms1=zero
    do is=1, numscalar
       if (uset(is) .eq. zero) cycle
       do k=ystart(3),yend(3)
          do i=ystart(1),yend(1)
             dms(is)=dms(is)+dep2(i,k,is)*area2(i,k)
          end do
       end do
    enddo

    call MPI_REDUCE(dms,dms1,numscalar,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

  end subroutine depositrate

  subroutine front ( phisum1, xp )

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: phisum1
    real(mytype),intent(out) :: xp(1:2,1:3)

    real(mytype) :: xp1(1:2)
    integer :: i, j ,k, code

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

  subroutine front2d (phim3, xp)

    USE decomp_2d_io
    USE MPI

    real(mytype),intent(in),dimension(zstart(1):zend(1),zstart(2):zend(2)) :: phim3
    real(mytype),intent(out) :: xp(1:2,1:2)
    real(mytype) :: xp1(1:2),y
    integer :: i, j, code

    xp(2,:) = real(nrank,mytype)
    xp(1,:) = zero
    xp1=zero
    jloop: do j=zstart(2),zend(2)
       y = real( j - 1, mytype )*dy
       iloop: do i=zend(1), zstart(1), -1
          if (phim3(i,j).ge.0.01_mytype) then
             xp(1,1:2) = (/ real (i-1,mytype)*dx, real (j-1,mytype)*dy /)
             exit jloop
          end if
       end do iloop
    end do jloop

    call MPI_ALLREDUCE(xp(:,1),xp1,1,real2_type,MPI_MAXLOC,MPI_COMM_WORLD,code)
    call MPI_Bcast(xp(1,:), 2,real_type, int(xp1(2)), MPI_COMM_WORLD,code)

  end subroutine front2d

  subroutine set_fluid_properties_lockexch(rho1, mu1)

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

    if(.not.ilmn) return

    mu1(:,:,:) = rho1(:,:,:)

  endsubroutine set_fluid_properties_lockexch

end module lockexch
