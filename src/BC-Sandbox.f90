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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-Sandbox.f90
!!!      AUTHOR: Felipe N. Schuch <felipe.schuch@edu.pucrs.br>
!!! DESCRIPTION: This work aims to break many of the barriers to entry in a
!!!              Navier-Stokes solver by coupling it to a Jupyter sandbox
!!!              environment. For students in computational fluid dynamics, it
!!!              provides direct hands-on experience and a safe place for
!!!              practising and learning, while for advanced users and code
!!!              developers, it works as a rapid prototyping tool.
!!!
!!!             Try it together with the Xcompact3d_toolbox:
!!!             https://github.com/fschuch/Xcompact3d_toolbox
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module sandbox

  use decomp_2d, only : mytype, real_type, real2_type
  use decomp_2d, only : xsize, ysize, zsize
  use decomp_2d, only : xstart, ystart, zstart
  use decomp_2d, only : xend, yend, zend
  use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, &
      transpose_z_to_y, transpose_y_to_x

  use variables, only : numscalar

  use var, only : xnu, ri, uset, sc, Fr
  use var, only : gravx, gravy, gravz
  use var, only : nrank

  use var, only : zero, half, one, two, five, twelve, thirteen
  use var, only : dx, dy, dz, nx, ny, nz

  use param, only : itime, ioutput, iprocessing
  use param, only : t

  use visu, only : ixdmf

  implicit none

  character(len=120)  :: filename
  integer, parameter :: filenum = 67

#ifdef DOUBLE_PREC
#if defined(__GFORTRAN__)
  integer, parameter :: filerecl = 8
#else
  integer, parameter :: filerecl = 2
#endif
#else
#if defined(__GFORTRAN__)
  integer, parameter :: filerecl = 4
#else
  integer, parameter :: filerecl = 1
#endif
#endif

  logical, save :: init = .FALSE.
  real(mytype), save, allocatable :: vol1_frc(:,:,:)
  real(mytype), save, allocatable :: bxx1_sb(:,:), bxy1_sb(:,:), bxz1_sb(:,:), noise_mod_x1(:,:)
  real(mytype), save, allocatable :: bxphi1(:,:,:), byphi1(:,:,:), byphin(:,:,:)
  real(mytype), save, allocatable :: vol1(:,:,:)
  integer, save :: sz_a_i = 0
  real(mytype), save :: sz_a_length = 0.0
  character(len=1),parameter :: NL=char(10) !new line character

  private
  public :: init_sandbox, boundary_conditions_sandbox, postprocess_sandbox, &
  geomcomplex_sandbox, deposit, postprocessing_aux, budget

contains

  subroutine geomcomplex_sandbox(epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, yp, remp)

    use decomp_2d, only : mytype, xstart, xend
    use decomp_2d_io, only : decomp_2d_read_one
    use param, only : one, two
    use variables, only : nx, nz
    use complex_geometry, only : nxraf,nyraf,nzraf
    use ibm
    use MPI

    implicit none

    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp
    real(mytype)               :: dx,dz
    real(mytype)               :: remp
    integer                    :: code, ierror
    !
    if (nxi.eq.1.and.nxf.eq.nx.and.&
        nyi.eq.xstart(2).and.nyf.eq.xend(2).and.&
        nzi.eq.xstart(3).and.nzf.eq.xend(3)) then
        !
        call decomp_2d_read_one(1,epsi,'./data/geometry/epsilon.bin')
        !
    else
      print *,'Invalid parameters at geomcomplex_sandbox'
      call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
    endif
    !
    return
  end subroutine geomcomplex_sandbox

  subroutine boundary_conditions_sandbox (ux, uy, uz, phi1)

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    if (nclx1 .eq. 2) call inflow (phi1)
    if (nclxn .eq. 2) call outflow (ux,uy,uz,phi1)
    if (nclx) call flow_rate_control(ux)
    if (sz_a_i .gt. 1) call flow_rate_control_SZA(ux, uy, uz)
    if (iscalar .eq. 1) call deposit(phi1)

    return
  end subroutine boundary_conditions_sandbox
  !********************************************************************
  subroutine deposit (phi1)
    !================================================================================
    !
    !  SUBROUTINE: deposit
    ! DESCRIPTION:
    !      AUTHOR: Felipe N. Schuch <felipe.schuch@edu.pucrs.br>
    !
    !================================================================================
    USE param
    USE variables
    USE decomp_2d
    USE MPI
    USE var, only : phi2, phi3, ta2, di2

    implicit none

    integer  :: i,j,k,is
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype) :: temp

    !Robin on top BC, following Schuch (2020)
    if (nclySn .eq. 2 .and. xend(2) .eq. ny) then
      j = xsize(2)
      do is=1, numscalar
        if (uset(is) .eq. zero .and. itype.eq.itype_sandbox) then !For no-flux use nclySn = 1
          do k=1, xsize(3)
            do i=1, xsize(1)
              phi1(i,j,k,is) = byphin(i,k,is)
            enddo
          enddo
        else
          temp = one / (six*re*sc(is)*dy*uset(is) + eleven)
          do k=1, xsize(3)
            do i=1, xsize(1)
              phi1(i,j,k,is)= temp*(two*phi1(i,j-3,k,is)&
              -nine*phi1(i,j-2,k,is)&
              +nine*two*phi1(i,j-1,k,is)&
              )
            enddo
          enddo
        endif
      enddo
    endif

    !Deposit on bottom BC
    if (nclyS1 .eq. 2) then
      do is=1, numscalar
        if (uset(is) .eq. zero) cycle
        call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))
        call deryS (ta2, phi2, di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, iibmS)
          do k=1, ysize(3)
            do i=1, ysize(1)
              if (ta2(i,1,k) .lt. zero) then
                phi2(i,1,k,is) = (-four*dy*ta2(i,2,k)+four*phi2(i,2,k,is)+phi2(i,3,k,is))/five
              endif
            enddo
          enddo
        call transpose_y_to_x(phi2(:,:,:,is),phi1(:,:,:,is))
      enddo
    endif

    return
  end subroutine deposit
  !********************************************************************
  subroutine inflow (phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    integer  :: j,k,is
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)

    bxo = (two*bxo - one)*inflow_noise*noise_mod_x1
    byo = (two*byo - one)*inflow_noise*noise_mod_x1
    bzo = (two*bzo - one)*inflow_noise*noise_mod_x1

    bxx1(:,:) = bxx1_sb(:,:) + bxo(:,:)
    bxy1(:,:) = bxy1_sb(:,:) + byo(:,:)
    bxz1(:,:) = bxz1_sb(:,:) + bzo(:,:)

    if (iscalar.eq.1 .and. nclxS1.eq.2) then
      phi(1,:,:,:) = bxphi1(:,:,:)
    endif

    return
  end subroutine inflow
  !********************************************************************
  subroutine outflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax,uxmin1,uxmax1

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609.
    uxmin=1609.
    do k=1,xsize(3)
      do j=1,xsize(2)
        if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
        if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
      enddo
    enddo

    call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    if (u1.eq.zero) then
      cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
    elseif (u1.eq.one) then
      cx=uxmax1*gdt(itr)*udx
    elseif (u1.eq.two) then
      cx=u2*gdt(itr)*udx    !works better
    else
      stop
    endif

    do k=1,xsize(3)
      do j=1,xsize(2)
        bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
        bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
        bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
      enddo
    enddo

    if (iscalar==1) then
      if (u2.eq.zero) then
        cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
      elseif (u2.eq.one) then
        cx=uxmax1*gdt(itr)*udx
      elseif (u2.eq.two) then
        cx=u2*gdt(itr)*udx    !works better
      else
        stop
      endif

      do k=1,xsize(3)
        do j=1,xsize(2)
          phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
        enddo
      enddo
    endif

    if (nrank==0) write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)

    return
  end subroutine outflow
  !********************************************************************
  subroutine flow_rate_control (ux)
    !================================================================================
    !
    !  SUBROUTINE: flow_rate_control
    ! DESCRIPTION:
    !      AUTHOR: Felipe N. Schuch <felipe.schuch@edu.pucrs.br>
    !
    !================================================================================
    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer  :: i,j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)):: ux, tmp
    real(mytype) :: int, int1

    int = zero; int1 = zero

    tmp(:,:,:) = vol1_frc(:,:,:) * ux(:,:,:)

    int = sum(tmp)

    call MPI_ALLREDUCE(int,int1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    if (nrank==0) write(*,*) "Integration at frc : ",int1

    ux = ux / int1

    return
  end subroutine flow_rate_control
  !********************************************************************
  subroutine flow_rate_control_SZA (ux, uy, uz)
    !================================================================================
    !
    !  SUBROUTINE: flow_rate_control_SZA
    ! DESCRIPTION: flow rate control at Sponge Zone A, for Recycling Turbulent
    !              inflow condition, following Schuch (2020)
    !      AUTHOR: Felipe N. Schuch <felipe.schuch@edu.pucrs.br>
    !
    !================================================================================
    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer  :: i,j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)):: ux, uy, uz, tmp
    real(mytype) :: int, int1

    int = zero; int1 = zero

    tmp(:,:,:) = vol1_frc(:,:,:) * ux(:,:,:)

    int = sum(tmp)

    call MPI_ALLREDUCE(int,int1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    if (nrank==0) write(*,*) "Integration at frc SZA : ",int1

    !Flow rate control
    do i=1, sz_a_i
      ux(i,:,:) = ux(i,:,:) / int1
    enddo

    !Recycling four points just to avoid the diferenciation schemes near boundary
    do i=1, 4
      ux(i,:,:) = ux(sz_a_i+i-1,:,:)
      uy(i,:,:) = uy(sz_a_i+i-1,:,:)
      uz(i,:,:) = uz(sz_a_i+i-1,:,:)
    enddo

    !Get Boundary Condition

    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)

    bxo = (two*bxo - one)*inflow_noise*noise_mod_x1
    byo = (two*byo - one)*inflow_noise*noise_mod_x1
    bzo = (two*bzo - one)*inflow_noise*noise_mod_x1

    bxx1(:,:) = ux(1,:,:) + bxo
    bxy1(:,:) = uy(1,:,:) + byo
    bxz1(:,:) = uz(1,:,:) + bzo

    return
  end subroutine flow_rate_control_SZA
  !********************************************************************
  subroutine init_sandbox (ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    USE var, only : mu1

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    integer :: i, j, k, is, pos

    ! NAMELIST /CASE/ sz_a_length
    !
    ! !! Read parameters
    ! open(10, file=InputFN)
    ! read(10, nml=CASE); rewind(10) !! Read case-specific variables
    ! close(10)

    !! Sponge Zone A for recycling inflow condition
    if (sz_a_length .gt. 0.0) sz_a_i = int(sz_a_length / dx) + 1

    !Read phi
    if (iscalar .ne. 0) then
      do is = 1, numscalar
        if (nrank.eq.0) write(*,*) 'reading : ', './data/phi'//char(is+48)//'.bin'
        call decomp_2d_read_one(1, phi1(:,:,:,is), './data/phi'//char(is+48)//'.bin')
      enddo
    endif

    !Read velocity field
    if (nrank.eq.0) write(*,*) 'reading : ', './data/ux.bin'
    call decomp_2d_read_one(1,ux1,'./data/ux.bin')
    if (nrank.eq.0) write(*,*) 'reading : ', './data/uy.bin'
    call decomp_2d_read_one(1,uy1,'./data/uy.bin')
    if (nrank.eq.0) write(*,*) 'reading : ', './data/uz.bin'
    call decomp_2d_read_one(1,uz1,'./data/uz.bin')

    !Read integration operator for flow_rate_control
    if (nclx .or. sz_a_i .gt. 1) then
      call alloc_x(vol1_frc)
      vol1_frc = zero
      filename = './data/vol_frc.bin'
      if (nrank.eq.0) write(*,*) 'reading : ', filename
      call decomp_2d_read_one(1,vol1_frc,filename)
    endif

    !Read inflow profile
    if (nclx1 .eq. 2) then
      !
      allocate(bxx1_sb(xsize(2),xsize(3)))
      allocate(bxy1_sb(xsize(2),xsize(3)))
      allocate(bxz1_sb(xsize(2),xsize(3)))
      allocate(noise_mod_x1(xsize(2),xsize(3)))
      !
      if (nrank.eq.0) write(*,*) 'reading : ', './data/bxx1.bin'
      if (nrank.eq.0) write(*,*) 'reading : ', './data/bxy1.bin'
      if (nrank.eq.0) write(*,*) 'reading : ', './data/bxz1.bin'
      if (nrank.eq.0) write(*,*) 'reading : ', './data/noise_mod_x1.bin'
      !
      OPEN(filenum+1,FILE='./data/bxx1.bin',FORM='UNFORMATTED',ACCESS='DIRECT', RECL=filerecl, STATUS='OLD')
      OPEN(filenum+2,FILE='./data/bxy1.bin',FORM='UNFORMATTED',ACCESS='DIRECT', RECL=filerecl, STATUS='OLD')
      OPEN(filenum+3,FILE='./data/bxz1.bin',FORM='UNFORMATTED',ACCESS='DIRECT', RECL=filerecl, STATUS='OLD')
      OPEN(filenum+4,FILE='./data/noise_mod_x1.bin',FORM='UNFORMATTED',ACCESS='DIRECT', RECL=filerecl, STATUS='OLD')
      !
      do k=1,xsize(3)
        do j=1,xsize(2)
          pos = (k - 1 + xstart(3) - 1) * ny + j + xstart(2) - 1
          read(filenum+1,rec=pos) bxx1_sb(j,k)
          read(filenum+2,rec=pos) bxy1_sb(j,k)
          read(filenum+3,rec=pos) bxz1_sb(j,k)
          read(filenum+4,rec=pos) noise_mod_x1(j,k)
        enddo
      enddo
      close(filenum+1)
      close(filenum+2)
      close(filenum+3)
      close(filenum+4)
    endif

    !Read phi inflow profile
    if (iscalar .ne. 0 .and. nclxS1 .eq. 2) then
      allocate(bxphi1(xsize(2),xsize(3),numscalar))
      bxphi1 = zero
      do is=1, numscalar
        !
        filename = './data/bxphi1'//char(is+48)//'.bin'
        if (nrank.eq.0) write(*,*) 'reading : ', filename
        OPEN(filenum,FILE=filename,FORM='UNFORMATTED',ACCESS='DIRECT', RECL=filerecl, STATUS='OLD')
        !
        do k=1,xsize(3)
          do j=1,xsize(2)
            pos = (k - 1 + xstart(3) - 1) * ny + j + xstart(2) - 1
            read(filenum,rec=pos) bxphi1(j,k,is)
          enddo
        enddo
        close(filenum)
        !
      enddo
    endif
    !Read phi bottom BC
    if (iscalar .ne. 0 .and. nclyS1 .eq. 2) then
      allocate(byphi1(xsize(1),xsize(3),numscalar))
      byphi1 = zero
      do is=1, numscalar
        !
        if (uset(is) .ne. zero) cycle !in this case we use deposition BC
        !
        filename = './data/byphi1'//char(is+48)//'.bin'
        if (nrank.eq.0) write(*,*) 'reading : ', filename
        OPEN(filenum,FILE=filename,FORM='UNFORMATTED',ACCESS='DIRECT', RECL=filerecl, STATUS='OLD')
        !
        do k=1,xsize(3)
          do i=1,xsize(1)
            pos = (k - 1 + xstart(3) - 1) * nx + i + xstart(1) - 1
            read(filenum,rec=pos) byphi1(i,k,is)
          enddo
        enddo
        close(filenum)
        !
      enddo
    endif
    !Read phi top BC
    if (iscalar .ne. 0 .and. nclySn .eq. 2) then
      allocate(byphin(xsize(2),xsize(3),numscalar))
      byphin = zero
      do is=1, numscalar
        !
        if (uset(is) .ne. zero) cycle !in this case we use deposition BC
        !
        filename = './data/byphin'//char(is+48)//'.bin'
        if (nrank.eq.0) write(*,*) 'reading : ', filename
        OPEN(filenum,FILE=filename,FORM='UNFORMATTED',ACCESS='DIRECT', RECL=filerecl, STATUS='OLD')
        !
        do k=1,xsize(3)
          do i=1,xsize(1)
            pos = (k - 1 + xstart(3) - 1) * nx + i + xstart(1) - 1
            read(filenum,rec=pos) byphin(i,k,is)
          enddo
        enddo
        close(filenum)
        !
      enddo
    endif

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

    return
  end subroutine init_sandbox

  subroutine postprocess_sandbox(ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

    use param, only : iscalar
    use visu, only : filenamedigits, ifilenameformat, &
        write_xdmf_header, write_field, write_xdmf_footer
    use var, only : rho1
    use var, only : ux2,uy2,uz2,phi2
    use var, only : ux3,uy3,uz3,phi3

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: diss1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: diss2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: diss3

    integer :: is
    character(len=30) :: num

    if (mod(itime,iprocessing).ne.0) return

    if (filenamedigits .eq. 0) then
      WRITE(num, ifilenameformat) itime
    else
      WRITE(num, ifilenameformat) itime/iprocessing
    endif

    call write_xdmf_header(filenum+1, num, './data/xdmf/xy_planes', nx, ny, 1)
    call write_xdmf_header(filenum+2, num, './data/xdmf/xz_planes', nx, 1, nz)

    ! <ToDo> Control this call with an if
    call budget(rho1,ux1,uy1,uz1,phi1,vol1, diss1)
    call postprocessing_aux(diss1,diss2,diss3,'diss',num)

    call postprocessing_aux(ux1,ux2,ux3,'ux',num)
    call postprocessing_aux(uy1,uy2,uy3,'uy',num)
    call postprocessing_aux(uz1,uz2,uz3,'uz',num)

    if (iscalar.eq.1) then
      do is=1, numscalar
        call postprocessing_aux(phi1(:,:,:,is),phi2(:,:,:,is),phi3(:,:,:,is),'phi'//char(is+48),num)
      enddo
    endif

    call write_xdmf_footer(filenum+1)
    call write_xdmf_footer(filenum+2)


    return
  end subroutine postprocess_sandbox

  subroutine postprocessing_aux(u1,u2,u3,name,num) !By Felipe Schuch

    USE decomp_2d
    USE decomp_2d_io
    USE param, only : zero, one
    USE var, only : ep1
    USE variables, only: nx, ny, nz, prec
    use tools, only : mean_plane_y, mean_plane_z

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: u1
    real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u2
    real(mytype),intent(out),dimension(zsize(1),zsize(2),zsize(3)) :: u3
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: num
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tmp1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tmp2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: tmp3
    !
    tmp1 = u1*(one-ep1)
    tmp2 = zero
    tmp3 = zero
    call transpose_x_to_y (tmp1,u2)
    call transpose_y_to_z (u2,u3)
    !depth averaged
    call mean_plane_y(u2,ysize(1),ysize(2),ysize(3),tmp2(:,1,:))
    filename = './data/xz_planes/'//trim(name)//'-'//trim(num)//'.bin'
    call decomp_2d_write_plane(2,tmp2,2,1,filename)
    !
    if (nrank.eq.0 .and. ixdmf) then
      write(filenum+2,*)'        <Attribute Name="'//trim(name)//'" Center="Node">'
      write(filenum+2,*)'           <DataItem Format="Binary"'
      write(filenum+2,*)'            DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
      write(filenum+2,*)'            Dimensions="',nz,1,nx,'">'
      write(filenum+2,*)'              ../xz_planes/'//trim(name)//'-'//trim(num)//'.bin'
      write(filenum+2,*)'           </DataItem>'
      write(filenum+2,*)'        </Attribute>'
    endif
    !spanwise averaged
    call mean_plane_z(u3,zsize(1),zsize(2),zsize(3),tmp3(:,:,1))
    filename = './data/xy_planes/'//trim(name)//'-'//trim(num)//'.bin'
    call decomp_2d_write_plane(3,tmp3,3,1,filename)
    !
    if (nrank.eq.0 .and. ixdmf) then
      write(filenum+1,*)'        <Attribute Name="'//trim(name)//'" Center="Node">'
      write(filenum+1,*)'           <DataItem Format="Binary"'
      write(filenum+1,*)'            DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
      write(filenum+1,*)'            Dimensions="',1,ny,nx,'">'
      write(filenum+1,*)'              ../xy_planes/'//trim(name)//'-'//trim(num)//'.bin'
      write(filenum+1,*)'           </DataItem>'
      write(filenum+1,*)'        </Attribute>'
    endif
    !center plane
    filename = './data/xy_planes/'//trim(name)//'c-'//trim(num)//'.bin'
    call decomp_2d_write_plane(3,u3,3,nz/2,filename)
    !
    if (nrank.eq.0 .and. ixdmf) then
      write(filenum+1,*)'        <Attribute Name="'//trim(name)//'c'//'" Center="Node">'
      write(filenum+1,*)'           <DataItem Format="Binary"'
      write(filenum+1,*)'            DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
      write(filenum+1,*)'            Dimensions="',1,ny,nx,'">'
      write(filenum+1,*)'              ../xy_planes/'//trim(name)//'c-'//trim(num)//'.bin'
      write(filenum+1,*)'           </DataItem>'
      write(filenum+1,*)'        </Attribute>'
    endif
    !bot plane
    filename = './data/xz_planes/'//trim(name)//'b-'//trim(num)//'.bin'
    call decomp_2d_write_plane(2,u2,2,1,filename)
    !
    if (nrank.eq.0 .and. ixdmf) then
      write(filenum+2,*)'        <Attribute Name="'//trim(name)//'b'//'" Center="Node">'
      write(filenum+2,*)'           <DataItem Format="Binary"'
      write(filenum+2,*)'            DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
      write(filenum+2,*)'            Dimensions="',nz,1,nx,'">'
      write(filenum+2,*)'              ../xz_planes/'//trim(name)//'b-'//trim(num)//'.bin'
      write(filenum+2,*)'           </DataItem>'
      write(filenum+2,*)'        </Attribute>'
    endif
    !top plane
    filename = './data/xz_planes/'//trim(name)//'t-'//trim(num)//'.bin'
    call decomp_2d_write_plane(2,u2,2,ny,filename)
    !
    if (nrank.eq.0 .and. ixdmf) then
      write(filenum+2,*)'        <Attribute Name="'//trim(name)//'t'//'" Center="Node">'
      write(filenum+2,*)'           <DataItem Format="Binary"'
      write(filenum+2,*)'            DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
      write(filenum+2,*)'            Dimensions="',nz,1,nx,'">'
      write(filenum+2,*)'              ../xz_planes/'//trim(name)//'t-'//trim(num)//'.bin'
      write(filenum+2,*)'           </DataItem>'
      write(filenum+2,*)'        </Attribute>'
    endif

  end subroutine postprocessing_aux

  subroutine budget(rho1,ux1,uy1,uz1,phi1,vol1, diss1)

    USE decomp_2d
    USE decomp_2d_io
    USE MPI

    use variables, only : derx, dery, derys, derz
    use variables, only : derxx, deryy, derzz
    use variables, only : derxxs, deryys, derzzs

    use param, only : iibm, iibmS
    use param, only : ilmn

    use var, only : ffx, ffxp, fsx, fsxp, fwx, fwxp, sfxps, ssxps, swxps, sx, &
         sfxp, ssxp, swxp
    use var, only : ffy, ffyp, ffys, fsy, fsyp, fsys, fwy, fwyp, fwys, ppy, sfyps, ssyps, swyps, sy, &
         sfyp, ssyp, swyp
    use var, only : ffz, ffzp, fsz, fszp, fwz, fwzp, sfzps, sszps, swzps, sz, &
         sfzp, sszp, swzp

    use var, only : prandtl
    use var, only : nrhotime
    use var, only : mu1
    use var, only : rho2, ux2, uy2, uz2, phi2
    use var, only : rho3, ux3, uy3, uz3, phi3

    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    use var, only : ta2,tb2,tc2,td2,te2,tf2,di2
    use var, only : ta3,tb3,tc3,di3

    use tools, only : mean_plane_z

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,vol1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: diss1, dphixx1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: dphiy2, dphixx2, dphiyy2, dphizz2, ddphi2, vol2, temp2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: dphizz3, temp3

    real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A

    real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu

    real(8) :: ek,ek1,dek,dek1,ep,ep1,dep,dep1,xvol
    integer :: ijk,i,j,k,l,m,is,code
    character(len=30) :: filename

    real(mytype) :: y

    ek=zero;ek1=zero;dek=zero;dek1=zero;ep=zero;ep1=zero;dep=zero;dep1=zero;diss1=zero

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,iibm)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,iibm)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,iibm)
    !y-derivatives
    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,iibm)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,iibm)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,iibm)
    !!z-derivatives
    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,iibm)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,iibm)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,iibm)
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

    if (ilmn) then

      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               do m=1,3
                  do l=1,3
                     diss1(i,j,k) = diss1(i,j,k) &
                          + (two * xnu * mu1(i, j, k)) &
                          * (half * (A(l,m,i,j,k) + A(m,l,i,j,k)))**two
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ijk=1,xsize(1)*xsize(2)*xsize(3)
         xvol=real(vol1(ijk,1,1),8)
         ek = ek + half * xvol * rho1(ijk,1,1,1) * (ux1(ijk,1,1)**two+uy1(ijk,1,1)**two+uz1(ijk,1,1)**two)
         dek = dek + xvol * diss1(ijk,1,1)
      enddo

    else

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

    endif

    call transpose_x_to_y(vol1,vol2)

    ! if (ivirt==2) then
    !    ilag=0
    ! endif
    do is=1, numscalar
       if (ri(is) .eq. zero) cycle
       call derxxS (dphixx1,phi1(:,:,:,is),di1,sx,sfxpS,ssxpS,swxpS,xsize(1),xsize(2),xsize(3),1,iibmS)

       call transpose_x_to_y(dphixx1,dphixx2)

       call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))

       call deryS (dphiy2,phi2(:,:,:,is),di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),1,iibmS)

       call deryyS (dphiyy2,phi2(:,:,:,is),di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1,iibmS)

       call transpose_y_to_z(phi2(:,:,:,is),phi3(:,:,:,is))

       call derzzS (dphizz3,phi3(:,:,:,is),di3,sz,sfzpS,sszpS,swzpS,zsize(1),zsize(2),zsize(3),1,iibmS)

       call transpose_z_to_y(dphizz3,dphizz2)

       ddphi2(:,:,:)=dphixx2(:,:,:)+dphiyy2(:,:,:)+dphizz2(:,:,:)

       do k=1,ysize(3)
          do j=1,ysize(2)
             y = (j + ystart(2) - 2) * dy
             do i=1,ysize(1)
                xvol=real(vol2(i,j,k),8)
                ep = ep - xvol * ri(is) * phi2(i,j,k,is) * (gravy * y)
                dep = dep &
                     - xvol * ri(is) * (ddphi2(i,j,k)*xnu/sc(is) &
                     - uset(is) * (gravy * dphiy2(i,j,k))) &
                     * (gravy * y)
             enddo
          enddo
       enddo
    enddo

    if (ilmn.and.((Fr**2).gt.zero)) then
       call derxx(ta1, rho1(:,:,:,1), di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1, iibm)
       call transpose_x_to_y(ta1, ta2)
       call transpose_x_to_y(rho1(:,:,:,1), rho2)

       call deryy(tb2, rho2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1, iibm)
       call transpose_y_to_z(rho2, rho3)

       call derzz(ta3, rho3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1, iibmS)
       call transpose_z_to_y(ta3, tc2)

       do k = 1, ysize(3)
          do j = 1, ysize(2)
             y = (j + ystart(2) - 2) * dy
             do i = 1, ysize(1)
                xvol = real(vol2(i, j, k), 8)
                ep = ep - xvol * (one / Fr**2) * rho2(i, j, k) * (gravy * y)
                dep = dep - xvol * ((xnu / prandtl / (Fr**2)) &
                     * (ta2(i, j, k) + tb2(i, j, k) + tc2(i, j, k))) &
                     * (gravy * y)
             enddo
          enddo
       enddo
    endif
    ! if (ivirt==2) then
    !    ilag=1
    ! endif

    call MPI_REDUCE(ek,ek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(dek,dek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(ep,ep1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(dep,dep1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)

    if (nrank .eq. 0) then
       open(67,file='./data/budget.csv',status='unknown',form='formatted',&
            access='direct',recl=71) !71=5*14+1
       write(67,"(5E14.6,A)",rec=itime/iprocessing+1) t,ek1,dek1,ep1,dep1,NL
       close(67)
    end if

  end subroutine budget

end module sandbox
