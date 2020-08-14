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

module visu

  implicit none

  logical :: ixdmf = .True.
  integer :: filenamedigits = 0
  character(len=20) :: ifilenameformat = '(I9.9)'

  private
  public :: visu_init, write_snapshot, filenamedigits, ifilenameformat, &
    write_xdmf_header, write_field, write_xdmf_footer, ixdmf

contains
  !############################################################################
  !############################################################################
  subroutine visu_init()

    use param, only : ilmn, iscalar, ilast, ifirst, ioutput
    use variables, only : nx,ny,nz,numscalar,prec,beta
    use decomp_2d, only : nrank, mytype
    !
    implicit none
    !
    logical :: dir_exists
    integer :: noutput, nsnapout
    real(mytype) :: memout
    !
    !###################################################################
    ! Check if all folders exist
    !###################################################################
    if (nrank==0) then
      inquire(file="data", exist=dir_exists)
      if (.not.dir_exists) then
         call system("mkdir data 2> /dev/null")
      end if
      inquire(file="data/3d_snapshots", exist=dir_exists)
      if (.not.dir_exists) then
         call system("mkdir data/3d_snapshots 2> /dev/null")
      end if
      inquire(file="data/xy_planes", exist=dir_exists)
      if (.not.dir_exists) then
         call system("mkdir data/xy_planes 2> /dev/null")
      end if
      inquire(file="data/xz_planes", exist=dir_exists)
      if (.not.dir_exists) then
         call system("mkdir data/xz_planes 2> /dev/null")
      end if
      inquire(file="data/yz_planes", exist=dir_exists)
      if (.not.dir_exists) then
         call system("mkdir data/yz_planes 2> /dev/null")
      end if
      if (ixdmf) then
        inquire(file="data/xdmf", exist=dir_exists)
        if (.not.dir_exists) then
           call system("mkdir data/xdmf 2> /dev/null")
        end if
      endif
    end if
    !###################################################################
    ! Deprecisation warning
    !###################################################################
    if (nrank==0) then
      if (filenamedigits .eq. 1) then
        print *,'==========================================================='
        print *,'Deprecation warning:'
        print *,'    The classic enumeration system for output may be removed'
        print *,'    in a further release'
        print *,'==========================================================='
      endif
    end if
    !###################################################################
    ! Memory usage of visu module
    !###################################################################
    if (nrank==0) then
      noutput = (ilast - ifirst +1)/ioutput
      !
      nsnapout = 4
      if (ilmn)         nsnapout = nsnapout + 1
      if (iscalar.ne.0) nsnapout = nsnapout + numscalar
      !
      memout = nx*ny*nz*prec*nsnapout*noutput
      print *,'==========================================================='
      print *,'Visu module requires ',real(memout*1e-9,4),'GB'
      print *,'==========================================================='
    end if

  end subroutine visu_init
  !############################################################################
  !############################################################################
  subroutine write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime)

    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : mytype, xsize, ysize, zsize
    use decomp_2d, only : nrank
    use decomp_2d, only : fine_to_coarsev
    use decomp_2d_io, only : decomp_2d_write_one

    use param, only : ivisu, ioutput, nrhotime, ilmn, iscalar, iibm, istret
    use param, only : xlx, yly, zlz

    use variables, only : derx, dery, derz
    use variables, only : ffx, ffxp, fsx, fsxp, fwx, fwxp
    use variables, only : ffy, ffyp, fsy, fsyp, fwy, fwyp, ppy
    use variables, only : ffz, ffzp, fsz, fszp, fwz, fwzp
    use variables, only : sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    use variables, only : sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use variables, only : sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use variables, only : nx, ny, nz, beta, numscalar

    use var, only : one
    use var, only : uvisu
    use var, only : pp1, ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1, nxmsize
    use var, only : pp2, ta2, tb2, tc2, td2, te2, tf2, ppi2, di2, dip2, ph2, nymsize
    use var, only : ppi3, ta3, tb3, tc3, td3, te3, tf3, di3, dip3, ph3, nzmsize
    use var, only : npress
    use var, only : dt,t

    use tools, only : rescale_pressure

    use MPI

    implicit none

    character(len=80) :: filename
    !! inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize,npress), intent(in) :: pp3
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    integer, intent(in) :: itime

    character(len=32) :: fmt2,fmt3,fmt4
    integer :: code,ierror, is, io

    character(len=30) :: num

    real :: tstart, tend
    !###################################################################
    if (nrank.eq.0) then
      call cpu_time(tstart)
      print *,'Writing snapshots =>',itime/ioutput
    end if
    !
    if (filenamedigits .eq. 0) then
      ! New enumeration system, it works integrated with xcompact3d_toolbox
      WRITE(num, ifilenameformat) itime
    elseif (filenamedigits .eq. 1) then
      ! Classic enumeration system, may be removed in the near future
      WRITE(num, ifilenameformat) itime/ioutput
    else
       print *,'Invalid value for ilenamedigits, it should be 0 or 1'
       call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
    endif
    !
    io = 67
    call write_xdmf_header(io, num, './data/xdmf/3d_snapshots', nx, ny, nz)
    !###################################################################
    !! Write velocity
    !###################################################################
    call write_field(ux1, io, num, 'ux', nx, ny, nz, 1)
    call write_field(uy1, io, num, 'uy', nx, ny, nz, 1)
    call write_field(uz1, io, num, 'uz', nx, ny, nz, 1)
    !###################################################################
    !! Write pressure
    !###################################################################
    !WORK Z-PENCILS
    call interzpv(ppi3,pp3(:,:,:,1),dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
    (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
            (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)

    call rescale_pressure(ta1)

    call write_field(ta1, io, num, 'pp', nx, ny, nz, 1)
    !###################################################################
    !! LMN - write out density
    !###################################################################
      if (ilmn) call write_field(rho1(:,:,:,1), io, num, 'rho', nx, ny, nz, 1)
    !###################################################################
    !! Scalars
    !###################################################################
    if (iscalar.ne.0) then
      do is = 1, numscalar
        call write_field(phi1(:,:,:,is), io, num, 'phi'//CHAR(48+is), nx, ny, nz, 1)
      enddo
    endif

    call write_xdmf_footer(io)

    if (.not. ixdmf) then
      !###################################################################
      ! Write ini-file
      !###################################################################
      write(fmt2,'("(A,I16)")')
      write(fmt3,'("(A,F16.4)")')
      write(fmt4,'("(A,F16.12)")')
      if (nrank==0) then
         write(filename,"('./data/snap',I7.7,'.ini')") itime/ioutput
         !
         write(fmt2,'("(A,I16)")')
         write(fmt3,'("(A,F16.4)")')
         write(fmt4,'("(A,F16.12)")')

         open (844,file=filename,action='write',status='replace')
         write(844,'(A)')'[domain]'
         write(844,fmt2) 'nx=      ',nx
         write(844,fmt2) 'ny=      ',ny
         write(844,fmt2) 'nz=      ',nz
         write(844,fmt2) 'istret=  ',istret
         write(844,fmt4) 'beta=    ',beta
         write(844,fmt3) 'Lx=      ',xlx
         write(844,fmt3) 'Ly=      ',yly
         write(844,fmt3) 'Lz=      ',zlz
         write(844,'(A)')'[time]'
         write(844,fmt2) 'itime=   ',itime
         write(844,fmt3) 'dt=      ',dt
         write(844,fmt3) 't =      ',t
         close(844)
         !###################################################################
         call cpu_time(tend)
         write(*,'(" Time for writing snapshots (s): ",F12.8)') tend-tstart
         !###################################################################
    endif
  endif

  end subroutine write_snapshot
  !############################################################################
  subroutine write_xdmf_header(io, num, filename, nx, ny, nz)
    !############################################################################
    !!  SUBROUTINE: write_xdmf_header
    !!      AUTHOR: Felipe N. Schuch
    !! DESCRIPTION: Writes the header of the xdmf file, for data visualization
    !############################################################################
    use variables, only : nvisu, prec, yp
    use param, only : dx,dy,dz,istret
    USE decomp_2d, only : mytype, nrank

    implicit none

    integer, intent(in) :: io, nx, ny, nz
    character(len=*), intent(in) :: num, filename
    integer :: i,k

    real(mytype) :: xp(nx),zp(nz)

    if (.not. ixdmf) return

    if (nrank.eq.0) then
      OPEN(io,file=filename//'-'//trim(num)//'.xdmf')

      WRITE(io,'(A22)')'<?xml version="1.0" ?>'
      WRITE(io,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      WRITE(io,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
      WRITE(io,*)'<Domain>'
      if (istret.ne.0) then
        do i=1,nx
          xp(i) = real(i-1,mytype)*dx
        enddo
        do k=1,nz
          zp(k) = real(k-1,mytype)*dz
        enddo
        write(io,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
        write(io,*)'        Dimensions="',nz,ny,nx,'">'
        write(io,*)'    </Topology>'
        write(io,*)'    <Geometry name="geo" Type="VXVYVZ">'
        write(io,*)'        <DataItem Dimensions="',nx,'" NumberType="Float" Precision="4" Format="XML">'
        write(io,*)'        ',xp(:)
        write(io,*)'        </DataItem>'
        write(io,*)'        <DataItem Dimensions="',ny,'" NumberType="Float" Precision="4" Format="XML">'
        write(io,*)'        ',yp(:)
        write(io,*)'        </DataItem>'
        write(io,*)'        <DataItem Dimensions="',nz,'" NumberType="Float" Precision="4" Format="XML">'
        write(io,*)'        ',zp(:)
        write(io,*)'        </DataItem>'
        write(io,*)'    </Geometry>'
      else
        WRITE(io,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
        WRITE(io,*)'        Dimensions="',nz,ny,nx,'">'
        WRITE(io,*)'    </Topology>'
        WRITE(io,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
        WRITE(io,*)'        <!-- Origin -->'
        WRITE(io,*)'        <DataItem Format="XML" Dimensions="3">'
        WRITE(io,*)'        0.0 0.0 0.0'
        WRITE(io,*)'        </DataItem>'
        WRITE(io,*)'        <!-- DxDyDz -->'
        WRITE(io,*)'        <DataItem Format="XML" Dimensions="3">'
        WRITE(io,*)'        ',dz,dy,dx
        WRITE(io,*)'        </DataItem>'
        WRITE(io,*)'    </Geometry>'
      endif
      WRITE(io,*)'    <Grid Name="'//trim(num)//'" GridType="Uniform">'
      WRITE(io,*)'        <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      WRITE(io,*)'        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    endif
  end subroutine write_xdmf_header
  !######################################################################################
  subroutine write_xdmf_footer(io)
    !############################################################################
    !!  SUBROUTINE: write_xdmf_footer
    !!      AUTHOR: Felipe N. Schuch
    !! DESCRIPTION: Writes the footer of the xdmf file, for data visualization
    !############################################################################
    USE decomp_2d, only : nrank

    implicit none

    integer, intent(in) :: io

    if (ixdmf) then
      if (nrank.eq.0) then
        WRITE(io,'(/)')
        WRITE(io,*)'    </Grid>'
        WRITE(io,*)'</Domain>'
        WRITE(io,'(A7)')'</Xdmf>'
        close(io)
      endif
    endif

  end subroutine write_xdmf_footer
  !######################################################################################
  subroutine write_field(f1, io, num, filename, nx, ny, nz, clean_ibm)
    !############################################################################
    !!  SUBROUTINE: write_field
    !!      AUTHOR: Felipe N. Schuch
    !! DESCRIPTION: Writes the body of the xdmf file, in addition to write the
    !!              binary fields to disc, in a DRY (Don't repeat yourself) way
    !############################################################################
    use variables, only : nvisu, prec
    use var, only : ta1, ep1
    use var, only : zero, one
    use var, only : uvisu
    use param, only : iibm
    use decomp_2d, only : mytype, xsize, nrank
    use decomp_2d, only : fine_to_coarsev
    use decomp_2d_io, only : decomp_2d_write_one

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: f1
    integer, intent(in) :: io, nx, ny, nz, clean_ibm
    character(len=*), intent(in) :: num, filename
    !
    if (ixdmf) then
      !
      if (nrank.eq.0) then
        write(io,*)'        <Attribute Name="'//filename//'" Center="Node">'
        write(io,*)'           <DataItem Format="Binary"'
        write(io,*)'            DataType="Float" Precision="'//CHAR(48+prec)//'" Endian="little" Seek="0"'
        write(io,*)'            Dimensions="',nz,ny,nx,'">'
        write(io,*)'              ../3d_snapshots/'//filename//'-'//trim(num)//'.bin'
        write(io,*)'           </DataItem>'
        write(io,*)'        </Attribute>'
      endif
      !
    endif
    uvisu = zero
    ta1(:,:,:) = f1(:,:,:)
    if (clean_ibm.ne.0.and.iibm.ne.0) ta1(:,:,:) = (one - ep1(:,:,:)) * ta1(:,:,:)
    !
    call fine_to_coarseV(1,ta1,uvisu)
    !
    call decomp_2d_write_one(1,uvisu,'./data/3d_snapshots/'//filename//'-'//trim(num)//'.bin')

  end subroutine write_field
  !############################################################################
end module visu
