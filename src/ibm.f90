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

module ibm

  public

contains
  !############################################################################
  subroutine corgp_IBM (ux,uy,uz,px,py,pz,nlock)
    USE param
    USE decomp_2d
    USE variables
    implicit none
    integer :: i,j,k,nlock
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz
    if (nlock.eq.1) then
       if (nz.gt.1) then
          do k = 1, xsize(3)
             do j = 1, xsize(2)
                do i = 1, xsize(1)
                   ux(i,j,k)=-px(i,j,k)+ux(i,j,k)
                   uy(i,j,k)=-py(i,j,k)+uy(i,j,k)
                   uz(i,j,k)=-pz(i,j,k)+uz(i,j,k)
                enddo
             enddo
          enddo
       else
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ux(i,j,k)=-px(i,j,k)+ux(i,j,k)
                uy(i,j,k)=-py(i,j,k)+uy(i,j,k)
             enddo
          enddo
       endif
    endif
    if (nlock.eq.2) then
       if (nz.gt.1) then
          do k = 1, xsize(3)
             do j = 1, xsize(2)
                do i = 1, xsize(1)
                   ux(i,j,k)=px(i,j,k)+ux(i,j,k)
                   uy(i,j,k)=py(i,j,k)+uy(i,j,k)
                   uz(i,j,k)=pz(i,j,k)+uz(i,j,k)
                enddo
             enddo
          enddo
       else
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ux(i,j,k)=px(i,j,k)+ux(i,j,k)
                uy(i,j,k)=py(i,j,k)+uy(i,j,k)
             enddo
          enddo
       endif
    endif

    return
  end subroutine corgp_IBM
  !############################################################################
  !############################################################################
  subroutine body(ux1,uy1,uz1,ep1)
    USE param, only : zero, one, dx, dz
    USE decomp_2d, only : xstart, xend, xsize, mytype
    !USE decomp_2d_io
    USE variables, only : ny
    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    integer :: i,j,k

#ifdef DEBG
    if (nrank .eq. 0) print *,'# body start'
#endif

    do k = 1, xsize(3)
      do j = 1, xsize(2)
         do i = 1, xsize(1)
            ux1(i,j,k)=(one-ep1(i,j,k))*ux1(i,j,k)
            uy1(i,j,k)=(one-ep1(i,j,k))*uy1(i,j,k)
            uz1(i,j,k)=(one-ep1(i,j,k))*uz1(i,j,k)
         enddo
      enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# body done'
#endif

    return
  end subroutine body
  !############################################################################
  !############################################################################
  subroutine lagpolx(u)
    !
    USE param
    USE complex_geometry
    USE decomp_2d
    USE variables
    !
    implicit none
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u
    integer                                            :: i,j,k
    real(mytype)                                       :: x,y,z
    integer                                            :: ix              !=   position du point "zappé"
    integer                                            :: ipif,ipol,nxpif
    integer                                            :: ipoli,ipolf     !=   positions Initiales et Finales du POLynôme considéré
    real(mytype)                                       :: xpol,ypol,dypol   !|variables concernant les polynômes
    real(mytype),dimension(10)                         :: xa,ya           !|de   Lagrange. A mettre impérativement en
    integer                                            :: ia,na             !|double précision
    !
    do k=1,xsize(3)
       do j=1,xsize(2)
          if(nobjx(j,k).ne.0)then
             ia=0
             do i=1,nobjx(j,k)          !boucle sur le nombre d'objets par (j,k)
                !1ère frontière
                nxpif=npif
                ia=ia+1
                xa(ia)=xi(i,j,k)
                ya(ia)=0.
                if(xi(i,j,k).gt.0.)then!objet immergé
                   ix=xi(i,j,k)/dx+1
                   ipoli=ix+1
                   if(nxipif(i,j,k).lt.npif)nxpif=nxipif(i,j,k)
                   do ipif=1,nxpif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=(ix-1)*dx-ipif*dx
                         ya(ia)=u(ix-ipif,j,k)
                      else             !no zapping
                         xa(ia)=(ix-1)*dx-(ipif-1)*dx
                         ya(ia)=u(ix-ipif+1,j,k)
                      endif
                   enddo
                else                   !objet semi-immergé
                   ipoli=1
                endif
                !2ème frontière
                nxpif=npif
                ia=ia+1
                xa(ia)=xf(i,j,k)
                ya(ia)=0.
                if(xf(i,j,k).lt.xlx)then!objet immergé
                   ix=(xf(i,j,k)+dx)/dx+1
                   ipolf=ix-1
                   if(nxfpif(i,j,k).lt.npif)nxpif=nxfpif(i,j,k)
                   do ipif=1,nxpif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=(ix-1)*dx+ipif*dx
                         ya(ia)=u(ix+ipif,j,k)
                      else             !no zapping
                         xa(ia)=(ix-1)*dx+(ipif-1)*dx
                         ya(ia)=u(ix+ipif-1,j,k)
                      endif
                   enddo
                else                   !objet semi-immergé
                   ipolf=nx
                endif
                !calcul du polynôme
                na=ia
                do ipol=ipoli,ipolf
                   xpol=dx*(ipol-1)
                   call polint(xa,ya,na,xpol,ypol,dypol)
                   u(ipol,j,k)=ypol
                enddo
                ia=0
             enddo
          endif
       enddo
    enddo
    !
    return
  end subroutine lagpolx
  !############################################################################
  !############################################################################
  subroutine lagpoly(u)
    !
    USE param
    USE complex_geometry
    USE decomp_2d
    USE variables
    !
    implicit none
    !
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
    integer                                            :: i,j,k
    real(mytype)                                       :: x,y,z
    integer                                            :: jy              !=   position du point "zappé"
    integer                                            :: jpif,jpol,nypif
    integer                                            :: jpoli,jpolf     !=   positions Initiales et Finales du POLynôme considéré
    real(mytype)                                       :: xpol,ypol,dypol   !|variables concernant les polynômes
    real(mytype),dimension(10)                         :: xa,ya           !|de   Lagrange. A mettre impérativement en
    integer                                            :: ia,na             !|double précision
    !
    do k=1,ysize(3)
       do i=1,ysize(1)
          if(nobjy(i,k).ne.0)then
             ia=0
             do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
                !1ère frontière
                nypif=npif
                ia=ia+1
                xa(ia)=yi(j,i,k)
                ya(ia)=0.
                if(yi(j,i,k).gt.0.)then!objet immergé
                   jy=1!jy=yi(j,i,k)/dy+1
                   do while(yp(jy).lt.yi(j,i,k))
                      jy=jy+1
                   enddo
                   jy=jy-1
                   jpoli=jy+1
                   if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
                   do jpif=1,nypif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=yp(jy-jpif)!(jy-1)*dy-jpif*dy
                         ya(ia)=u(i,jy-jpif,k)
                      else             !no zapping
                         xa(ia)=yp(jy-jpif+1)!(jy-1)*dy-(jpif-1)*dy
                         ya(ia)=u(i,jy-jpif+1,k)
                      endif
                   enddo
                else                   !objet semi-immergé
                   jpoli=1
                endif
                !2ème frontière
                nypif=npif
                ia=ia+1
                xa(ia)=yf(j,i,k)
                ya(ia)=0.
                if(yf(j,i,k).lt.yly)then!objet immergé
                   jy=1!jy=(yf(j,i,k)+dy)/dy+1
                   do while(yp(jy).lt.yf(j,i,k))  !there was a bug here yi<-->yf
                      jy=jy+1
                   enddo
                   jpolf=jy-1
                   if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
                   do jpif=1,nypif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=yp(jy+jpif)!(jy-1)*dy+jpif*dy
                         ya(ia)=u(i,jy+jpif,k)
                      else             !no zapping
                         xa(ia)=yp(jy+jpif-1)!(jy-1)*dy+(jpif-1)*dy
                         ya(ia)=u(i,jy+jpif-1,k)
                      endif
                   enddo
                else                   !objet semi-immergé
                   jpolf=ny
                endif
                !calcul du polynôme
                na=ia
                do jpol=jpoli,jpolf
                   xpol=yp(jpol)!dy*(jpol-1)
                   call polint(xa,ya,na,xpol,ypol,dypol)
                   u(i,jpol,k)=ypol
                enddo
                ia=0
             enddo
          endif
       enddo
    enddo
    !
    return
  end subroutine lagpoly
  !############################################################################
  !############################################################################
  subroutine lagpolz(u)
    !
    USE param
    USE complex_geometry
    USE decomp_2d
    USE variables
    !
    implicit none
    !
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
    integer                                            :: i,j,k
    real(mytype)                                       :: x,y,z
    integer                                            :: kz              !=   position du point "zappé"
    integer                                            :: kpif,kpol,nzpif
    integer                                            :: kpoli,kpolf     !=   positions Initiales et Finales du POLynôme considéré
    real(mytype)                                       :: xpol,ypol,dypol   !|variables concernant les polynômes
    real(mytype),dimension(10)                         :: xa,ya           !|de   Lagrange. A mettre imérativement en
    integer                                            :: ia,na             !|double précision
    !
    do j=1,zsize(2)
       do i=1,zsize(1)
          if(nobjz(i,j).ne.0)then
             ia=0
             do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple   (i,j)
                !1ère frontière
                nzpif=npif
                ia=ia+1
                xa(ia)=zi(k,i,j)
                ya(ia)=0.
                if(zi(k,i,j).gt.0.)then!objet immergé
                   kz=zi(k,i,j)/dz+1
                   kpoli=kz+1
                   if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
                   do kpif=1,nzpif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=(kz-1)*dz-kpif*dz
                         ya(ia)=u(i,j,kz-kpif)
                      else             !no zapping
                         xa(ia)=(kz-1)*dz-(kpif-1)*dz
                         ya(ia)=u(i,j,kz-kpif+1)
                      endif
                   enddo
                else                   !objet semi-immergé
                   kpoli=1
                endif
                !2ème frontière
                nzpif=npif
                ia=ia+1
                xa(ia)=zf(k,i,j)
                ya(ia)=0.
                if(zf(k,i,j).lt.zlz)then!objet immergé
                   kz=(zf(k,i,j)+dz)/dz+1
                   kpolf=kz-1
                   if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
                   do kpif=1,nzpif
                      ia=ia+1
                      if(izap.eq.1)then!zapping
                         xa(ia)=(kz-1)*dz+kpif*dz
                         ya(ia)=u(i,j,kz+kpif)
                      else             !no zapping
                         xa(ia)=(kz-1)*dz+(kpif-1)*dz
                         ya(ia)=u(i,j,kz+kpif-1)
                      endif
                   enddo
                else                   !objet semi-immergé
                   kpolf=nz
                endif
                !calcul du polynôme
                na=ia
                do kpol=kpoli,kpolf
                   xpol=dz*(kpol-1)
                   call polint(xa,ya,na,xpol,ypol,dypol)
                   u(i,j,kpol)=ypol
                enddo
                ia=0
             enddo
          endif
       enddo
    enddo
    !
    return
  end subroutine lagpolz
  !############################################################################
  !############################################################################
  subroutine polint(xa,ya,n,x,y,dy)
    !
    USE decomp_2d
    !
    implicit none
    !
    integer,parameter            :: nmax=30
    integer                      :: n,i,m,ns
    real(mytype)                 :: dy,x,y,den,dif,dift,ho,hp,w
    real(mytype),dimension(nmax) :: c,d
    real(mytype),dimension(n)    :: xa,ya
    ns=1
    dif=abs(x-xa(1))
    do i=1,n
       dift=abs(x-xa(i))
       if(dift.lt.dif)then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
    enddo
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          !         if(den.eq.0)read(*,*)
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       enddo
       if (2*ns.lt.n-m)then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    enddo
    return
  end subroutine polint
  !############################################################################
  !############################################################################
  subroutine lagpolx_phi(u)
    !================================================================================
    !
    !  SUBROUTINE: lagpolx_phi
    ! DESCRIPTION: Modified from the original Lagrangian interpolator.
    !              It imposes no-flux condition for the scalar field at the
    !              solid/fluid interface. Recomended only if the object's face
    !              is perpendicular with one of the spatial coordinates.
    !              Following Schuch (2020).
    !      AUTHOR: Felipe N. Schuch <felipe.schuch@edu.pucrs.br>
    !
    !================================================================================
    USE param
    USE complex_geometry
    USE decomp_2d
    USE variables
    !
    implicit none
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u
    integer                                            :: i,j,k
    real(mytype)                                       :: x,y,z
    integer                                            :: ix              != position du point "zappé"
    integer                                            :: ipif,ipol,nxpif
    integer                                            :: ipoli,ipolf     != positions Initiales et Finales du POLynôme considéré
    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en
    integer                                            :: ia,na           !|double précision
    !
    do k=1,xsize(3)
      do j=1,xsize(2)
        if(nobjx(j,k).ne.0)then
          ia=0
          do i=1,nobjx(j,k)          !boucle sur le nombre d'objets par (j,k)
            !1ère frontière
            nxpif=npif
            ! ia=ia+1
            ! xa(ia)=xi(i,j,k)
            ! ya(ia)=0.
            if(xi(i,j,k).gt.0.)then!objet immergé
              ix=xi(i,j,k)/dx+1
              ipoli=ix+1
              if(nxipif(i,j,k).lt.npif)nxpif=nxipif(i,j,k)
              do ipif=1,nxpif
                ia=ia+1
                if(izap.eq.1)then!zapping
                  xa(ia)=(ix-1)*dx-ipif*dx
                  ya(ia)=u(ix-ipif,j,k)
                  !
                  ia=ia+1
                  xa(ia)=two*xi(i,j,k)-((ix-1)*dx-ipif*dx)
                  ya(ia)=u(ix-ipif,j,k)
                else             !no zapping
                  xa(ia)=(ix-1)*dx-(ipif-1)*dx
                  ya(ia)=u(ix-ipif+1,j,k)
                  !
                  ia=ia+1
                  xa(ia)=two*xi(i,j,k)-((ix-1)*dx-(ipif-1)*dx)
                  ya(ia)=u(ix-ipif+1,j,k)
                endif
              enddo
            else                   !objet semi-immergé
              ipoli=1
              !
              ia=ia+1
              xa(ia)=zero
              ya(ia)=zero
            endif
            !2ème frontière
            nxpif=npif
            ! ia=ia+1
            ! xa(ia)=xf(i,j,k)
            ! ya(ia)=0.
            if(xf(i,j,k).lt.xlx)then!objet immergé
              ix=(xf(i,j,k)+dx)/dx+1
              ipolf=ix-1
              if(nxfpif(i,j,k).lt.npif)nxpif=nxfpif(i,j,k)
              do ipif=1,nxpif
                ia=ia+1
                if(izap.eq.1)then!zapping
                  xa(ia)=(ix-1)*dx+ipif*dx
                  ya(ia)=u(ix+ipif,j,k)
                  !
                  ia=ia+1
                  xa(ia)=two*xf(i,j,k)-((ix-1)*dx+ipif*dx)
                  ya(ia)=u(ix+ipif,j,k)
                else             !no zapping
                  xa(ia)=(ix-1)*dx+(ipif-1)*dx
                  ya(ia)=u(ix+ipif-1,j,k)
                  !
                  ia=ia+1
                  xa(ia)=two*xf(i,j,k)-((ix-1)*dx+(ipif-1)*dx)
                  ya(ia)=u(ix+ipif-1,j,k)
                endif
              enddo
            else                   !objet semi-immergé
              ipolf=nx
              !
              ia=ia+1
              xa(ia)=xlx
              ya(ia)=zero
            endif
            !calcul du polynôme
            na=ia
            do ipol=ipoli,ipolf
              xpol=dx*(ipol-1)
              call polint(xa,ya,na,xpol,ypol,dypol)
              u(ipol,j,k)=ypol
            enddo
            ia=0
          enddo
        endif
      enddo
    enddo
    !
    return
  end subroutine lagpolx_phi
  !############################################################################
  !############################################################################
  subroutine lagpoly_phi(u)
    !================================================================================
    !
    !  SUBROUTINE: lagpoly_phi
    ! DESCRIPTION: Modified from the original Lagrangian interpolator.
    !              It imposes no-flux condition for the scalar field at the
    !              solid/fluid interface. Recomended only if the object's face
    !              is perpendicular with one of the spatial coordinates.
    !              Following Schuch (2020).
    !      AUTHOR: Felipe N. Schuch <felipe.schuch@edu.pucrs.br>
    !
    !================================================================================
    USE param
    USE complex_geometry
    USE decomp_2d
    USE variables
    !
    implicit none
    !
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u
    integer                                            :: i,j,k
    real(mytype)                                       :: x,y,z
    integer                                            :: jy              != position du point "zappé"
    integer                                            :: jpif,jpol,nypif
    integer                                            :: jpoli,jpolf     != positions Initiales et Finales du POLynôme considéré
    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre impérativement en
    integer                                            :: ia,na           !|double précision
    !
    do k=1,ysize(3)
      do i=1,ysize(1)
        if(nobjy(i,k).ne.0)then
          ia=0
          do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
            !1ère frontière
            nypif=npif
            ! ia=ia+1
            ! xa(ia)=yi(j,i,k)
            ! ya(ia)=0.
            if(yi(j,i,k).gt.0.)then!objet immergé
              jy=1!jy=yi(j,i,k)/dy+1
              do while(yp(jy).lt.yi(j,i,k))
                jy=jy+1
              enddo
              jy=jy-1
              jpoli=jy+1
              if(nyipif(j,i,k).lt.npif)nypif=nyipif(j,i,k)
              do jpif=1,nypif
                ia=ia+1
                if(izap.eq.1)then!zapping
                  xa(ia)=yp(jy-jpif)!(jy-1)*dy-jpif*dy
                  ya(ia)=u(i,jy-jpif,k)
                  !
                  ia=ia+1
                  xa(ia)=two*yi(j,i,k)-yp(jy-jpif)
                  ya(ia)=u(i,jy-jpif,k)
                else             !no zapping
                  xa(ia)=yp(jy-jpif+1)!(jy-1)*dy-(jpif-1)*dy
                  ya(ia)=u(i,jy-jpif+1,k)
                  !
                  ia=ia+1
                  xa(ia)=two*yi(j,i,k)-yp(jy-jpif+1)
                  ya(ia)=u(i,jy-jpif+1,k)
                endif
              enddo
            else                   !objet semi-immergé
              jpoli=1
              !
              ia=ia+1
              xa(ia)=zero
              ya(ia)=zero
            endif
            !2ème frontière
            nypif=npif
            ! ia=ia+1
            ! xa(ia)=yf(j,i,k)
            ! ya(ia)=0.
            if(yf(j,i,k).lt.yly)then!objet immergé
              jy=1!jy=(yf(j,i,k)+dy)/dy+1
              do while(yp(jy).lt.yf(j,i,k))  !there was a bug here yi<-->yf
                jy=jy+1
              enddo
              jpolf=jy-1
              if(nyfpif(j,i,k).lt.npif)nypif=nyfpif(j,i,k)
              do jpif=1,nypif
                ia=ia+1
                if(izap.eq.1)then!zapping
                  xa(ia)=yp(jy+jpif)!(jy-1)*dy+jpif*dy
                  ya(ia)=u(i,jy+jpif,k)
                  !
                  ia=ia+1
                  xa(ia)=two*yi(j,i,k)-yp(jy+jpif)
                  ya(ia)=u(i,jy+jpif,k)
                else             !no zapping
                  xa(ia)=yp(jy+jpif-1)!(jy-1)*dy+(jpif-1)*dy
                  ya(ia)=u(i,jy+jpif-1,k)
                  !
                  ia=ia+1
                  xa(ia)=two*yi(j,i,k)-yp(jy+jpif-1)
                  ya(ia)=u(i,jy+jpif-1,k)
                endif
              enddo
            else                   !objet semi-immergé
              jpolf=ny
              !
              ia=ia+1
              xa(ia)=yly
              ya(ia)=zero
            endif
            !calcul du polynôme
            na=ia
            do jpol=jpoli,jpolf
              xpol=yp(jpol)!dy*(jpol-1)
              call polint(xa,ya,na,xpol,ypol,dypol)
              u(i,jpol,k)=ypol
            enddo
            ia=0
          enddo
        endif
      enddo
    enddo
    !
    return
  end subroutine lagpoly_phi
  !############################################################################
  !############################################################################
  subroutine lagpolz_phi(u)
    !================================================================================
    !
    !  SUBROUTINE: lagpolz_phi
    ! DESCRIPTION: Modified from the original Lagrangian interpolator.
    !              It imposes no-flux condition for the scalar field at the
    !              solid/fluid interface. Recomended only if the object's face
    !              is perpendicular with one of the spatial coordinates.
    !              Following Schuch (2020).
    !      AUTHOR: Felipe N. Schuch <felipe.schuch@edu.pucrs.br>
    !
    !================================================================================
    USE param
    USE complex_geometry
    USE decomp_2d
    USE variables
    !
    implicit none
    !
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u
    integer                                            :: i,j,k
    real(mytype)                                       :: x,y,z
    integer                                            :: kz              != position du point "zappé"
    integer                                            :: kpif,kpol,nzpif
    integer                                            :: kpoli,kpolf     != positions Initiales et Finales du POLynôme considéré
    real(mytype)                                       :: xpol,ypol,dypol !|variables concernant les polynômes
    real(mytype),dimension(10)                         :: xa,ya           !|de Lagrange. A mettre imérativement en
    integer                                            :: ia,na           !|double précision
    !
    do j=1,zsize(2)
      do i=1,zsize(1)
        if(nobjz(i,j).ne.0)then
          ia=0
          do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple (i,j)
            !1ère frontière
            nzpif=npif
            ! ia=ia+1
            ! xa(ia)=zi(k,i,j)
            ! ya(ia)=0.
            if(zi(k,i,j).gt.0.)then!objet immergé
              kz=zi(k,i,j)/dz+1
              kpoli=kz+1
              if(nzipif(k,i,j).lt.npif)nzpif=nzipif(k,i,j)
              do kpif=1,nzpif
                ia=ia+1
                if(izap.eq.1)then!zapping
                  xa(ia)=(kz-1)*dz-kpif*dz
                  ya(ia)=u(i,j,kz-kpif)
                  !
                  ia=ia+1
                  xa(ia)=two*zi(k,i,j)-((kz-1)*dz-kpif*dz)
                  ya(ia)=u(i,j,kz-kpif)
                else             !no zapping
                  xa(ia)=(kz-1)*dz-(kpif-1)*dz
                  ya(ia)=u(i,j,kz-kpif+1)
                  !
                  ia=ia+1
                  xa(ia)=two*zi(k,i,j)-((kz-1)*dz-(kpif-1)*dz)
                  ya(ia)=u(i,j,kz-kpif+1)
                endif
              enddo
            else                   !objet semi-immergé
              kpoli=1
              !
              ia=ia+1
              xa(ia)=zero
              ya(ia)=zero
            endif
            !2ème frontière
            nzpif=npif
            ! ia=ia+1
            ! xa(ia)=zf(k,i,j)
            ! ya(ia)=0.
            if(zf(k,i,j).lt.zlz)then!objet immergé
              kz=(zf(k,i,j)+dz)/dz+1
              kpolf=kz-1
              if(nzfpif(k,i,j).lt.npif)nzpif=nzfpif(k,i,j)
              do kpif=1,nzpif
                ia=ia+1
                if(izap.eq.1)then!zapping
                  xa(ia)=(kz-1)*dz+kpif*dz
                  ya(ia)=u(i,j,kz+kpif)
                  !
                  ia=ia+1
                  xa(ia)=two*zi(k,i,j)-((kz-1)*dz+kpif*dz)
                  ya(ia)=u(i,j,kz+kpif)
                else             !no zapping
                  xa(ia)=(kz-1)*dz+(kpif-1)*dz
                  ya(ia)=u(i,j,kz+kpif-1)
                  !
                  ia=ia+1
                  xa(ia)=two*zi(k,i,j)-((kz-1)*dz+(kpif-1)*dz)
                  ya(ia)=u(i,j,kz+kpif-1)
                endif
              enddo
            else                   !objet semi-immergé
              kpolf=nz
              !
              ia=ia+1
              xa(ia)=zlz
              ya(ia)=zero
            endif
            !calcul du polynôme
            na=ia
            do kpol=kpoli,kpolf
              xpol=dz*(kpol-1)
              call polint(xa,ya,na,xpol,ypol,dypol)
              u(i,j,kpol)=ypol
            enddo
            ia=0
          enddo
        endif
      enddo
    enddo
    !
    return
  end subroutine lagpolz_phi
  !############################################################################
end module ibm
