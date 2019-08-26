module ice_ridging_mod
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of SIS2.                                        *
!*                                                                     *
!* SIS2 is free software; you can redistribute it and/or modify it and *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* SIS2 is distributed in the hope that it will be useful, but WITHOUT *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! replaced T. Martin code with wrapper for Icepack ridging function - mw 1/18  !
!                                                                              !
! Prioritized to do list as of 6/4/19 (mw):                                    !
!                                                                              !
! 1) implement new snow_to_ocean diagnostic to record this flux.               !
! 2) implement ridging_rate diagnostics: ridging_shear, ridging_conv           !
! 3) implement "do_j" style optimization as in "compress_ice" or               !
!    "adjust_ice_categories" (SIS_transport.F90) if deemed necessary           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type
use SIS_tracer_registry, only : SIS_tracer_registry_type, SIS_tracer_type
use icepack_kinds
use icepack_itd, only: icepack_init_itd, cleanup_itd
use icepack_mechred, only: ridge_ice
use icepack_warnings, only: icepack_warnings_flush, icepack_warnings_aborted, &
                            icepack_warnings_setabort
use icepack_tracers, only: icepack_init_tracer_indices
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING

implicit none ; private

#include <SIS_memory_macros.h>
#include <MOM_memory_macros.h>

real, parameter :: Rho_ice = 905.0, Rho_snow = 330.0 ! can we eliminate these below?

! future namelist parameters?
integer (kind=int_kind), parameter :: &
     krdg_partic = 0  , & ! 1 = new participation, 0 = Thorndike et al 75
     krdg_redist = 0      ! 1 = new redistribution, 0 = Hibler 80

! e-folding scale of ridged ice, krdg_partic=1 (m^0.5)
real(kind=dbl_kind), parameter ::  mu_rdg = 3.0

public :: ice_ridging, ridge_rate

contains

!
! ice_ridging is a wrapper for the icepack ridging routine ridge_ice
!
subroutine ice_ridging(part_sz, mca_ice, mca_snow, mca_pond, &
                        mH_ice, mH_snow, mH_pond, TrReg, G, IG, dt, uc, vc, &
                        snow_to_ocn, enth_to_ocn, water_to_ocn)
  type(SIS_hor_grid_type),                       intent(inout) :: G
  type(ice_grid_type),                           intent(inout) :: IG
  real, dimension(SZI_(G),SZJ_(G),0:SZCAT_(IG)), intent(inout) :: part_sz
  ! the following are masses per grid area (i.e. already concentration weighted)
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),   intent(inout) :: mca_ice, mca_snow, mca_pond
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),   intent(inout) :: mH_ice, mH_snow, mH_pond
  type(SIS_tracer_registry_type),                pointer       :: TrReg
  real (kind=dbl_kind),                          intent(in)    :: dt
  real, dimension(SZIB_(G),SZJ_(G)),            intent(in)    :: uc ! velocities for
  real, dimension(SZI_(G),SZJB_(G)),            intent(in)    :: vc ! ridging rate
  real, dimension(SZI_(G),SZJ_(G)),             intent(inout)   :: snow_to_ocn ! accumulating
  real, dimension(SZI_(G),SZJ_(G)),             intent(inout)   :: enth_to_ocn ! accumulating
  real, dimension(SZI_(G),SZJ_(G)),             intent(inout)   :: water_to_ocn ! accumulating

! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end,
!                      in/out.
!  (inout)   mca_ice - The mass per unit grid-cell area of the ice in each
!                      category in H (often kg m-2).
!  (inout)   mca_snow - The mass per unit grid-cell area of the snow atop the
!                       ice in each category in H.
!  (inout)   mca_pond - The mass per unit grid-cell area of the melt ponds atop
!                       the ice in each category in H.
!  (inout)   mH_ice - The thickness of the ice in each category in H.
!  (inout)   mH_snow - The thickness of the snow atop the ice in each category
!                     in H.
!  (inout)   mH_pond - The thickness of the pond atop the ice in each category
!                     in H.
!  (inout)   TrReg - The registry of registered SIS ice and snow tracers.
!  (in)      G - The ocean's grid structure.
!  (in)      IG - The sea-ice-specific grid structure.
!  (in)      dt      - The amount of time over which the ice dynamics are to be
!                      advanced, in s.

  ! these strain metrics are calculated here from the velocities used for advection
  real :: sh_Dt ! sh_Dt is the horizontal tension (du/dx - dv/dy) including
                ! all metric terms, in s-1.
  real :: sh_Dd ! sh_Dd is the flow divergence (du/dx + dv/dy) including all
                ! metric terms, in s-1.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    sh_Ds       ! sh_Ds is the horizontal shearing strain (du/dy + dv/dx)
                ! including all metric terms, in s-1.


  integer :: i, j, k ! loop vars
  integer isc, iec, jsc, jec ! loop bounds
  integer :: halo_sh_Ds  ! The halo size that can be used in calculating sh_Ds.
  integer :: nIlyr, nSlyr, nCat ! # of ice/snow layers & thickness categories

  integer (kind=int_kind) :: &
     ndtd = 1  , & ! number of dynamics subcycles
     n_aero = 0, & ! number of aerosol tracers
     ntrcr = 0     ! number of tracer level


  real(kind=dbl_kind), dimension(0:IG%CatIce) :: hin_max ! category limits (m)

  real(kind=dbl_kind) :: &
     del_sh        , & ! shear strain measure
     rdg_conv = 0.0, & ! normalized energy dissipation from convergence (1/s)
     rdg_shear= 0.0    ! normalized energy dissipation from shear (1/s)
 
  real(kind=dbl_kind), dimension(IG%CatIce) :: &
     aicen, & ! concentration of ice
     vicen, & ! volume per unit area of ice          (m)
     vsnon    ! volume per unit area of snow         (m)
 
  ! ice tracers; ntr*(NkIce+NkSnow) guaranteed to be enough for all (intensive)
  real(kind=dbl_kind), dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow),IG%CatIce) :: trcrn

  real(kind=dbl_kind) :: aice0          ! concentration of open water

  integer (kind=int_kind), dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow)) :: &
     trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon (weighting to use)
     n_trcr_strata  ! number of underlying tracer layers

  real(kind=dbl_kind), dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow),3) :: &
     trcr_base      ! = 0 or 1 depending on tracer dependency
                    ! argument 2:  (1) aice, (2) vice, (3) vsno

  integer, dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow),IG%CatIce) :: &
     nt_strata      ! indices of underlying tracer layers

  type(SIS_tracer_type), dimension(:), pointer :: Tr=>NULL() ! SIS2 tracers
  integer :: m, n ! loop vars for tracer; n is tracer #; m is tracer layer

  nSlyr = IG%NkSnow
  nIlyr = IG%NkIce
  nCat  = IG%CatIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  ! copy strain calculation code from SIS_C_dynamics; might be a more elegant way ...
  !
  halo_sh_Ds = min(isc-G%isd, jsc-G%jsd, 2)
  do J=jsc-halo_sh_Ds,jec+halo_sh_Ds-1 ; do I=isc-halo_sh_Ds,iec+halo_sh_Ds-1
      ! This uses a no-slip boundary condition.
      sh_Ds(I,J) = (2.0-G%mask2dBu(I,J)) * &
          (G%dxBu(I,J)*G%IdyBu(I,J)*(uc(I,j+1)*G%IdxCu(I,j+1) - uc(I,j)*G%IdxCu(I,j)) + &
           G%dyBu(I,J)*G%IdxBu(I,J)*(vc(i+1,J)*G%IdyCv(i+1,J) - vc(i,J)*G%IdyCv(i,J)))
  enddo ; enddo

  ! set category limits; Icepack has a max on the largest, unlimited, category (why?)
  do k=1,nCat
     hin_max(k) = IG%mH_cat_bound(k)/(Rho_ice*IG%kg_m2_to_H)
  end do
  hin_max(nCat+1) = 1e5; ! not sure why this is needed, set big

  trcr_base = 0.0; n_trcr_strata = 0; nt_strata = 0; ! init some tracer vars
  ! When would we use icepack tracer "strata"?

  ! set icepack tracer index "nt_lvl" to (last) pond tracer so it gets dumped when
  ! ridging in ridge_ice (this is what happens to "level" ponds); first add up ntrcr;
  ! then set nt_lvl to ntrcr+1; could move this to an initializer - mw
  ntrcr = 0
  if (TrReg%ntr>0) then ! sum tracers
    Tr => TrReg%Tr_snow
    do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
      ntrcr = ntrcr + 1
    enddo ; enddo
    Tr => TrReg%Tr_ice
    do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
          ntrcr = ntrcr + 1
    enddo ; enddo
  endif
  call icepack_init_tracer_indices(nt_vlvl_in=ntrcr+1); ! pond will be last tracer

  do j=jsc,jec; do i=isc,iec;
  if ((G%mask2dT(i,j) .gt. 0.0) .and. (sum(part_sz(i,j,1:nCat)) .gt. 0.0)) then
  ! feed locations to Icepack's ridge_ice
  
    ! start like we're putting ALL the snow in the ocean
    snow_to_ocn(i,j) = snow_to_ocn(i,j) + sum(mca_snow(i,j,:));
    enth_to_ocn(i,j) = enth_to_ocn(i,j) + sum(mca_snow(i,j,:)*TrReg%Tr_snow(1)%t(i,j,:,1));
    water_to_ocn(i,j) = water_to_ocn(i,j) + sum(mca_pond(i,j,:));
    aicen = part_sz(i,j,1:nCat);

    if (sum(aicen) .le. 0.0) then ! no ice -> no ridging
      part_sz(i,j,0) = 1.0
    else

      ! set up ice and snow volumes
      vicen = mca_ice(i,j,:) /Rho_ice
      vsnon = mca_snow(i,j,:)/Rho_snow

      sh_Dt = (G%dyT(i,j)*G%IdxT(i,j)*(G%IdyCu(I,j) * uc(I,j) - &
                                       G%IdyCu(I-1,j)*uc(I-1,j)) - &
               G%dxT(i,j)*G%IdyT(i,j)*(G%IdxCv(i,J) * vc(i,J) - &
                                       G%IdxCv(i,J-1)*vc(i,J-1)))
      sh_Dd = (G%IareaT(i,j)*(G%dyCu(I,j) * uc(I,j) - &
                              G%dyCu(I-1,j)*uc(I-1,j)) + &
               G%IareaT(i,j)*(G%dxCv(i,J) * vc(i,J) - &
                              G%dxCv(i,J-1)*vc(i,J-1)))

      del_sh = sqrt(sh_Dd**2 + 0.25 * (sh_Dt**2 + &
                   (0.25 * ((sh_Ds(I-1,J-1) + sh_Ds(I,J)) + &
                            (sh_Ds(I-1,J) + sh_Ds(I,J-1))))**2 ) ) ! H&D eqn 9
      rdg_conv  = -min(sh_Dd,0.0)              ! energy dissipated by convergence ...
      rdg_shear = 0.5*(del_sh-abs(sh_Dd))      ! ... and by shear

!rdg_shear = 86400.0/100.0 ! ... for column model testing

      aice0 = 1.0+dt*sh_Dd-sum(aicen)
     
      ntrcr = 0
      if (TrReg%ntr>0) then ! load tracer array

        Tr => TrReg%Tr_snow
        do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
          ntrcr = ntrcr + 1
          trcrn(ntrcr,:) = Tr(n)%t(i,j,:,m)
          trcr_depend(ntrcr) = 2; ! 2 means snow-based tracer
          trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,3) = 1.0; ! 3rd index for snow
        enddo ; enddo

        Tr => TrReg%Tr_ice
        do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
          ntrcr = ntrcr + 1
          trcrn(ntrcr,:) = Tr(n)%t(i,j,:,m)
          trcr_depend(ntrcr) = 1; ! 1 means ice-based tracer
          trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,2) = 1.0; ! 2nd index for ice
        enddo ; enddo
       
      endif ! have tracers to load

      ! load pond on top of stack
      ntrcr = ntrcr + 1
      trcrn(ntrcr,:) = mH_pond(i,j,:)
      trcr_depend(ntrcr) = 0; ! 1 means ice area-based tracer
      trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,1) = 1.0; ! 1st index for ice area



      ! call Icepack routine; how are ponds treated?
      call ridge_ice (dt, ndtd, ncat, n_aero, nilyr, nslyr, ntrcr, hin_max,   &
                      rdg_conv, rdg_shear, aicen, trcrn, vicen, vsnon,        &
                      aice0, trcr_depend, trcr_base, n_trcr_strata, nt_strata,&
                      krdg_partic, krdg_redist, mu_rdg, tr_brine=.false.)

      if ( icepack_warnings_aborted() ) then
        call icepack_warnings_flush(0);
        call icepack_warnings_setabort(.false.)
        call SIS_error(WARNING,'icepack ridge_ice error');
      endif

      ! pop pond off top of stack
      mH_pond(i,j,:) = trcrn(ntrcr,:)
      mca_pond(i,j,:) = mH_pond(i,j,:)*aicen

      if (TrReg%ntr>0) then
        ! unload tracer array reversing order of load -- stack-like fashion

        Tr => TrReg%Tr_ice
        do n=TrReg%ntr,1,-1 ; do m=Tr(n)%nL,1,-1
          ntrcr = ntrcr - 1
          Tr(n)%t(i,j,:,m) = trcrn(ntrcr,:)
        enddo ; enddo
       
        Tr => TrReg%Tr_snow
        do n=TrReg%ntr,1,-1 ; do m=Tr(n)%nL,1,-1
          ntrcr = ntrcr - 1
          Tr(n)%t(i,j,:,m) = trcrn(ntrcr,:)
        enddo ; enddo

      endif ! have tracers to unload

      ! output: snow/ice masses/thicknesses
      do k=1,nCat
        if (aicen(k) > 0.0) then
          part_sz(i,j,k)  = aicen(k)
          mca_ice(i,j,k)  = vicen(k)*Rho_ice
          mH_ice(i,j,k)   = vicen(k)*Rho_ice/aicen(k)
          mca_snow(i,j,k) = vsnon(k)*Rho_snow
          mH_snow(i,j,k)  = vsnon(k)*Rho_snow/aicen(k)
        else
          part_sz(i,j,k) = 0.0
          mca_ice(i,j,k)  = 0.0
          mH_ice(i,j,k) = 0.0
          mca_snow(i,j,k) = 0.0
          mH_snow(i,j,k) = 0.0
        endif
        ! How to treat ponds?
      enddo

      ! negative part_sz(i,j,0) triggers compress_ice clean_up later
      part_sz(i,j,0) = 1.0 - sum(part_sz(i,j,1:nCat))

    endif
    ! subtract new snow/pond mass and energy on ice to sum net fluxes to ocean
    snow_to_ocn(i,j) = snow_to_ocn(i,j) - sum(mca_snow(i,j,:));
    enth_to_ocn(i,j) = enth_to_ocn(i,j) - sum(mca_snow(i,j,:)*TrReg%Tr_snow(1)%t(i,j,:,1));
    water_to_ocn(i,j) = water_to_ocn(i,j) - sum(mca_pond(i,j,:));

  endif; enddo; enddo ! part_sz, j, i

end subroutine ice_ridging

!TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ridge_rate - deformation rate or                                             !
!              total energy dissipation rate due to ridging                    !
!              (Flato and Hibler, 1995, JGR) or                                !
!              net area loss in riding (CICE documentation)                    !
!              depending on the state of the ice drift                         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function ridge_rate(del2, div) result (rnet)
  real, intent(in)  :: del2, div
  real              :: del, rnet, rconv, rshear
  !TOM> cs is now set in namelist:
!Niki: this was commented out
  real, parameter   :: cs=0.25 !(CICE documentation)

  del=sqrt(del2)

  rconv  = -min(div,0.0)           ! energy dissipated by convergence ...
  rshear = 0.5*(del-abs(div))      ! ... and by shear
  rnet   = rconv + cs*rshear       ! net energy contains only part of the
                                   !  shear energy as only a fraction is
           !  dissipated in the ridging process
  return
end function ridge_rate

end module ice_ridging_mod
