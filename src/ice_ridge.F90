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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type
use SIS_tracer_registry, only : SIS_tracer_registry_type, SIS_tracer_type
use icepack_kinds
use icepack_itd, only: icepack_init_itd, cleanup_itd
use icepack_mechred, only: ridge_ice

implicit none ; private

#include <SIS_memory_macros.h>
#include <MOM_memory_macros.h>

real, parameter :: Rho_ice = 905.0, Rho_snow = 330.0 ! can we eliminate these below?

! future namelist parameters?
integer (kind=int_kind), parameter :: &
     krdg_partic = 1  , & ! 1 = new participation, 0 = Thorndike et al 75
     krdg_redist = 1      ! 1 = new redistribution, 0 = Hibler 80

! e-folding scale of ridged ice, krdg_partic=1 (m^0.5)
real(kind=dbl_kind), parameter ::  mu_rdg = 3.0

public :: ice_ridging

contains

!
! need to get right and document input/output in terms of intensive/extensive quantities
!
subroutine ice_ridging(part_sz, mca_ice, mca_snow, mca_pond, &
                        mH_ice, mH_snow, mH_pond, TrReg, G, IG, dt, snow_to_ocn)
  type(SIS_hor_grid_type),                       intent(inout) :: G
  type(ice_grid_type),                           intent(inout) :: IG
  real, dimension(SZI_(G),SZJ_(G),0:SZCAT_(IG)), intent(inout) :: part_sz
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),   intent(inout) :: mca_ice, mca_snow, mca_pond
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),   intent(inout) :: mH_ice, mH_snow, mH_pond
  type(SIS_tracer_registry_type),                pointer       :: TrReg
  real (kind=dbl_kind),                          intent(in)    :: dt
  real, dimension(SZI_(G),SZJ_(G)),              intent(out)   :: snow_to_ocn

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

  integer :: i, j, k ! loop vars
  integer isc, iec, jsc, jec ! loop bounds
  integer :: nIlyr, nSlyr, nCat ! # of ice/snow layers & thickness categories

  integer (kind=int_kind) :: &
     ndtd = 1  , & ! number of dynamics subcycles
     n_aero = 0, & ! number of aerosol tracers
     ntrcr = 0     ! number of tracer level


  real(kind=dbl_kind), dimension(0:IG%CatIce) :: hin_max ! category limits (m)

  real(kind=dbl_kind) :: &
     rdg_conv = 0.0, & ! normalized energy dissipation from convergence (1/s)
     rdg_shear= 0.0    ! normalized energy dissipation from shear (1/s)
 
  real(kind=dbl_kind), dimension(IG%CatIce) :: &
     aicen, & ! concentration of ice
     vicen, & ! volume per unit area of ice          (m)
     vsnon    ! volume per unit area of snow         (m)
 
  ! ice tracers; ntr*(NkIce+NkSnow) guaranteed to be enough for all
  real(kind=dbl_kind), dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow),IG%CatIce) :: trcrn

  real(kind=dbl_kind) :: aice0          ! concentration of open water

  integer (kind=int_kind), dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow)) :: &
     trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
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

  ! set category limits; Icepack has a max on the largest, unlimited, category (why?)
  do k=1,nCat
     hin_max(k) = IG%mH_cat_bound(k)/(Rho_ice*IG%kg_m2_to_H)
  end do
  hin_max(nCat+1) = 1e5; ! not sure why this is needed, set big

  trcr_base = 0.0; n_trcr_strata = 0; nt_strata = 0; ! init some tracer vars
  ! When will we use icepack tracer "strata"?

  do j=jsc,jec; do i=isc,iec; ! feed locations to Icepack's ridge_ice
  
    aice0 = part_sz(i,j,0); aicen = part_sz(i,j,1:nCat);
    snow_to_ocn(i,j) = sum(mca_snow(i,j,:)); ! record original snow mass
    if (sum(aicen) .le. 0.0) then ! no ice -> no ridging
      part_sz(i,j,0) = 1.0
    else

      ! set up ice and snow volumes
      vicen = mca_ice(i,j,:) /Rho_ice
      vsnon = mca_snow(i,j,:)/Rho_snow
      ! later need to get these from dynamics
      rdg_conv = max((aice0+sum(aicen)-1)/dt,0.0);
      rdg_shear = 0.0/dt
     
      if (TrReg%ntr>0) then ! load tracer array
        ntrcr = 0

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

      ! call Icepack routine; how are ponds treated?
      ! call Icepack routine; snow gets dumped to ocean; how do we track this?
      call ridge_ice (dt, ndtd, ncat, n_aero, nilyr, nslyr, ntrcr, hin_max,   &
                      rdg_conv, rdg_shear, aicen, trcrn, vicen, vsnon,        &
                      aice0, trcr_depend, trcr_base, n_trcr_strata, nt_strata,&
                      krdg_partic, krdg_redist, mu_rdg, tr_brine=.false.)

      if (TrReg%ntr>0) then
        ! unload tracer array reversing order of load in stack-like fashion
        ntrcr = ntrcr + 1

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
          part_sz(i,j,k) = aicen(k)
          mca_ice(i,j,k)  = vicen(k)*Rho_ice
          mH_ice(i,j,k) = mca_ice(i,j,k)/part_sz(i,j,k)
          mca_snow(i,j,k) = vsnon(k)*Rho_snow
          mH_snow(i,j,k) = mca_snow(i,j,k)/part_sz(i,j,k)
        else
          part_sz(i,j,k) = 0.0
          mca_ice(i,j,k)  = 0.0
          mH_ice(i,j,k) = 0.0
          mca_snow(i,j,k) = 0.0
          mH_snow(i,j,k) = 0.0
        endif
        ! How to treat pond?  Maybe load pond volume as Icepack area tracer
        ! so they get reduced fractionally with compression in ridging (?)
      enddo

    endif
    ! subtract new snow mass to obtain flux to ocean
    snow_to_ocn(i,j) = snow_to_ocn(i,j) - sum(mca_snow(i,j,:));
  enddo; enddo ! j, i

end subroutine ice_ridging

end module ice_ridging_mod
