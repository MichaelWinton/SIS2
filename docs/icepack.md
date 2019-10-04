# SIS2 [Icepack](https://github.com/CICE-Consortium/Icepack/tree/Icepack1.0) usage

## General

SIS2 uses individual Icepack routines bypassing the Icepack driver.  Currently 
only ```icepack_therm_mushy.F90``` and ```icepack_mechred.F90``` are used directly.  Other
Icepack modules used by these are included in the SIS2 source unmodified from
their Icepack 1.0 versions.  ```icepack_mechred.F90``` is used without modification 
to implement ridging in SIS2.  Changes needed for SIS2 to use 
```icepack_therm_mushy.F90``` are documented below.  SIS2 currently uses the 
Briegleb/Light delta-Eddington code (in ``ice_shortwave_dEdd.F90``) copied from an 
earlier CICE version rather than the Icepack file, ```icepack_shortwave.F90```, 
to implement sea ice radiative transfer.

## ```icepack_therm_mushy.F90``` changes

SIS2 coupling modifications were made to icepack_therm_mushy.F90 for
three reasons:

1. **coupling to GFDL vertical atmosphere/sea-ice temperature solver**.  In
contrast to Icepack, The GFDL coupler calculates fluxes and their temperature 
sensitivities outside of the ice model and passes them in to be used unaltered
in the temperature timestepping.  To implement this in icepack_therm_mushy.F90,
optional arguments were added for the total flux and its temperature sensitivity 
and a logical variable to activate this behavior.  These variables needed to be
passed down to the picard solver where the flux calculation is skipped and its
output replaced with the passed-in values.

2. **coupling to the GFDL melt pond scheme**.  This scheme is under development
as of Fall 2019.  A novel aspect of the scheme is that the freezing sink of
pond mass is treated the same as seawater freezing onto the ice bottom.  If
pond is present, the surface temperature is held at the freezing point and
any imbalance between air-ice and conductive heat fluxes at the surface is
made up by melting, a pond mass source, or freezing, a sink.  The logical
argument ```gfdl_pond``` indicates the presence of pond and the activation of this
treatment.

3. **negligible snow and convergence**.  With the GFDL pond scheme it turned out
that negligible snow (ignored by Icepack) did not have negligible shortwave
absorption when the energy conservation convergence criterion was applied in the
iterative solver.  A change was made to add snow shortwave absorption to the surface 
energy budget when the snow is so small as to be neglected.

At the top level, these modifications are implemented as optional parameters to
subroutine ```temperature_changes_salinity```:

````
    logical, intent(in), optional :: gfdl_bc
    real(kind=dbl_kind), intent(in), optional :: fsurf0, dfsurf0_dT
    logical, intent(in), optional :: gfdl_pond
````
