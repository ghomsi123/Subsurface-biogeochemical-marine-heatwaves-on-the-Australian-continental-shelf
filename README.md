# Subsurface-biogeochemical-marine-heatwaves-on-the-Australian-continental-shelf

[![DOI](https://zenodo.org/badge/1246554592.svg)](https://doi.org/10.5281/zenodo.20344424)

This repository contains processed glider datasets and MATLAB scripts for analyzing subsurface marine heatwaves along the Australian continental shelf. The data are organized by region:

1. Queensland (QLD) – northeastern Australia  
2. New South Wales (NSW) – southeastern Australia  
3. Western Australia South (WAS) – southwest Australia  
4. Tasmania (TAS) – eastern Tasmania  

The included MATLAB script, `glider_qc_regrid_vf.m`, provides tools for quality control, regridding, and processing of glider observations, enabling investigation of subsurface biogeochemical characteristics during marine heatwaves.

---

# Study Reference

Mawren, D., Araujo, J., Le Gendre, R., Benthuysen, J. A., Ghomsi, F. E. K., Saranya, J. S., and Schaeffer, A. (2025): *Gliding through marine heatwaves: Subsurface biogeochemical characteristics on the Australian continental shelf*, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-6045.

© Author(s) 2025. Distributed under the Creative Commons Attribution 4.0 License.

---

# Repository Contents

1. `Archive_gliders_18092025-20260226T231052Z-1-001.zip`  
   Contains glider profiles from NSW, QLD, WAS, and TAS.

2. `glider_qc_regrid_vf.m`  
   MATLAB processing script for glider data.

3. `glider_mean_profiles_mhw.nc`  
   Mean vertical profiles separated by region, season, and Marine Heatwave (MHW) conditions.

---

# Glider Mean Profiles (MHW vs Non-MHW)

This dataset contains mean vertical ocean profiles derived from glider missions around Australia, separated by region, season, and Marine Heatwave (MHW) conditions.

## Data File

`glider_mean_profiles_mhw.nc`

---

## Regions

- NSW (New South Wales)
- QLD (Queensland)
- TAS (eastern Tasmania)
- WAS (Western Australia South)

---

## Seasons

- winter (Jun–Aug)
- spring (Sep–Nov)
- summer (Dec–Feb)
- autumn (Mar–May)

---

## MHW Conditions

- `no_mhw` = background conditions
- `mhw` = marine heatwave conditions (severity > 1)

---

## Variables

- `TEMP` – temperature
- `PSAL` – salinity
- `CPHL` – chlorophyll
- `DOX2` – dissolved oxygen

---

## Data Structure

The dataset dimensions are organized as:

```text
(variable, region, condition, depth)
```

where `condition = season + MHW class`, for example:

- `winter_no_mhw`
- `summer_mhw`

---

## Example in Python

```python
import xarray as xr

ds = xr.open_dataset(
    "/Volumes/CR4D/final_netcdf_files_gliders/glider_mean_profiles_mhw.nc"
)

temp = ds["mean_profile"].sel(
    variable="TEMP",
    region="NSW",
    condition="summer_mhw"
)
```

---

# Notes

Glider profile positions for each region are illustrated in the publication: QLD, NSW, WAS, and TAS.

This repository provides regional subsets of data for reproducibility; the full dataset and high-resolution bathymetry are available from the original sources.

The scripts and datasets support exploration of subsurface biogeochemical metrics during marine heatwaves, as discussed in the referenced study.
