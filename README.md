# EPR_code
This repository contains some code I used during my Ph.D., including:
- Calculate winter extreme precipitation regimes (EPRs) in eastern North America as defined in [Low et al. 2022](https://doi.org/10.1175/MWR-D-21-0255.1). To reproduce this, run the following scripts in order:
  - precip_daily.py
  - daily_pcp_pctile.py
  - pcp_vol.py
  - pcp_vol_pctile.py
  - detect_EPRs.py

- Calculate moisture budget of eastern North America used to define EPRs:
  - moisture_budget.py

- Categorize a set of synoptic-scale weather patterns characterized by the 1000-500 hPa thickness and sea-level pressure (SLP) fields using self-organizing maps (SOM) and plot SOM node composites:
  - SOM.py

- Calculate maximum equivalent potential temperature (theta_e) in the 700-900 hPa layer:
  - max_700-900_thetae.py
 
**Note: ERA5 data is used for all scripts. Required data files containing which variables at which spatial and temporal resolution are listed at the beginning of each script.** The scripts also assume that the NumPy, Xarray, MetPy, Matplotlib, and Cartopy packages are installed.
