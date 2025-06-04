#!/usr/bin/env python
# coding: utf-8

# # Convert HEALPix Zarr data to netCDF with variable subsets for TempestExtremes
# 
# # And compute the uivt and vivt based on ua, va, and hus
# # Attention: If the data of existing variables in an existing file needed to be updated, please delete the file first.
# #            Otherwise, nothing will be updated, as the script would do nothing if the variable already exists in the file.
# 
# ## Author:
# - Zhe Feng || zhe.feng@pnnl.gov 
# - Ziming Chen || ziming.chen@pnnl.gov


import os
import json
import datetime
import numpy as np
import xarray as xr
import intake
import easygems.healpix as egh

# =================== Configuration =======================
current_location = "NERSC"
s_Model = "icon_d3hp003" # "um_glm_n2560_RAL3p3" #  "nicam_gl11" # 
s_TimeRes = "PT6H"
zoom = 8
out_dir = f'/pscratch/sd/w/wcmca1/scream-cess-healpix/data4TE/Tmp/{s_Model}_{s_TimeRes}/'

# {Variable Names in the Raw Data : Output Variable Names}
varout_dict = {
    'time': 'time', 'lat': 'lat', 'lon': 'lon', 'pressure': 'lev', 'ua': 'ua', 'va': 'va', 'hus': 'hus',
    'orog': 'ELEV', 'pr': 'pr', 'prs': 'prs', 'ps': 'ps', 'psl': 'psl', 'uas': 'uas', 'vas': 'vas',
    'sfcWind': 'sfcWind', 'zg': 'zg', 'ta': 'ta', 'tas': 'tas', 'rlut': 'rlut'
}

# =============== Utility Functions =======================
def vertical_mass_integration(hus, ps, plev):
    '''
    The function computes the vertical mass integration (Eq.(1) in Chen et al. (2020 EF https://doi.org/ 10.1029/2019EF001435)) 
    for the input variable, hus
    Input para:
        hus: dataset for vertical mass integration in which the name of vertical coordinate must be "lev"
        ps: surface pressure datasest with same dimension sizes to hus except for the "lev" 
        plev: vertical pressure coordinate with same dimension sizes to hus lev dimension
    '''
    # Check vertical direction
    if not np.all(np.diff(plev.values) > 0):
        plev = plev[::-1]
    if not np.all(np.diff(hus["lev"].values) > 0):
        hus = hus.sel(lev=hus["lev"][::-1])
    #
    if ps.max() < 1200:
        ps = ps * 100
        ps.name = "Pa"
        ps.attrs["units"] = "Pa"
        ps.attrs["long_name"] = "Estimated surface pressure"
    if plev.max() < 1200:
        plev = plev * 100
        plev.name = "Pa"
        plev.attrs["units"] = "Pa"
        plev.attrs["long_name"] = "Pressure level"
    #
    ## Ziming Chen added: 2025-05-23
    # convert hus to use the same pressure coordinates as plev
    plev        = plev.assign_coords(lev=plev)
    hus         = hus.assign_coords(lev=plev)
    ## Ziming Chen added: 2025-05-23
    #
    plev_3d = plev * xr.ones_like(hus)
    ps_3d = ps * xr.ones_like(hus)
    #
    hus_masked = hus.where(plev_3d <= ps_3d)
    dp = np.gradient(plev.values) 
    dp = xr.DataArray(dp, coords={"lev": plev}, dims=["lev"])
    dp_3d = dp * xr.ones_like(hus)
    integrand = hus_masked * dp_3d / 9.8
    #
    return integrand.sum(dim="lev")

def clean_attrs_for_netcdf(ds):
    '''
    The function cleans the global attributes in dataset, ds, to avoid errors during outputing the data
    '''
    for key in list(ds.attrs):
        val = ds.attrs[key]
        if isinstance(val, dict):
            del ds.attrs[key]
        elif isinstance(val, bool):
            ds.attrs[key] = int(val)
    for var in ds.variables:
        for key, val in list(ds[var].attrs.items()):
            if isinstance(val, dict):
                ds[var].attrs[key] = json.dumps(val)
            elif isinstance(val, bool):
                ds[var].attrs[key] = int(val)
    return ds

def process_month_data(ds_month_year, year_val, month, original_history, new_history):
    if "sfcWind" not in ds_month_year and "uas" in ds_month_year and "vas" in ds_month_year:
        ds_month_year["sfcWind"] = np.sqrt(ds_month_year["uas"]**2 + ds_month_year["vas"]**2)
        ds_month_year["sfcWind"].attrs.update({"long_name": "surface Wind Speed", "units": "m s-1"})
    if all(k in ds_month_year for k in ["hus", "ua", "va"]):
        lev = ds_month_year["lev"]
        #
        if "ps" not in ds_month_year.data_vars:
            ## Ziming Chen added: 2025-05-23
            # The vertical pressure coordinate must range from 1000 hPa to 0 hPa
            # Otherwise, the ps would be the toppest level of the model
            if not np.all(np.diff(lev.values) < 0):
                lev = lev[::-1]
                lev = lev.assign_coords(lev=lev)
            hus_tmp = ds_month_year["hus"]
            if not np.all(np.diff(hus_tmp["lev"].values) < 0):
                hus_tmp = hus_tmp.sel(lev=hus_tmp["lev"][::-1])
                hus_tmp = hus_tmp.assign_coords(lev=lev)
            ##
            print("Vertical Pressure Coordinate")
            print(hus_tmp.coords['lev'].values)
            mask_valid = ~np.isnan(hus_tmp)
            # Use argmax to find the *first valid level from top* (lev should be sorted from high to low)
            first_valid_idx = mask_valid.argmax(dim="lev").compute()

            # Convert the index to the actual pres value
            ps         = lev.isel(lev=first_valid_idx)
            #
        else:
            ps         = ds_month_year["ps"]
        #
        ds_month_year["uivt"] = vertical_mass_integration(ds_month_year["hus"] * ds_month_year["ua"], ps, lev)
        ds_month_year["vivt"] = vertical_mass_integration(ds_month_year["hus"] * ds_month_year["va"], ps, lev)
        ds_month_year["uivt"].attrs.update({"long_name": "ZonalVapFlux", "units": "m^-1 s^-1 kg"})
        ds_month_year["vivt"].attrs.update({"long_name": "MeridionalVapFlux", "units": "m^-1 s^-1 kg"})
    else:
        print("Warning: Not Enough Variables (hus, ua, or va) to Compute the uivt or vivt")

    # ds_month_year.attrs = ds.attrs
    out_file = os.path.join(out_dir, f"{s_Model}_ivt_hp{zoom}_{s_TimeRes}.{year_val}{str(month).zfill(2)}.nc")

    if os.path.exists(out_file):
        existing_ds = xr.open_dataset(out_file)
        new_vars = set(ds_month_year.data_vars) - set(existing_ds.data_vars)
        existing_ds.close()
        for var in new_vars:
            v = clean_attrs_for_netcdf(ds_month_year[var])
            v.attrs.update({"history": f"{new_history}; {original_history}",
                            "source_model": s_Model,
                            "time_resolution": s_TimeRes,
                            "healpix_zoom": zoom,
                            "processing_script": "convert_zarr2nc_4TempestExtremes_compute_uivt_vivt.py"})
            v.to_netcdf(out_file, mode='a', engine='h5netcdf')
    else:
        ds_month_year = clean_attrs_for_netcdf(ds_month_year)
        ds_month_year.attrs.update({"history": f"{new_history}; {original_history}",
                                    "source_model": s_Model,
                                    "time_resolution": s_TimeRes,
                                    "healpix_zoom": zoom,
                                    "processing_script": "convert_zarr2nc_4TempestExtremes_compute_uivt_vivt.py"})
        ds_month_year.to_netcdf(out_file)
        print(f"Created new file: {out_file}")

# =============== Main Pipeline =======================
def main():
    os.makedirs(out_dir, exist_ok=True)

    cat = intake.open_catalog("https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml")[current_location]
    ds0 = cat[s_Model](zoom=zoom, time=s_TimeRes).to_dask()
    if "icon" in s_Model:
        ds0 = cat[s_Model](zoom=zoom, time=s_TimeRes, time_method='inst').to_dask()
    ds0 = ds0.pipe(egh.attach_coords)

    ds = {varout_dict[var]: ds0[var] for var in varout_dict if var in ds0.data_vars}
    ds = xr.Dataset(ds)
    ds.attrs = ds0.attrs
    for var in ds0.coords:
        if var in varout_dict and var != varout_dict[var]:
            ds = ds.rename_dims({var: varout_dict[var]})
            ds = ds.rename_vars({var: varout_dict[var]})
    ds['time'] = xr.decode_cf(ds).indexes['time']
    #
    original_history = ds.attrs.get('history', '')
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    new_history = f"{timestamp}: Data processed for TempestExtremes"

    for month, ds_month in ds.groupby('time.month'):
        for year, ds_month_year in ds_month.groupby('time.year'):
            year_val = ds_month_year.coords['time'].dt.year[0].values
            #
            # Filter data for just this month/year before processing
            start_date = f"{year}-{month:02d}-01"
            end_date = f"{year}-{month+1:02d}-01" if month < 12 else f"{year+1}-01-01"
                
            # Your processing code for this month
            process_month_data(ds_month_year, year_val, month, original_history, new_history)

if __name__ == "__main__":
    main()
