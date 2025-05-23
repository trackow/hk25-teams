
import numpy as np
import pandas as pd
import xarray as xr
import scipy.stats

import healpy

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

from tools.LvL import LvL

# 10x10 deg stcu domains from Klein and Hartmann 1993

domains10x10 = {
    "peruvian":     np.array([-90, -80, -20, -10]) ,
    "namibian":     np.array([0, 10, -20, -10]),
    "californian":  np.array([-130, -120, 20, 30]),
    "canarian":     np.array([-35, -25, 15, 25])
}


# Integrate a variable in vertical along pressure dimension

def integrate_wrt_pressure(da):
    if 'pressure' in da.dims:
        if da['pressure'].attrs['units'] == 'hPa':
            pfactor = 100
        else:
            pfactor = 1
        return da.sortby('pressure').integrate('pressure')/9.81*pfactor


# Max velocity in a column below a given pressure level [in Pa]

def reduce_below(da,plevel=900e2,fun=np.max):
    if 'pressure' in da.dims:
        if da['pressure'].attrs['units'] == 'hPa':
            pfactor = 100
        else:
            pfactor = 1
        # return da.isel(pressure=(da.pressure*pfactor>=plevel)).max(dim='pressure')
        return da.isel(pressure=(da.pressure*pfactor>=plevel)).reduce(fun,dim='pressure')


# Remapping function from easygems

def get_nn_lon_lat_index(nside, lons, lats):
    lons2, lats2 = np.meshgrid(lons, lats)
    return xr.DataArray(
        healpy.ang2pix(nside, lons2, lats2, nest=True, lonlat=True),
        coords=[("lat", lats), ("lon", lons)],
    )


# Read EarthCare stats table prepared by Johanna

def read_earthcare_csv(file):
    earthcare = pd.read_csv(file, index_col='date_time',
                            usecols=('date_time','lwp_mean','lwp_std','lwp_skew',
                                     'cloud_cover','LvL_KS1','LvL_KS2'),
                            parse_dates=True, date_format='%Y%m%dT%H%M%SZ') \
                .rename_axis('time').sort_values(by='time', ascending=True)  \
                .rename(columns={'lwp_skew':'lwp_skw','LvL_KS1':'ks_cloud','LvL_KS2':'ks_void'})

    earthcare.index = earthcare.index.map(lambda dt: dt.replace(year=2020))
    earthcare['lwp_hom'] = earthcare['lwp_mean']/earthcare['lwp_std']
    
    return xr.Dataset.from_dataframe(earthcare)


# Compute basic stats of a given variable

def basic_stats(ds,var):
    ds[var+'_mean'] = ds[var].mean(dim='cell')
    ds[var+'_std']  = ds[var].std(dim='cell')
    ds[var+'_skw']  = ds[var].reduce(scipy.stats.skew,dim='cell',nan_policy='omit')
    ds[var+'_hom']  = (ds[var+'_mean']/ds[var+'_std'])
    return ds


# Compute LvL metrics and assign them properly into dataset

def LvL2dataset(ds,dim='time'):
    
    if dim not in ds.dims:
        ds = ds.expand_dims(dim)
    Nt = ds.dims[dim]
    Lmax = max(ds.dims[d] for d in ('lat','lon'))

    ds = ds.assign_coords( chord = ('chord', np.arange(1,Lmax+1,1)) ) \
           .assign( ks_cloud     = (dim, np.zeros(Nt)),
                    ks_void      = (dim, np.zeros(Nt)),
                    cnt_cloud    = ([dim,'chord'], np.zeros((Nt,Lmax))),
                    cnt_cloud_r  = ([dim,'chord'], np.zeros((Nt,Lmax))),
                    cnt_void     = ([dim,'chord'], np.zeros((Nt,Lmax))),
                    cnt_void_r   = ([dim,'chord'], np.zeros((Nt,Lmax)))
                  )

    for i, t in enumerate(ds[dim]):
        ksc, ksv, c, cr, v, vr = LvL(ds['cloud_mask'].sel({dim:t}).data)
        Nc = min(len(c),Lmax)
        Nv = min(len(v),Lmax)
        ds['ks_cloud'][i],      ds['ks_void'][i]         = ksc,    ksv
        ds['cnt_cloud'][i,:Nc], ds['cnt_cloud_r'][i,:Nc] = c[:Nc], cr[:Nc]
        ds['cnt_void'][i,:Nv],  ds['cnt_void_r'][i,:Nv]  = v[:Nv], vr[:Nv]

    return ds


# Collapse LvL stats in groupby groups

def LvL2groupby(ds):
    lvl = ds.sum(dim='time')[['cnt_cloud','cnt_void']] \
            .merge( ds.mean(dim='time')[['cloud_cover','ks_cloud','ks_void']] )
    lvl['Nt'] = ds.count(dim='time')['cloud_cover']
    
    Npix = ds.dims['lon']*ds.dims['lat']
    lvl['cnt_cloud_r'] = lvl['Nt']*Npix*(1-lvl['cloud_cover'])**2*lvl['cloud_cover']**lvl.chord
    lvl['cnt_void_r']  = lvl['Nt']*Npix*lvl['cloud_cover']**2*(1-lvl['cloud_cover'])**lvl.chord

    return lvl


# Plot LvL distributions of coud and void chord lengths

def plot_LvL_dist(ax,ds):
    if 'resolution' in ds.attrs:
        x = ds.chord*ds.attrs['resolution']
        ax.set_xlabel('Chord length [km]')
    else:
        x = ds.chord
        ax.set_xlabel('Chord length [pix]')
    ax.plot(x, ds['cnt_cloud']/np.sum(ds['cnt_cloud']), label='cloud obs')
    ax.plot(x, ds['cnt_void']/np.sum(ds['cnt_void']), label='void obs')
    ax.plot(x, ds['cnt_cloud_r']/np.sum(ds['cnt_cloud_r']), '--', label='cloud random')
    ax.plot(x, ds['cnt_void_r']/np.sum(ds['cnt_void_r']), '--', label='cloud observed')
    ax.grid()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title(f"CF={ds['cloud_cover']:.2f} KScloud={ds['ks_cloud']:.2f} KSvoid={ds['ks_void']:.2f}")


# Draw coastlines and grid on a map

def draw_map(ax,map_domain):
    ax.set_extent(map_domain, crs=ccrs.PlateCarree())
    ax.add_feature(cf.NaturalEarthFeature('physical', 'land', '10m'),
               facecolor='none', edgecolor='red', linewidth=1)
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    return gl


# Add colorbar and annnotate variable name and datetime on a map

def annotate_map(ax,da,im):
    for att in ['standard_name','name','long_name']:
        if att in da.attrs:
            name = da.attrs[att]
        else:
            name = da.name
    if 'units' in da.attrs:
        unit = f" [{da.attrs['units']:s}]"
    else:
        unit = ""
        
    cb = plt.colorbar(im, ax=ax, shrink=0.9, aspect=30, pad=0.02, label=name+unit)
    
    datestr = da.time.values.astype('datetime64[h]').item().strftime('%Y-%m-%dT%H')
    if 'crs' in da.coords and 'healpix_nside' in da.crs.attrs:
        zoom = np.log2(da.crs.healpix_nside).astype(int)
        ax.set_title(datestr+f"  zoom={zoom:d}")
    else:
        ax.set_title(datestr)