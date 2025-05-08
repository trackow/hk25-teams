import xarray as xr
import numpy as np
import healpy as hp
import easygems.healpix as egh

# ------------------------------------------------------------------------------
# Functions to determine triggering area of MCSs
# ------------------------------------------------------------------------------
def add_circular_trigger_areas(
        mcs_trigger_locs: xr.DataArray,
        RADII: np.ndarray,
        hp_grid: str,
        ) -> xr.DataArray:
    """
    Add circular trigger areas to the MCS trigger locations.

    This function calculates the circular trigger areas around the given 
    trigger locations based on the specified radii and Healpix grid 
    configuration. The trigger areas are represented as indices of the 
    Healpix cells within the specified radius around each trigger location.

    Parameters
    ----------
    mcs_trigger_locs : xr.DataArray
        DataArray containing the MCS trigger locations. It must include 
        the 'trigger_idx' coordinate, which specifies the Healpix 
        cell indices of the trigger locations.
    RADII : np.ndarray
        Array of radii (in the same units as the Healpix grid) for which 
        the trigger areas are to be calculated.
    hp_grid : str
        String specifying the Healpix grid configuration. This is used to 
        determine the nside and nesting scheme of the grid.

    Returns
    -------
    xr.DataArray
        Updated DataArray with an additional coordinate 'trigger_area_idxs', 
        which contains the indices of the Healpix cells within the trigger 
        areas for each radius in RADII.

    Notes
    -----
    - The Healpix grid attributes (nside and nest) are determined using 
        helper functions from the `egh` module.
    - The trigger areas are initialized and stored in a new coordinate 
        'trigger_area_idxs', which is a 3D array with dimensions corresponding 
        to the number of tracks, the maximum number of cells in any trigger 
        area, and the number of radii in RADII.
    """
    # Determine attributes of the healpix grid
    nside = egh.get_nside(hp_grid)
    nest = True if egh.get_nest(hp_grid) else False

    for i, radius in enumerate(RADII):
        # Get all pixels within radius around triggering location, i.e. get the
        # 'triggering area'
        trigger_area_idxs = [
            _get_trigger_area_idxs(nside, nest, cell_idx, radius)
            for cell_idx in mcs_trigger_locs['trigger_idx'].values
            ]

        if i == 0:
            # Initialize data_array to save cell indices of the trigger area
            mcs_trigger_locs['trigger_area_idxs'] = \
                _init_trigger_area_idxs_array(
                    mcs_trigger_locs['tracks'].values, trigger_area_idxs,
                    RADII
                    )

        for j, idxs in enumerate(trigger_area_idxs):
            mcs_trigger_locs['trigger_area_idxs'][j, :idxs.shape[0], i] = idxs

    return mcs_trigger_locs


def _get_trigger_area_idxs(
        nside: int,
        nest: bool,
        cell_idx: int,
        radius_deg: float,
        ) -> np.ndarray:
    """
    Get the indices of the trigger area for a given cell index and radius.
    The trigger area is defined as the area within a certain radius from the
    healpix cell that marks the MCS center at its triggering time.
    Parameters
    ----------
    nside : int
        The nside parameter for the HEALPix map.
    nest : bool
        If True, use nested indexing. If False, use ring indexing.
    cell_idx : int
        The index of the cell for which to get the trigger area indices.
    radius_deg : float
        The radius in degrees for the trigger area.
    Returns
    -------
    np.ndarray
        The indices of the trigger area.
    """
    return hp.query_disc(
        nside, hp.pix2vec(nside, cell_idx, nest=nest), np.radians(radius_deg),
        inclusive=False, nest=nest,
    )

def _init_trigger_area_idxs_array(
        tracks: np.ndarray,
        trigger_area_idxs: list,
        radii: np.ndarray,
        ) -> xr.DataArray:
    """
    Initialize a DataArray to store trigger area indices.
    This function creates an xarray DataArray with dimensions corresponding to 
    tracks, cells, and radius degrees. The array is initialized with NaN values 
    and is intended to store indices of trigger areas for each track and radius.
    Parameters
    ----------
    tracks : np.ndarray
        An array representing the tracks for which the trigger area indices 
        are being initialized.
    trigger_area_idxs_list : list
        A list of arrays, where each array contains the indices of trigger areas 
        for a specific track.
    radii_degree : np.ndarray
        An array of radius values (in degrees) for which the trigger area
        indices are calculated.
    Returns
    -------
    xr.DataArray
        A DataArray with dimensions ['tracks', 'cell', 'radius'],
        initialized with NaN values, and coordinates corresponding to the input
        parameters.
    """    
    cells = np.arange(
        0, np.array([l.shape[0] for l in trigger_area_idxs]).max()
        )

    trigger_area_idxs_array = xr.DataArray(
        data=np.full(
            shape=(tracks.size, cells.size, radii.size),
            fill_value=np.nan,
        ),
        dims=['tracks', 'cell', 'radius'],
        coords={'tracks': tracks, 'cell': cells, 'radius': radii},
        )
    trigger_area_idxs_array['radius'].attrs['units'] = 'degree'
    return trigger_area_idxs_array


# ------------------------------------------------------------------------------
# Functions to subsample MCSs
# ------------------------------------------------------------------------------
def remove_land_triggers(
        mcs_trigger_locs: xr.DataArray,
        ocean_mask: xr.DataArray,
        ) -> xr.DataArray:
    """
    Removes land-based triggers from the given MCS (Mesoscale Convective
    Systems) trigger locations based on an ocean mask and returns only the
    all ocean-based triggers.

    Parameters:
    -----------
    mcs_trigger_locs : xr.DataArray
        An xarray DataArray containing the trigger locations for MCS. It is
        expected to have a dimension 'tracks' and a variable
        'trigger_idx' representing the Healpix cell index of the
        triggers.
    ocean_mask : xr.DataArray
        An xarray DataArray representing the ocean mask. It is used to determine
        whether the trigger areas are entirely over the ocean.

    Returns:
    --------
    xr.DataArray
        A filtered xarray DataArray containing only the MCS trigger locations
        that are entirely over the ocean. The returned DataArray has the same
        structure as the input but excludes land-based triggers.
    """
    mcs_trigger_locs['is_trigger_area_all_ocean'] = (
        'tracks', _is_max_trigger_area_all_ocean(
            mcs_trigger_locs, ocean_mask
            )
        )
    mcs_trigger_locs_ocean = mcs_trigger_locs.where(
        mcs_trigger_locs['is_trigger_area_all_ocean'] == True, drop=True
        )
    mcs_trigger_locs_ocean = mcs_trigger_locs_ocean.drop(
        'is_trigger_area_all_ocean'
        )
    mcs_trigger_locs_ocean['trigger_idx'] = \
        mcs_trigger_locs_ocean['trigger_idx'].astype(int)
    
    return mcs_trigger_locs_ocean


def _is_max_trigger_area_all_ocean(
        mcs_trigger_locs: xr.DataArray,
        ocean_mask: xr.DataArray
        ) -> xr.DataArray:
    is_max_trigger_area_all_ocean = []
    for track in mcs_trigger_locs['tracks']:
        # Get the trigger area cell indices for current track without filling
        # NaNs
        trigger_area_idxs = mcs_trigger_locs['trigger_area_idxs'].sel(
            radius=mcs_trigger_locs['radius'].max().values,
            tracks=track,
            )
        trigger_area_idxs = trigger_area_idxs[~np.isnan(trigger_area_idxs)]
        
        # Check if all cells in the trigger area are ocean
        is_max_trigger_area_all_ocean.append(
            all(~np.isnan(ocean_mask.sel(cell=trigger_area_idxs)))
            )
    return is_max_trigger_area_all_ocean


# ------------------------------------------------------------------------------
# Functions to subsample variables in MCS trigger area
# ------------------------------------------------------------------------------
def get_var_in_trigger_area(
        mcs_trigger_locs: xr.DataArray,
        data_field: xr.DataArray,
        ) -> xr.DataArray:
    """
    Select the values of a field in the triggering area of MCSs.

    This function selects the values of a given field within the 
    triggering area of Mesoscale Convective Systems (MCSs) for each track 
    and radius.

    Parameters
    ----------
    mcs_trigger_locs : xr.DataArray
        DataArray containing the MCS trigger locations. It must include 
        the 'trigger_area_idxs' coordinate, which specifies the Healpix 
        cell indices of the trigger areas.
    data_field : xr.DataArray
        DataArray containing the simulation data. It must include a 
        'cell' dimension corresponding to the Healpix grid.

    Returns
    -------
    xr.DataArray
        DataArray containing the mean values of the variable in the 
        triggering area for each track and radius. The dimensions are 
        ['tracks', 'radius'].
    """
    
    var_in_trigger_area = _init_var_in_trigger_area(mcs_trigger_locs)

    for j, track in enumerate(mcs_trigger_locs['tracks']):
        start_basetime = mcs_trigger_locs.sel(tracks=track)['start_basetime']
        var_before_triggering = data_field.sel(
            time=start_basetime, method='pad',
            ).compute()

        for i, radius in enumerate(mcs_trigger_locs['radius']):
            # Get trigger area cell indicees
            trigger_area_idxs = mcs_trigger_locs['trigger_area_idxs']\
                .sel(tracks=track, radius=radius)
            trigger_area_idxs = trigger_area_idxs[~np.isnan(trigger_area_idxs)]

            var_in_trigger_area[j, :trigger_area_idxs.shape[0], i] = \
                var_before_triggering.sel(cell=trigger_area_idxs).data
    
    return var_in_trigger_area


def _init_var_in_trigger_area(
        mcs_trigger_locs: xr.DataArray,
        ) -> xr.DataArray:
    tracks = mcs_trigger_locs['tracks']
    cells = mcs_trigger_locs['cell']
    radii = mcs_trigger_locs['radius']
    return xr.DataArray(
        data=np.full(
            shape=(tracks.size, cells.size, radii.size), fill_value=np.nan,
            ),
        dims=['tracks', 'cell', 'radius'],
        coords={'tracks': tracks, 'cell': cells, 'radius': radii},
        )