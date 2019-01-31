#!/usr/bin/env python
"""
This script outputs ECHAM-HAM fields in a CMOR compliant fashion.

See https://code.mpimet.mpg.de/projects/afterburner/wiki for a table of afterburner variable codes

(c) Duncan watson-parris 2017
"""
import argparse
from cmor_var import cmor_var, select_vars
from calculate_EF_ACI import *
import cis
from cf_units import Unit
from functools import partial
import itertools

unitless = Unit(None)

ACCELERATION_DUE_TO_GRAVITY = 9.8


def calc_od440aer(infile, product):
    tau_550, ang_550 = cis.read_data_list(infile, ["TAU_2D_550nm", "ANG_550nm_865nm"], product)
    return tau_550 * (550 / 440) ** ang_550.data


def calc_od870aer(infile, product):
    tau_865, ang_550 = cis.read_data_list(infile, ["TAU_2D_865nm", "ANG_550nm_865nm"], product)
    return tau_865 * (865 / 870) ** ang_550.data


def calc_od550lt1aer(infile, product):
    from cis.data_io.gridded_data import make_from_cube
    fine_mode_taus = cis.read_data_list(infile, ['TAU_MODE_K?_550nm', 'TAU_MODE_A?_550nm'], product)
    total_modes = make_from_cube(sum(fine_mode_taus))
    od550, = total_modes.collapsed('hybrid level at layer midpoints', 'sum')
    return od550


def calc_ext550aer(infile, product):
    """
    I want units of m-1, Tau *should* be in units of 1 so I just need to divide by layer thickness
    :param infile:
    :param product:
    :return:
    """
    from cis.data_io.gridded_data import make_from_cube
    from utils import get_stream_file
    deltaz = cis.read_data(get_stream_file(infile, 'vphysc'), 'grheightm1', product)
    tau_3d = make_from_cube(sum(cis.read_data_list(infile, ['TAU_MODE_??_550nm'], product)))
    ext_3d = tau_3d / deltaz
    return ext_3d


def sum_variables(infile, variables, product):
    return sum(cis.read_data_list(infile, variables, product))


def diff_variables(infile, var1, var2, product):
    """
    var1 - var2
    """
    import operator
    return operator.sub(*cis.read_data_list(infile, [var1, var2], product))


def ratio_variables(infile, var1, var2, product):
    """
    var1 / var2
    """
    import operator
    return operator.truediv(*cis.read_data_list(infile, [var1, var2], product))


def multiply_sum_by_air_density(infile, variables, product):
    """
    Read multiple variables, sum them together then multiply by air pressure
    (for e.g. calculating aerosol concentrations
    :param infile: CIS filename(s)
    :param variables: CIS variables
    :return:
    """
    from utils import get_stream_file
    summed_var = sum(cis.read_data_list(infile, variables, product))
    vphysc_file = get_stream_file(infile, 'vphysc')
    air_density = cis.read_data(vphysc_file, "rhoam1", product)
    res = summed_var * air_density.data
    res.units = summed_var.units * air_density.units
    return res


def calc_pbl_height(infile, product):
    """
    I want units of m-1, Tau *should* be in units of 1 so I just need to divide by layer thickness
    :param infile:
    :param product:
    :return:
    """
    from .utils import get_stream_file
    altitude = cis.read_data(get_stream_file(infile, 'vphysc'), 'geom1', product) / ACCELERATION_DUE_TO_GRAVITY
    pbl = cis.read_data(get_stream_file(infile, 'vphysc'), 'pbl', product)
    pbl_height = altitude[pbl.data]
    return pbl_height


core = [
    cmor_var('area', 'gboxarea', stream='rad', long_name='horizontal area of grid-box', standard_name='cell_area',
             units=Unit('m^2'), vertical_coord_type='Surface'),
    cmor_var('orog', 'geosp', stream='echam', standard_name='surface_altitude',
             units=Unit('m'), vertical_coord_type='Surface', scaling=1./ACCELERATION_DUE_TO_GRAVITY),
    cmor_var('landf', 'slm', stream='echam', standard_name='land_binary_mask', units=Unit('1'),
             vertical_coord_type='Surface'),

    cmor_var('zmla', calc_pbl_height, stream='vphysc', standard_name='atmosphere_boundary_layer_thickness',
             units=Unit('m'), vertical_coord_type='Surface'),
    cmor_var('layer_thick', 'grheightm1', stream='vphysc', long_name='Layer thickness', standard_name='cell_thickness',
             units=Unit('m'), vertical_coord_type='ModelLevel'),
    cmor_var('zgeo', 'geom1', stream='vphysc', long_name='Geopotential', standard_name='geopotential',
             units=Unit('m')),
    cmor_var('zh', 'geom1', stream='vphysc', long_name='Geopotential Height', standard_name='geopotential_height',
             units=Unit('m'), scaling=1. / ACCELERATION_DUE_TO_GRAVITY),

    cmor_var('ps', 'aps', stream='echam', long_name='Surface air pressure', standard_name='surface_air_pressure',
             units=Unit('Pa'), vertical_coord_type='Surface'),
    cmor_var('psl', 'var151', stream='after', long_name='Sea Level Pressure', standard_name='air_pressure_at_sea_level',
             units=Unit('Pa'), vertical_coord_type='Surface'),
    cmor_var('rho', 'rhoam1', stream='vphysc', long_name='Air density', standard_name='air_density',
             units=Unit('kg m-3')),
    cmor_var('airmass', 'grmassm1', stream='vphysc', standard_name='atmosphere_mass_of_air_per_unit_area',
             units=Unit('kg')),

    cmor_var('ts', 'tslm1', stream='echam', long_name='surface temperature of land',
             standard_name='surface_temperature', units=Unit('K')),
    cmor_var('ta', 'st', stream='echam', standard_name='air_temperature', units=Unit('K'),
             vertical_coord_type='ModelLevel'),

    cmor_var('hus', 'q', stream='echam', long_name='specific humidity', standard_name='specific_humidity',
             units=Unit('1'), vertical_coord_type='ModelLevel'),

    cmor_var('evspsbl', 'evap', stream='echam', long_name='Evaporation', standard_name='water_evaporation_flux',
             units=Unit('kg m-2 s-1'), vertical_coord_type='ModelLevel'),

    # Scale so that it's upward
    cmor_var('hfls', 'ahfl', stream='echam', long_name='Surface Upward Latent Heat Flux',
             standard_name='surface_upward_latent_heat_flux', units=Unit('W m-2'), scaling=-1.0),
    # Scale so that it's upward
    cmor_var('hfss', 'ahfs', stream='echam', long_name='Surface Upward Sensible Heat Flux',
             standard_name='surface_upward_sensible_heat_flux', units=Unit('W m-2'), scaling=-1.0),

    cmor_var('hur', 'relhum', stream='echam',
             long_name='Relative Humidity', standard_name='relative_humidity',
             units=Unit('%'), vertical_coord_type='ModelLevel'),

    cmor_var('hurs', 'relhum', stream='echam', product="ECHAM_HAM_surface_only",
             long_name='Near-Surface Relative Humidity', standard_name='relative_humidity',
             units=Unit('%'), vertical_coord_type='Surface'),

    cmor_var('huss', 'q', stream='echam', long_name='specific humidity', product="ECHAM_HAM_surface_only",
             standard_name='specific_humidity', units=Unit('1'), vertical_coord_type='Surface'),

    cmor_var('sfcWind', 'wind10', stream='echam', long_name='Near-Surface Wind Speed', vertical_coord_type='Surface'),
    cmor_var('tas', 'temp2', stream='echam', long_name='Near-Surface Air Temperature', standard_name='air_temperature',
             units=Unit('K'), vertical_coord_type='Surface'),
    cmor_var('uas', 'u10', stream='echam', long_name='Eastward Near-Surface Wind', standard_name='eastward_wind',
             vertical_coord_type='Surface'),
    cmor_var('vas', 'v10', stream='echam', long_name='Northward Near-Surface Wind', standard_name='northward_wind',
             vertical_coord_type='Surface'),
]

# This uses the netCDF_Gridded product since the vertical coordinates are straight pressure levels
pdrmip_stratified_fields = [
    cmor_var('hur', 'relhum', stream='after', long_name='Relative Humidity',
             standard_name='relative_humidity', units=Unit('%'), product='multi_netcdf'),
    cmor_var('ta', 'st', stream='after', long_name='Air Temperature',
             standard_name='air_temperature', units=Unit('K'), product='multi_netcdf'),
    cmor_var('wap', 'var135', stream='after', long_name='omega (=dp/dt)',
             standard_name='lagrangian_tendency_of_air_pressure', units=Unit('Pa s-1'), product='multi_netcdf'),

    cmor_var('zg', 'var156', stream='after', long_name='Geopotential Height', standard_name='geopotential_height',
             units=Unit('m'), product='multi_netcdf'),
    cmor_var('ua', 'var131', stream='after', long_name='Eastward Wind', standard_name='eastward_wind',
             units=Unit('m s-1'), product='multi_netcdf'),
    cmor_var('va', 'var132', stream='after', long_name='Northward Wind', standard_name='northward_wind',
             units=Unit('m s-1'), product='multi_netcdf'),
]

forcing = [
    cmor_var('R_CLEAN_SKY_SW_TOA', calc_clean_sky_toa_sw_cloud_forcing, stream='echam',
             long_name='"Clean Sky" shortwave cloud radiative forcing', units=Unit('W m-2'), 
             vertical_coord_type='TOA')
]

cloud = [
    cmor_var('ci', 'CONV_TIME', stream='conv', long_name='Fraction of Time Convection Occurs', units=Unit('1'),
             standard_name='convection_time_fraction', vertical_coord_type='ModelLevel'),
    cmor_var('pr', partial(sum_variables, variables=["aprl", "aprc"]), stream='echam',
             long_name='Precipitation',
             standard_name='precipitation_flux', units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('prc', "aprc", stream='echam', long_name='Convective Precipitation',
             standard_name='convective_precipitation_flux', units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('prsn', "aprs", stream='echam', long_name='Snowfall Flux',
             standard_name='snowfall_flux', units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('prw', "qvi", stream='echam', long_name='Water Vapor Path',
             standard_name='atmosphere_water_vapor_content', units=Unit('kg m-2'), vertical_coord_type='Surface'),
    cmor_var('prl', "aprl", stream='echam', long_name='Large-scale Precipitation',
             standard_name='large_scale_precipitation_flux', units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),

    cmor_var('clt', 'aclcov', stream='echam', long_name='Total Cloud Fraction', standard_name='cloud_area_fraction',
             units=Unit('%'), vertical_coord_type='Column'),
    cmor_var('cl', 'aclcac', stream='echam', long_name='Cloud Area Fraction', units=Unit('%'),
             standard_name='cloud_area_fraction_in_atmosphere_layer', vertical_coord_type='ModelLevel'),
    cmor_var('cli', 'xi', stream='echam', long_name='Mass Fraction of Cloud Ice', units=Unit('1'),
             standard_name='mass_fraction_of_cloud_ice_in_air', vertical_coord_type='ModelLevel'),
    cmor_var('clw', 'xl', stream='echam', long_name='Mass Fraction of Cloud Liquid Water', units=Unit('1'),
             standard_name='mass_fraction_of_cloud_liquid_water_in_air', vertical_coord_type='ModelLevel'),
    cmor_var('lwp', 'xlvi', stream='echam', standard_name='atmosphere_mass_content_of_cloud_liquid_water',
             vertical_coord_type='Column'),
    cmor_var('cdnc_cld_top', 'CDNC_CT', stream='activ', long_name='cloud droplet number concentration near cloud top',
             units=Unit('cm-3'), vertical_coord_type='CloudTop'),
    cmor_var('r_e', 'REFFL_CT', stream='activ', long_name='cloud top  effective radius, liquid', units=Unit('m'),
             vertical_coord_type='CloudTop')
]

pdrmip_daily = [
    cmor_var('tasmin', 't2min', stream='hifreq', long_name='Daily Minimum Near-Surface Air Temperature', standard_name='air_temperature',
             units=Unit('K'), vertical_coord_type='Surface'),
    cmor_var('tasmax', 't2max', stream='hifreq', long_name='Daily Maximum Near-Surface Air Temperature', standard_name='air_temperature',
             units=Unit('K'), vertical_coord_type='Surface'),
    cmor_var('pr', partial(sum_variables, variables=["aprl", "aprc"]), stream='hifreq', long_name='precipitation',
             standard_name='precipitation_flux', units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('sfcWind', 'wind10', stream='hifreq', long_name='Daily-Mean Near-Surface Wind Speed',
             vertical_coord_type='Surface', standard_name='wind_speed'),
    # Second table
    cmor_var('rhs', 'relhum', stream='hifreq', vertical_coord_type='Surface', product="ECHAM_HAM_surface_only"),
    cmor_var('prc', 'aprc', stream='hifreq', vertical_coord_type='Surface', standard_name='convective_precipitation_flux'),
    cmor_var('sfcWindmax', 'wimax', stream='hifreq', long_name='Daily Maximum Near-Surface Wind Speed',
             vertical_coord_type='Surface'),
]

pdrmip_fixed_daily = [
    cmor_var('tasmin', 't2min', stream='fixed_daily', long_name='Daily Minimum Near-Surface Air Temperature', standard_name='air_temperature',
             units=Unit('K'), vertical_coord_type='Surface', product='multi_netcdf'),
    cmor_var('tasmax', 't2max', stream='fixed_daily', long_name='Daily Maximum Near-Surface Air Temperature', standard_name='air_temperature',
             units=Unit('K'), vertical_coord_type='Surface', product='multi_netcdf'),
    cmor_var('pr', partial(sum_variables, variables=["aprl", "aprc"]), stream='fixed_daily', long_name='precpitation',
             standard_name='precipitation_flux', units=Unit('kg m-2 s-1'), vertical_coord_type='Surface', product='multi_netcdf'),
    # This hasn't been fixed, it doesn't appear to need scaling by number of days but I don't know why...
    #cmor_var('sfcWind', 'wind10', stream='fixed_daily', long_name='Daily-Mean Near-Surface Wind Speed',
    #         vertical_coord_type='Surface', standard_name='wind_speed', product='multi_netcdf'),
    cmor_var('prc', 'aprc', stream='fixed_daily', vertical_coord_type='Surface', standard_name='convective_precipitation_flux', product='multi_netcdf'),
    cmor_var('sfcWindmax', 'wimax', stream='fixed_daily', long_name='Daily Maximum Near-Surface Wind Speed',
             vertical_coord_type='Surface', product='multi_netcdf'),
]

pdrmip_fixed_monthly = [
    cmor_var('pr', partial(sum_variables, variables=["aprl", "aprc"]), stream='fixed_monthly', long_name='precipitation',
             standard_name='precipitation_flux', units=Unit('kg m-2 s-1'), vertical_coord_type='Surface', product="multi_netcdf"),
    cmor_var('prc', 'aprc', stream='fixed_monthly', vertical_coord_type='Surface',
             standard_name='convective_precipitation_flux', product="multi_netcdf"),
    # This hasn't been fixed, it doesn't appear to need scaling by number of days but I don't know why...
    #cmor_var('sfcWind', 'wind10', stream='fixed_monthly', long_name='Near-Surface Wind Speed', vertical_coord_type='Surface', product="multi_netcdf"),
    #         standard_name='snowfall_flux', units=Unit('kg m-2 s-1'), vertical_coord_type='Surface', product='multi_netcdf'),
    cmor_var('prsn', "aprs", stream='fixed_monthly', long_name='Snowfall Flux',
             standard_name='snowfall_flux', units=Unit('kg m-2 s-1'), vertical_coord_type='Surface', product="multi_netcdf"),
]

aer_rad = [
    cmor_var('od550aer', 'TAU_2D_550nm', stream='rad', long_name='AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column'),
    cmor_var('od440aer', calc_od440aer, stream='rad', long_name='AOD@440nm',
             standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column'),
    cmor_var('od870aer', calc_od870aer, stream='rad', long_name='AOD@870nm',
             standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column'),
    cmor_var('od550lt1aer', partial(sum_variables, variables=['TAU_MODE_K?_550nm', 'TAU_MODE_A?_550nm']), stream='rad',
             long_name='Fine mode AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_pm1_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column',
             comment='Ill-defined: total AOT from Aitken and accumulation modes.'),
    cmor_var('abs550aer', 'ABS_2D_550nm', stream='rad', long_name='Absorption AOD@550nm',
             standard_name='atmosphere_absorption_optical_thickness_due_to_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column'),
    cmor_var('od550aerh2o', 'TAU_COMP_WAT_550nm', stream='rad', long_name='Aerosol Water AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_water_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column',
             comment='Ill-defined: model allows 4 modes of varying mixtures.'),
    cmor_var('abs550bc', 'ABS_COMP_BC_550nm', stream='rad', long_name='Absorption AOD@550nm due to BC',
             standard_name='atmosphere_absorption_optical_thickness_due_to_black_carbon_ambient_aerosol',
             units=Unit('1'), vertical_coord_type='Column',
             comment='Ill-defined: model allows 4 modes of varying mixtures.'),
    cmor_var('od550so4', 'TAU_COMP_SO4_550nm', stream='rad', long_name='Sulfate AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_sulfate_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column',
             comment='Ill-defined: model allows 4 modes of varying mixtures.'),
    cmor_var('od550bc', 'TAU_COMP_BC_550nm', stream='rad', long_name='Black carbon AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_black_carbon_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column',
             comment='Ill-defined: model allows 4 modes of varying mixtures.'),
    cmor_var('od550oa', 'TAU_COMP_OC_550nm', stream='rad', long_name='POM AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_particulate_organic_matter_ambient_aerosol',
             units=Unit('1'), vertical_coord_type='Column',
             comment='Ill-defined: model allows 4 modes of varying mixtures.'),
    cmor_var('od550ss', 'TAU_COMP_SS_550nm', stream='rad', long_name='Sea Salt AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_seasalt_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column',
             comment='Ill-defined: model allows 4 modes of varying mixtures.'),
    cmor_var('od550dust', 'TAU_COMP_DU_550nm', stream='rad', long_name='Dust AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_dust_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column',
             comment='Ill-defined: model allows 4 modes of varying mixtures.'),
    cmor_var('ext550aer', calc_ext550aer, stream='rad', long_name='3D Aerosol Extinction @550nm',
             units=Unit('m-1'), standard_name='volume_extinction_coefficient_in_air_due_to_ambient_aerosol_particles')
]

rad = [
    # LW
    # 'LW surface down-welling' = 'net surface thermal radiation' - 'surface thermal radiation upward'
    # Note that 'net surface thermal radiation' is net downwards, and 'surface thermal radiation upward' is actually
    #  +ve down, but these signs cancel in the below operation to give the right sign overall.
    # Positive: down
    cmor_var('rlds', partial(diff_variables, var1='trads', var2='tradsu'), stream='echam',
             long_name='Surface Downwelling Longwave Radiation',
             standard_name='surface_downwelling_longwave_flux_in_air', units=Unit('W m-2'),
             vertical_coord_type='Surface'),
    # 'LW surface down-welling (clear sky)' = 'net surface thermal radiation (clear-sky)' - 'surface thermal radiation upward'
    # See above about the signs here
    # Positive: down
    cmor_var('rldscs', partial(diff_variables, var1='trafs', var2='tradsu'), stream='echam',
             long_name='Surface Downwelling Clear-Sky Longwave Radiation',
             standard_name='surface_downwelling_longwave_flux_in_air_assuming_clear_sky', units=Unit('W m-2'),
             vertical_coord_type='Surface'),
    # Positive: up
    # Scale this since the 'surface thermal radiation upward' is +ve down...!
    cmor_var('rlus', 'tradsu', stream='echam', long_name='LW surface up-welling ',
             standard_name='surface_upwelling_longwave_flux_in_air', units=Unit('W m-2'),
             vertical_coord_type='Surface', scaling=-1.),
    # Positive: up
    # Scale this since the 'top thermal radiation (OLR)' is +ve down...!
    cmor_var('rlut', 'trad0', stream='echam', long_name='TOA Outgoing Longwave Radiation',
             standard_name='toa_outgoing_longwave_flux', units=Unit('W m-2'),
             vertical_coord_type='TOA', scaling=-1.),
    # Positive: up
    # Scale this since the 'top thermal radiation (OLR)' is +ve down...!
    cmor_var('rlutcs', 'traf0', stream='echam', long_name='TOA Outgoing Clear-Sky Longwave Radiation',
             standard_name='toa_outgoing_longwave_flux_assuming_clear_sky', units=Unit('W m-2'),
             vertical_coord_type='TOA', scaling=-1.),
    # SW
    # Note that 'net surface solar radiation' is net downwards, and 'surface solar radiation upward' is actually
    #  +ve down, but these signs cancel in the below operation to give the right sign overall.
    #  Positive: down
    cmor_var('rsds', partial(diff_variables, var1='srads', var2='sradsu'), stream='echam',
             long_name='Surface Downwelling Shortwave Radiation',
             standard_name='surface_downwelling_shortwave_flux_in_air', units=Unit('W m-2'),
             vertical_coord_type='Surface'),

    # Because the SW radiation upward doesn't have a clear sky analogue, and because the upward solar radiation
    #  will be different in clear/all sky conditions we can't decompose the net clear sky into upward and downward
    #  fluxes directly. Physically the upward SW at the surface is just a reflection, wheras the LW we can use the all
    #  sky.
    # # 'SW surface down-welling (clear sky)' = 'net surface solar radiation (clear sky)' - 'surface solar radiation upward'
    # cmor_var('rsdscs', partial(diff_variables, var1='srafs', var2='sradsu'), stream='',
    #          long_name='Surface Downwelling Clear-Sky Shortwave Radiation',
    #          standard_name='surface_downwelling_shortwave_flux_in_air_assuming_clear_sky',
    #          units=Unit('W m-2'), vertical_coord_type='Surface'),
    # # SW upwelling solar flux (clear sky) = 'net surface solar radiation (clear sky)' - 'SW surface down-welling (clear sky)'
    # #                           rsuscs    =  srafs - (srafs - sradsu)
    # #                      ->   rsuscs    =  sradsu    ...
    # cmor_var('rsuscs', 'traf0', stream='echam', long_name='Surface Upwelling Clear-Sky Shortwave Radiation',
    #          standard_name='surface_upwelling_shortwave_flux_in_air_assuming_clear_sky', units=Unit('W m-2'),
    #          vertical_coord_type='Surface'),

    # Positive: up
    cmor_var('rsnscs', 'srafs', stream='echam', long_name='net surface solar radiation (clear sky)',
             standard_name='surface_net_downward_shortwave_flux_assuming_clear_sky', units=Unit('W m-2'),
             vertical_coord_type='Surface', scaling=-1.),
    # Positive: down
    cmor_var('rsdt', 'srad0d', stream='echam', long_name='TOA Incident Shortwave Radiation',
             standard_name='toa_incoming_shortwave_flux', units=Unit('W m-2'), vertical_coord_type='TOA'),
    # Positive: up
    cmor_var('rsus', 'sradsu', stream='echam', long_name='Surface Upwelling Shortwave Radiation ',
             standard_name='surface_upwelling_shortwave_flux_in_air', units=Unit('W m-2'),
             vertical_coord_type='Surface', scaling=-1.),
    # Positive: up
    cmor_var('rsut', 'srad0u', stream='echam', long_name='SW upwelling solar flux',
             standard_name='toa_outgoing_shortwave_flux', units=Unit('W m-2'),
             vertical_coord_type='TOA', scaling=-1.),

    # 'SW upwelling solar flux in clear sky regions' = 'top incoming solar radiation' - 'net downward top solar radiation (clear sky)'
    # Positive: up
    cmor_var('rsutcs', partial(diff_variables, var1='srad0d', var2='sraf0'), stream='echam',
             long_name='SW upwelling solar flux in clear sky regions',
             standard_name='toa_outgoing_shortwave_flux_assuming_clear_sky', units=Unit('W m-2'),
             vertical_coord_type='TOA'),
]

double_rad = [
    # ALL scaled by -1. to turn net Downward fluxes into net Upward fluxes
    cmor_var('rsntcs_irf', get_clean_clear_sky_toa_net_sw, stream='echam', units=Unit('W m-2'),
             long_name='net TOA solar radiation (clear sky, no aerosol)', vertical_coord_type='TOA',
             scaling=-1.),

    cmor_var('rsnscs_irf', get_clean_clear_sky_surface_net_sw, stream='echam', units=Unit('W m-2'),
             long_name='net surface solar radiation (clear sky, no aerosol)', vertical_coord_type='Surface',
             scaling=-1.),

    cmor_var('rsnt_irf', get_clean_all_sky_toa_net_sw, stream='echam', units=Unit('W m-2'),
             long_name='net TOA solar radiation (no aerosol)', vertical_coord_type='TOA',
             scaling=-1.),

    cmor_var('rsns_irf', get_clean_all_sky_surface_net_sw, stream='echam', units=Unit('W m-2'),
             long_name='net surface solar radiation (no aerosol)', vertical_coord_type='Surface',
             scaling=-1.),

    cmor_var('rlntcs_irf', get_clean_clear_sky_toa_net_lw, stream='echam', units=Unit('W m-2'),
             long_name='net TOA longwave radiation (clear sky, no aerosol)', vertical_coord_type='TOA',
             scaling=-1.),

    cmor_var('rlnscs_irf', get_clean_clear_sky_surface_net_lw, stream='echam', units=Unit('W m-2'),
             long_name='net longwave surface radiation (clear sky, no aerosol)', vertical_coord_type='Surface',
             scaling=-1.),

    cmor_var('rlnt_irf', get_clean_all_sky_toa_net_lw, stream='echam', units=Unit('W m-2'),
             long_name='net TOA longwave radiation (no aerosol)', vertical_coord_type='TOA',
             scaling=-1.),

    cmor_var('rlns_irf', get_clean_all_sky_surface_net_lw, stream='echam', units=Unit('W m-2'),
             long_name='net surface longwave radiation (no aerosol)', vertical_coord_type='Surface',
             scaling=-1.)
]


aerosol = [
    cmor_var('emioa', 'emi_OC', stream='emi', long_name='total emission of POM',
             standard_name='tendency_of_atmosphere_mass_content_of_particulate_organic_matter_dry_aerosol_due_to_net_chemical_production_and_emission',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('emibc', 'emi_BC', stream='emi', long_name='total emission of BC',
             standard_name='tendency_of_atmosphere_mass_content_of_black_carbon_dry_aerosol_due_to_emission',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('emiso2', 'emi_SO2', stream='emi', long_name='total emission of SO2',
             standard_name='tendency_of_atmosphere_mass_content_of_sulfur_dioxide_due_to_emission',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('emiso4', 'emi_SO4', stream='emi', long_name='total direct emission of SO4',
             standard_name='tendency_of_atmosphere_mass_content_of_sulfate_dry_aerosol_due_to_emission',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('emidms', 'emi_DMS', stream='emi', long_name='total emission of DMS',
             standard_name='tendency_of_atmosphere_mass_content_of_dimethyl_sulfide_due_to_emission',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('emiss', 'emi_SS', stream='emi', long_name='total emission of seasalt',
             standard_name='tendency_of_atmosphere_mass_content_of_seasalt_dry_aerosol_due_to_emission',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('emidust', 'emi_DU', stream='emi', long_name='total emission of dust',
             standard_name='tendency_of_atmosphere_mass_content_of_dust_dry_aerosol_due_to_emission',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('emibb', partial(sum_variables, variables='emi_OC_?fire'), stream='emi',
             long_name='total emission of Biomass Burning Aerosol',
             standard_name='tendency_of_atmosphere_mass_content_of_particulate_organic_matter_dry_aerosol_due_to_emission',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),

    cmor_var('wetoa', partial(sum_variables, variables='wdep_OC_??'), stream='wetdep', long_name='wet deposition of POM',
             standard_name='tendency_of_atmosphere_mass_content_of_particulate_organic_matter_dry_aerosol_due_to_wet_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('wetbc', partial(sum_variables, variables='wdep_BC_??'), stream='wetdep', long_name='wet deposition of BC',
             standard_name='tendency_of_atmosphere_mass_content_of_black_carbon_dry_aerosol_due_to_wet_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('wetso4', partial(sum_variables, variables='wdep_SO4_??'), stream='wetdep',
             long_name='wet deposition of SO4',
             standard_name='tendency_of_atmosphere_mass_content_of_sulfate_dry_aerosol_due_to_wet_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('wetso2', 'wdep_SO2', stream='wetdep', long_name='wet deposition of SO2',
             standard_name='tendency_of_atmosphere_mass_content_of_sulfur_dioxide_due_to_wet_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('wetss', partial(sum_variables, variables='wdep_SS_??'), stream='wetdep',
             long_name='wet deposition of seasalt',
             standard_name='tendency_of_atmosphere_mass_content_of_seasalt_dry_aerosol_due_to_wet_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('wetdust', partial(sum_variables, variables='wdep_DU_??'), stream='wetdep',
             long_name='wet deposition of dust',
             standard_name='tendency_of_atmosphere_mass_content_of_dust_dry_aerosol_due_to_wet_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),

    cmor_var('loadoa', 'burden_OC', stream='burden', long_name='Load of POM',
             standard_name='atmosphere_mass_content_of_particulate_organic_matter_dry_aerosol', units=Unit('kg m-2'),
             vertical_coord_type='Column'),
    cmor_var('loadbc', 'burden_BC', stream='burden', long_name='Load of BC',
             standard_name='atmosphere_mass_content_of_black_carbon_dry_aerosol', units=Unit('kg m-2'),
             vertical_coord_type='Column'),
    cmor_var('loadso4', 'burden_SO4', stream='burden', long_name='Load of SO4',
             standard_name='atmosphere_mass_content_of_sulfate_dry_aerosol', units=Unit('kg m-2'),
             vertical_coord_type='Column'),
    cmor_var('loaddust', 'burden_DU', stream='burden', long_name='Load of DUST',
             standard_name='atmosphere_mass_content_of_dust_dry_aerosol', units=Unit('kg m-2'),
             vertical_coord_type='Column'),
    cmor_var('loadss', 'burden_SS', stream='burden', long_name='Load of SS',
             standard_name='atmosphere_mass_content_of_seasalt_dry_aerosol', units=Unit('kg m-2'),
             vertical_coord_type='Column'),

    cmor_var('dryso2', 'ddep_SO2', stream='drydep', long_name='dry deposition of SO2',
             standard_name='tendency_of_atmosphere_mass_content_of_sulfur_dioxide_due_to_dry_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('dryso4', partial(sum_variables, variables='ddep_SO4_??'), stream='drydep',
             long_name='dry deposition of SO4',
             standard_name='tendency_of_atmosphere_mass_content_of_sulfate_due_to_dry_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('dryss', partial(sum_variables, variables='ddep_SS_??'), stream='drydep',
             long_name='dry deposition of seasalt',
             standard_name='tendency_of_atmosphere_mass_content_of_seasalt_dry_aerosol_due_to_dry_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('drydust', partial(sum_variables, variables='ddep_DU_??'), stream='drydep',
             long_name='dry deposition of dust',
             standard_name='tendency_of_atmosphere_mass_content_of_dust_dry_aerosol_due_to_dry_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('dryoa', partial(sum_variables, variables='ddep_OC_??'), stream='drydep',
             long_name='dry deposition of POM',
             standard_name='tendency_of_atmosphere_mass_content_of_particulate_organic_matter_dry_aerosol_due_to_dry_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),
    cmor_var('drybc', partial(sum_variables, variables='ddep_BC_??'), stream='drydep', long_name='dry deposition of BC',
             standard_name='tendency_of_atmosphere_mass_content_of_black_carbon_dry_aerosol_due_to_dry_deposition',
             units=Unit('kg m-2 s-1'), vertical_coord_type='Surface'),

    cmor_var('mmraerh2o', partial(sum_variables, variables='WAT_??'), stream='tracer', long_name='mmr of aerosol water',
             standard_name='mass_fraction_of_water_in_ambient_aerosol_in_air', units=Unit('1'),
             vertical_coord_type='ModelLevel'),
    cmor_var('mmroa', partial(sum_variables, variables='OC_??'), stream='tracer', long_name='POM',
             standard_name='mass_fraction_of_particulate_organic_matter_dry_aerosol_in_air', units=Unit('1'),
             vertical_coord_type='ModelLevel'),
    cmor_var('mmrbc', partial(sum_variables, variables='BC_??'), stream='tracer', long_name='BC',
             standard_name='mass_fraction_ of_black_carbon_dry_aerosol_in_air', units=Unit('1'),
             vertical_coord_type='ModelLevel'),
    cmor_var('mmrso2', 'SO2', stream='tracer', long_name='SO2',
             standard_name='mass_fraction_of_sulfur_dioxide_in_air', units=Unit('1'),
             vertical_coord_type='ModelLevel'),
    cmor_var('mmrso4', partial(sum_variables, variables='SO4_??'), stream='tracer', long_name='SO4',
             standard_name='mass_fraction_of_sulfate_dry_aerosol_in_air', units=Unit('1'),
             vertical_coord_type='ModelLevel'),
    cmor_var('mmrss', partial(sum_variables, variables='SS_??'), stream='tracer', long_name='Sea Salt',
             standard_name='mass_fraction_ of_seasalt_dry_aerosol_in_air', units=Unit('1'),
             vertical_coord_type='ModelLevel'),
    cmor_var('mmrdu', partial(sum_variables, variables='DU_??'), stream='tracer', long_name='Dust',
             standard_name='mass_fraction_ of_dust_dry_aerosol_in_air', units=Unit('1'),
             vertical_coord_type='ModelLevel'),

    cmor_var('vmrso2', 'SO2', stream='tracer', long_name='SO2', standard_name='mole_fraction_of_sulfur_dioxide_in_air',
             units=Unit('1'), vertical_coord_type='ModelLevel', scaling=(28.97 / 64.066)),
    cmor_var('vmrso4', 'H2SO4', stream='tracer', long_name='SO4', standard_name='mole_fraction_of_sulfate_in_air',
             units=Unit('1'), vertical_coord_type='ModelLevel', scaling=(28.97 / 96.063)),
    cmor_var('vmrdms', 'DMS', stream='tracer', long_name='DMS',
             standard_name='mole_fraction_of_dimethyl_sulfide_in_air', units=Unit('1'),
             vertical_coord_type='ModelLevel', scaling=(28.97 / 62.134)),

    cmor_var('cheaqpso4', 'D_PROD_SO4_LIQ', stream='ham', long_name='aqu phase production so4',
             standard_name='tendency_of_atmosphere_mass_content_of_sulfate_dry_aerosol_due_to_net_chemical_production_and_emission',
             units=Unit('kg m-2 s-1'), vertical_coord_type='ModelLevel'),

    cmor_var('concbc', partial(multiply_sum_by_air_density, variables='BC_??'), stream='tracer',
             long_name='Concentration of Black Carbon Aerosol', units=Unit('kg m-3'), vertical_coord_type='ModelLevel'),
    cmor_var('concso4', partial(multiply_sum_by_air_density, variables='SO4_??'), stream='tracer',
             long_name='Concentration of SO4', units=Unit('kg m-3'), vertical_coord_type='ModelLevel'),

    cmor_var('sconcoa', partial(multiply_sum_by_air_density, variables='OC_??'), stream='tracer',
             long_name='Surface concentration POM', product='ECHAM_HAM_surface_only',
             standard_name='mass_concentration_of_particulate_organic_matter_dry_aerosol_in_air', units=Unit('kg m-3'),
             vertical_coord_type='Surface'),
    cmor_var('sconcbc', partial(multiply_sum_by_air_density, variables='BC_??'), stream='tracer',
             long_name='Surface concentration BC', product='ECHAM_HAM_surface_only',
             standard_name='mass_concentration_of_black_carbon_dry_aerosol_in_air', units=Unit('kg m-3'),
             vertical_coord_type='Surface'),
    cmor_var('sconcso4', partial(multiply_sum_by_air_density, variables='SO4_??'), stream='tracer',
             long_name='Surface concentration SO4',
             standard_name='mass_concentration_of_sulfate_dry_aerosol_in_air', units=Unit('kg m-3'),
             vertical_coord_type='Surface'),
    cmor_var('sconcdust', partial(multiply_sum_by_air_density, variables='DU_??'), stream='tracer',
             long_name='Surface concentration DUST',
             standard_name='mass_concentration_of_dust_dry_aerosol_in_air', units=Unit('kg m-3'),
             vertical_coord_type='Surface'),
    cmor_var('sconcss', partial(multiply_sum_by_air_density, variables='SS_??'), stream='tracer',
             long_name='Surface concentration SS',
             standard_name='mass_concentration_of_seasalt_dry_aerosol_in_air', units=Unit('kg m-3'),
             vertical_coord_type='Surface'),

    cmor_var('conccnmodeNS', partial(sum_variables, variables='NUM_NS'), stream='tracer',
             long_name='number concentration of mode NS', units=Unit('m-3'), vertical_coord_type='ModelLevel'),
    cmor_var('conccnmodeKS', partial(sum_variables, variables='NUM_KS'), stream='tracer',
             long_name='number concentration of mode KS', units=Unit('m-3'), vertical_coord_type='ModelLevel'),
    cmor_var('conccnmodeAS', partial(sum_variables, variables='NUM_AS'), stream='tracer',
             long_name='number concentration of mode AS', units=Unit('m-3'), vertical_coord_type='ModelLevel'),
    cmor_var('conccnmodeCS', partial(sum_variables, variables='NUM_CS'), stream='tracer',
             long_name='number concentration of mode CS', units=Unit('m-3'), vertical_coord_type='ModelLevel'),
    cmor_var('conccnmodeKI', partial(sum_variables, variables='NUM_KI'), stream='tracer',
             long_name='number concentration of mode KI', units=Unit('m-3'), vertical_coord_type='ModelLevel'),
    cmor_var('conccnmodeAI', partial(sum_variables, variables='NUM_AI'), stream='tracer',
             long_name='number concentration of mode AI', units=Unit('m-3'), vertical_coord_type='ModelLevel'),
    cmor_var('conccnmodeCI', partial(sum_variables, variables='NUM_CI'), stream='tracer',
             long_name='number concentration of mode CI', units=Unit('m-3'), vertical_coord_type='ModelLevel'),

    cmor_var('ddrymodeNS', 'rdry_NS', stream='ham',
             long_name='dry diameter of mode NS', units=Unit('m'), vertical_coord_type='ModelLevel'),
    cmor_var('ddrymodeKS', 'rdry_KS', stream='ham',
             long_name='dry diameter of mode KS', units=Unit('m'), vertical_coord_type='ModelLevel'),
    cmor_var('ddrymodeAS', 'rdry_AS', stream='ham',
             long_name='dry diameter of mode AS', units=Unit('m'), vertical_coord_type='ModelLevel'),
    cmor_var('ddrymodeCS', 'rdry_CS', stream='ham',
             long_name='dry diameter of mode CS', units=Unit('m'), vertical_coord_type='ModelLevel'),
    cmor_var('ddrymodeKI', 'rdry_KI', stream='ham',
             long_name='dry diameter of mode KI', units=Unit('m'), vertical_coord_type='ModelLevel'),
    cmor_var('ddrymodeAI', 'rdry_AI', stream='ham',
             long_name='dry diameter of mode AI', units=Unit('m'),  vertical_coord_type='ModelLevel'),
    cmor_var('ddrymodeCI', 'rdry_CI', stream='ham',
             long_name='dry diameter of mode CI', units=Unit('m'), vertical_coord_type='ModelLevel')

]

all_vars = core + cloud + rad + aerosol

# ----------- Nick's RemSens experiment -----------
remsens_3hrly = [
    cmor_var('deltaz3d', 'grheightm1', stream='hifreq',
             long_name='Layer thickness', standard_name='cell_thickness',
             units=Unit('m'), vertical_coord_type='ModelLevel'),
    cmor_var('humidity3d', 'relhum', stream='hifreq'),

    # Needs ["TAU_2D_550nm", "ANG_550nm_865nm"]
    cmor_var('od440aer', calc_od440aer, stream='hifreq', long_name='AOD@440nm',
             standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column'),
    cmor_var('od550aer', 'TAU_2D_550nm', stream='hifreq', long_name='AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column'),
    cmor_var('od550dryaer', 'TAU_DRY_2D_550nm', stream='hifreq', long_name='Dry AOD@550nm', units=Unit('1'),
             vertical_coord_type='Column'),
    # Needs ["TAU_2D_865nm", "ANG_550nm_865nm"]
    cmor_var('od870aer', calc_od870aer, stream='hifreq', long_name='AOD@870nm',
             standard_name='atmosphere_optical_thickness_due_to_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column'),
    cmor_var('abs550aer', 'ABS_2D_550nm', stream='hifreq', long_name='Absorption AOD@550nm',
             standard_name='atmosphere_absorption_optical_thickness_due_to_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column'),
    cmor_var('od550lt1aer', calc_od550lt1aer, stream='hifreq', long_name='Fine mode AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_pm1_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column', comment='Ill-defined: total AOT from Aitken and accumulation modes.'),
    cmor_var('od550dust', 'TAU_COMP_DU_550nm', stream='hifreq', long_name='Dust AOD@550nm',
             standard_name='atmosphere_optical_thickness_due_to_dust_ambient_aerosol', units=Unit('1'),
             vertical_coord_type='Column',
             comment='Ill-defined: model allows 4 modes of varying mixtures.'),
    cmor_var('od550aer3d', partial(sum_variables, variables='TAU_MODE_??_550nm'), stream='hifreq',
             long_name='Layer Aerosol Optical Thickness @550nm',
             units=Unit('1'), vertical_coord_type='ModelLevel'),
]

# ---------------- PDRMIP specific variables -------------------
pdrmip_core = select_vars(core, ["evspsbl", "hfls", "hfss", "hurs", "hus", "huss", "ps", "psl", "sfcWind", 'tas', 'uas',
                                 "vas"])
pdrmip_aer_rad = select_vars(aer_rad, ['od550aer', 'abs550aer'])
pdrmip_rad = select_vars(rad, ["rlds", "rldscs", "rlus", "rlut", "rlutcs", "rsds", "rsutcs", "rsdt", "rsus", "rsut",
                               # "rsuscs", "rsdscs", IGNORE the up and down surface SW cs since we can't calculate it
                               "rsnscs"])  # Output the net instead
pdrmip_cloud = select_vars(cloud, ["pr", "prc", "prsn", "prw", "cl", "cli", "clw", "clt", "ci"])
pdrmip_aerosol = select_vars(aerosol, ["emibc", "emiso2", "emiso4", "loadbc", "loadso4", "concbc", "concso4"])
pdrmip = pdrmip_daily + pdrmip_aerosol + pdrmip_cloud + pdrmip_stratified_fields + pdrmip_core + pdrmip_rad + \
         pdrmip_aer_rad + double_rad

#  ------------ Aerocom CTRL (2D monthly) fields  -------------
#  Variables-Parameter:Speciation
# EMI-Emissions: BC, OA, SO2, DMS, NOx, VOC, DUST, SS
aerocom_emi = select_vars(aerosol, ['emibc', 'emioa', 'emiso2', 'emidms', 'emidust', 'emiss'])  # No VOC or NOx
# (column integrated, if emission at altitude, eg SOA is accumulated in OA emissions)
# LOAD-Column Loads: BC, OA, SO4, NO3, DUST, SS
aerocom_load = select_vars(aerosol, ['loadbc', 'loadoa', 'loadso4', 'loaddust', 'loadss'])  # No NO3
# SCONC-Surface concentrations: PM10, PM25, BC, OA, SO4, NO3, DUST, SS
aerocom_sconc = select_vars(aerosol, ['sconcoa', 'sconcbc', 'sconcso4', 'sconcdust', 'sconcss'])  # No PM10, PM25, NO3
# DEP-Total Deposition: BC, OA, SO4, NO3, DUST, SS
# This is a pain, I'll do as needed. Do weet and dry separately for now
aerocom_wet = select_vars(aerosol, ['drybc', 'dryoa', 'dryso4', 'drydust', 'dryss'])  # No NO3
aerocom_dry = select_vars(aerosol, ['wetbc', 'wetoa', 'wetso4', 'wetdust', 'wetss'])  # No NO3

# OD550-Aerosol optical depth @550nm: AER, fine mode AER, coarse mode AOD, (tier 2: BC, OA, SO4, NO3, DUST, SS)
aerocom_aer_rad = select_vars(aer_rad, ['od550aer', 'od550lt1aer', # No coarse mode AOD
                                        'od550bc', 'od550oa', 'od550so4', 'od550dust', 'od550ss']) # No NO3 AOD
# Total AOD effective in radiative forcing code (OD550AER) and clearsky AOD (od550csaer)
#   This is the same in ECHAM-HAM
# SWTOA-Top of Atm Fluxes clear-sky and all-sky : AER, fine mode AER, coarse mode
#                                                 (tier 2: AOD, BC, OA, SO4, NO3, DUST, SS)
aerocom_rad = select_vars(rad, ['rsutcs', 'rsut'])
# LWTOA-Fluxes clear-sky and all-sky : AER, fine mode AER, coarse mode
aerocom_rad += select_vars(rad, ['rlutcs', 'rlut'])
# CCN Number concentration @ 850 hPa
#   Ill-defined - what supersaturation?
# IN Number concentration @ 100 hPa
#   Don't know...
# Total Cloud cover
# Cloud Water Path
# Low level cloud cover
# Precipitation rate
aerocom_cloud = select_vars(cloud, ['clt', 'lwp', 'pr'])  # No Low-level cloud cover

aerocom_ctrl = aerocom_emi + aerocom_load + aerocom_sconc + aerocom_wet + aerocom_dry + aerocom_aer_rad\
               + aerocom_rad + aerocom_cloud

# -------------- AeroCom Holuhraun Experiment ---------------
# Cloud
holuhraun_cloud = select_vars(cloud, [
    # 'ccn_0p3pc_1km',         Not done yet, doable but probably not a priority
    'cdnc_cld_top',
    # 'temp_cld_top'           This could be tricky, need to pull out cloud top heights then get temperatures I guess
    'clt',  # (cld_frac)
    # 'cot',                   Can't seem to find...
    'lwp',
    'r_e',
    'prc',  # (rain_conv)
    'prl',  # (rain_ls)
    'cl',
])

holuhraun_aer_rad = select_vars(aer_rad, [
    'od550aer',
    'od440aer',
    'od870aer',
])

holuhraun_core = select_vars(core, [
    # 'LTS'                    lower troposheric stability: Theta(700hPa) - Theta(SRF), Not sure where to get this
    'ps',
    'ts',
    'layer_thick',
    # 'p',                     Can I just get this Iris cube??
    'pho',
    'hus',
    'hur'
])

holuhraun_aer = select_vars(aerosol, [
    # 'emiso2',  # I don't have emiso2_srf or emiso2_high. This is actually accumulated emissions too...
    'wetso2',
    # 'wetso4',                This isn't per mode in the holuhraun output...
    'cheaqpso4',
    'mmrso2',
    'mmrso4',

])

holuhraun = holuhraun_cloud + holuhraun_aer_rad + holuhraun_core + holuhraun_aer + rad + forcing

# -------------- AeroCom Trajectory Experiment ---------------


trajectory_2d = select_vars(core, [
    'orog',
    'landf',
    'area',

    'ps',
    'zmla',
    'u10',
    'v10',
    'tas',
    'hfss',
])

trajectory_cloud = select_vars(cloud, ['prl', 'prc'])

trajectory_3d = select_vars(core, [
    'zgeo',
    'ta',
    # 'ua',  # I'll have to run afterburner if I really need these (I suspect they're in the GRIB files though
    # 'va',
    'hus',
    # 'omega',
    'hur',
    'rho',
    # 'plev',  # This will be in the other files anyway
    'airmass',
    'zh',
    'layer_thick',
])

trajectory_aer = select_vars(aerosol, ['conccnmodeNS', 'conccnmodeKS', 'conccnmodeAS', 'conccnmodeCS',
                                       'conccnmodeKI', 'conccnmodeAI', 'conccnmodeCI',
                                       'mmroa', 'mmrbc', 'mmrso4', 'mmrss', 'mmrdu',
                                       'ddrymodeNS', 'ddrymodeKS', 'ddrymodeAS', 'ddrymodeCS',
                                       'ddrymodeKI', 'ddrymodeAI', 'ddrymodeCI'])

trajectory = trajectory_2d + trajectory_cloud + trajectory_3d + trajectory_aer


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("infile", help="input file(s)")
    parser.add_argument("outbase", help="output filename base: aerocom3_<!ModelName>_<!ExperimentName>")
    parser.add_argument("-t", "--time", help="time period to use for variables that don't have a unique one")
    parser.add_argument("--experiment", help="More detailed experiment information")
    parser.add_argument("--contact", help="Contact details")
    parser.add_argument("--institute", help="Institute details")
    parser.add_argument("-o", "--overwrite", help="Force overwrite of existing file (default False)",
                        action='store_true')
    parser.add_argument("--product", help="The CIS product to use")
    parser.add_argument("--pdrmip_format", help="Use the PDRMIP filename formatting style", action='store_true')

    # Parameter sets
    parser.add_argument('-c', '--core', action='append_const', const=core, dest='params')
    parser.add_argument('--cloud', action='append_const', const=cloud, dest='params')
    parser.add_argument('-r', '--rad', action='append_const', const=rad, dest='params')
    parser.add_argument('-a', '--aerosol', action='append_const', const=aerosol, dest='params')
    parser.add_argument('--pdrmip', action='append_const', const=pdrmip, dest='params')
    parser.add_argument('--aerocom', action='append_const', const=aerocom_ctrl, dest='params',
                        help="Core AeroCom diagnostics")
    parser.add_argument('--holuhraun', action='append_const', const=holuhraun, dest='params',
                        help="AeroCom Holuhraun experiment diagnostics")
    parser.add_argument('--remsens', action='append_const', const=remsens_3hrly, dest='params',
                        help="Diagnostics for Nick's remote sensing experiment")
    parser.add_argument('--fixed_pdrmip', action='append_const', const=pdrmip_fixed_daily+pdrmip_fixed_monthly, dest='params',
                                    help="Fixed PDRMIP diagnostics (due to hifreq output issue)")
    parser.add_argument('--all', action='append_const', const=all_vars, dest='params')
    parser.add_argument('--trajectory', action='append_const', const=trajectory, dest='params')

    freq_grp = parser.add_mutually_exclusive_group()
    freq_grp.add_argument("-m", "--monthly", action="store_true",
                          help="assume monthly output")
    freq_grp.add_argument("-d", "--daily", action="store_true",
                          help="assume daily output")

    args = parser.parse_args()

    # Flatten the parameters and remove any duplicates
    variables = set(itertools.chain.from_iterable(args.params))

    for v in variables:
        print("Processing {}...".format(v.cmor_var_name))
        c = v.load_var(args.infile, args.product)
        print("Global (un-weighted) mean: {}".format(c.data.mean()))
        v.write_var(c, args.time, args.outbase, args.daily, args.monthly, args.experiment, args.contact,
                    args.overwrite, args.pdrmip_format)
        print("..done")
