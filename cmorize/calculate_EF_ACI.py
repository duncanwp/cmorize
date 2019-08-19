"""
(c) Duncan watson-parris 2017
"""
import cis
from utils import get_stream_file

def get_clean_all_sky_toa_net_lw(infile, product=None):
    """
    This is $F_{lw,all}^{\top,0}$ as defined in A.7.1 of ECHAM6_userguide.pdf

    :param infile: An ECHAM-HAM echam filepath
    :param product:
    :return:
    """
    forcing_file = get_stream_file(infile, 'forcing')
    clean_all_sky_toa_sw_forcing = cis.read_data(infile, 'trad0', product) - \
                                   cis.read_data(forcing_file, 'FLW_TOTAL_TOP', product)
    return clean_all_sky_toa_sw_forcing


def get_clean_all_sky_surface_net_lw(infile, product=None):
    """
    This is $F_{lw,all}^{\bot,0}$ as defined in A.7.1 of ECHAM6_userguide.pdf

    :param infile: An ECHAM-HAM echam filepath
    :param product:
    :return:
    """
    forcing_file = get_stream_file(infile, 'forcing')
    clean_all_sky_toa_sw_forcing = cis.read_data(infile, 'trads', product) - \
                                   cis.read_data(forcing_file, 'FLW_TOTAL_SUR', product)
    return clean_all_sky_toa_sw_forcing


def get_clean_clear_sky_toa_net_lw(infile, product=None):
    """
    This is $F_{lw,clear}^{\top,0}$ as defined in A.7.1 of ECHAM6_userguide.pdf

    :param infile: An ECHAM-HAM echam filepath
    :param product:
    :return:
    """
    forcing_file = get_stream_file(infile, 'forcing')
    clean_clear_sky_toa_sw_forcing = cis.read_data(infile, 'traf0', product) - \
                                     cis.read_data(forcing_file, 'FLW_CLEAR_TOP', product)
    return clean_clear_sky_toa_sw_forcing


def get_clean_clear_sky_surface_net_lw(infile, product=None):
    """
    This is $F_{lw,clear}^{\bot,0}$ as defined in A.7.1 of ECHAM6_userguide.pdf

    :param infile: An ECHAM-HAM echam filepath
    :param product:
    :return:
    """
    forcing_file = get_stream_file(infile, 'forcing')
    clean_clear_sky_toa_sw_forcing = cis.read_data(infile, 'trafs', product) - \
                                     cis.read_data(forcing_file, 'FLW_CLEAR_SUR', product)
    return clean_clear_sky_toa_sw_forcing


def get_clean_all_sky_toa_net_sw(infile, product=None):
    """
    This is $F_{sw,all}^{\top,0}$ as defined in A.7.1 of ECHAM6_userguide.pdf

    :param infile: An ECHAM-HAM echam filepath
    :param product:
    :return:
    """
    forcing_file = get_stream_file(infile, 'forcing')
    clean_all_sky_toa_sw_forcing = cis.read_data(infile, 'srad0', product) - \
                                   cis.read_data(forcing_file, 'FSW_TOTAL_TOP', product)
    return clean_all_sky_toa_sw_forcing


def get_clean_all_sky_surface_net_sw(infile, product=None):
    """
    This is $F_{sw,all}^{\bot,0}$ as defined in A.7.1 of ECHAM6_userguide.pdf

    :param infile: An ECHAM-HAM echam filepath
    :param product:
    :return:
    """
    forcing_file = get_stream_file(infile, 'forcing')
    clean_all_sky_toa_sw_forcing = cis.read_data(infile, 'srads', product) - \
                                   cis.read_data(forcing_file, 'FSW_TOTAL_SUR', product)
    return clean_all_sky_toa_sw_forcing


def get_clean_clear_sky_toa_net_sw(infile, product=None):
    """
    This is $F_{sw,clear}^{\top,0}$ as defined in A.7.1 of ECHAM6_userguide.pdf

    :param infile: An ECHAM-HAM echam filepath
    :param product:
    :return:
    """
    forcing_file = get_stream_file(infile, 'forcing')
    clean_clear_sky_toa_sw_forcing = cis.read_data(infile, 'sraf0', product) - \
                                     cis.read_data(forcing_file, 'FSW_CLEAR_TOP', product)
    return clean_clear_sky_toa_sw_forcing


def get_clean_clear_sky_surface_net_sw(infile, product=None):
    """
    This is $F_{sw,clear}^{\bot,0}$ as defined in A.7.1 of ECHAM6_userguide.pdf

    :param infile: An ECHAM-HAM echam filepath
    :param product:
    :return:
    """
    forcing_file = get_stream_file(infile, 'forcing')
    clean_clear_sky_toa_sw_forcing = cis.read_data(infile, 'srafs', product) - \
                                     cis.read_data(forcing_file, 'FSW_CLEAR_SUR', product)
    return clean_clear_sky_toa_sw_forcing


def calc_clean_sky_toa_sw_cloud_forcing(infile, product=None):
    cloud_forcing = get_clean_all_sky_toa_net_sw(infile, product) - \
                    get_clean_clear_sky_toa_net_sw(infile, product)
    return cloud_forcing

