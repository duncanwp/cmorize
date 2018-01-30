"""
(c) Duncan watson-parris 2017
"""


def get_stream_file(infile, stream):
    """
    Given an ECHAM filepath get a different stream
    """
    # This is a bit of a hack really, but should work ~95% of the time
    return "_".join(infile.split('_')[:-1] + [stream+'.nc'])


def get_time_delta(time_coord):
    """
    Return the unique timestep from a time coordinate, or the non-unique timesteps for monthly data

    :param cis.data_id.coords.Coord time_coord:
    :return ndarry: Array of timesteps
    :raises ValueError when the step can't be determined or a non-regular period is found
    """
    import datetime
    import numpy as np
    time_delta = np.unique(np.diff(time_coord.units.num2date(time_coord.points)))
    if len(time_delta) == 0:
        if time_coord.has_bounds():
            time_delta, = np.diff(time_coord.units.num2date(time_coord.bounds))[0]
        else:
            raise ValueError("Cannot guess time step from a single one without bounds")
    elif len(time_delta) == 1:
        time_delta = time_delta[0]
    elif (np.amin(time_delta) >= datetime.timedelta(days=28) and
                np.amax(time_delta) <= datetime.timedelta(days=31)):
        # Monthly timedelta
        time_delta = datetime.timedelta(days=30)
    else:
        raise ValueError("Non-uniform period (%g to %g) between timesteps" % (
            np.amin(time_delta), np.amax(time_delta)))
    return time_delta


def filename_suffix(f, suffix):
    """
    Append a suffix to a filename, before the extension
    :param str f: Filename (and optionally path)
    :param str suffix: The suffix
    :return str: The full filename with new suffix
    """
    from os.path import splitext
    f, ext = splitext(f)
    return f + suffix + ext
