"""
A module for writing Iris cubes to CMOR / AeroCom compatible files

Adapted from script by Zak Kipling
(c) Duncan watson-parris 2017
"""
import cis
import iris
import numpy as np
from utils import get_time_delta


def get_time_freq(time_coord, pdrmip_format):
    """
    Take a time coordinate and return a string describing it's frequency

    :param iris.Coords.Coord time_coord:
    :param bool pdrmip_format: Is the string format required a PDRMIP one?
    :return:
    """
    time_delta = get_time_delta(time_coord)
    # The above call will either return a length one timestep - or raise a ValueError (which we allow to bubble up)
    if time_delta.microseconds:
        raise ValueError("Non-integer number of seconds (+%d us) between timesteps" % time_delta.microseconds)
    elif time_delta.seconds % 3600:
        raise ValueError("Non-integer number of hours (+%d s) between timesteps" % time_delta.seconds)
    elif time_delta.days == 30:
        time_freq = "Amon" if pdrmip_format else "monthly"
    elif time_delta.seconds == 3600:
        time_freq = "hourly"
    elif time_delta.seconds:
        time_freq = "%dhourly" % (time_delta.seconds / 3600)
    elif time_delta.days == 1:
        time_freq = 'day' if pdrmip_format else "daily"
    elif time_delta.days > 1:
        time_freq = "%ddaily" % time_delta.days
    else:
        raise ValueError("Repeating time steps")

    return time_freq


def select_vars(cmor_vars, var_names):
    # This slightly backward method ensures that an error is thrown if the var isn't in the list of cmor_vars
    cmor_var_names = [v.cmor_var_name for v in cmor_vars]
    return [cmor_vars[cmor_var_names.index(v)] for v in var_names]


def get_daily_cubes(c):
    import iris.coord_categorisation
    iris.coord_categorisation.add_day_of_month(c, 'time')
    iris.coord_categorisation.add_month_number(c, 'time')
    
    # Create the aggregation group-by instance.
    groupby = iris.analysis._Groupby([c.coord('day_of_month'), c.coord('month_number')])
    dimension_to_groupby = c.coord_dims(c.coord('day_of_month'))[0]
    cube_slice = [slice(None, None)] * len(c.shape)

    for groupby_slice in groupby.group():
    # for day in c.slices_over('day_of_month'):
        # Put the groupby slice in the right place
        cube_slice[dimension_to_groupby] = groupby_slice
        day = c[tuple(cube_slice)]
        yield day


def get_monthly_cubes(c):
    import iris.coord_categorisation
    iris.coord_categorisation.add_year(c, 'time')
    iris.coord_categorisation.add_month_number(c, 'time')

    # Create the aggregation group-by instance.
    groupby = iris.analysis._Groupby([c.coord('month_number'), c.coord('year')])
    dimension_to_groupby = c.coord_dims(c.coord('month_number'))[0]
    cube_slice = [slice(None, None)] * len(c.shape)

    for groupby_slice in groupby.group():
    # for day in c.slices_over('day_of_month'):
        # Put the groupby slice in the right place
        cube_slice[dimension_to_groupby] = groupby_slice
        month = c[tuple(cube_slice)]
        yield month

def output_monthly_cubes(cube, filename_template='out_{}'):
    for c in get_monthly_cubes(cube):
        # Just take the first value since they're all the same day...
        date = c.coord('time').units.num2date(c.coord('time').points[0])
        date_fmt = date.strftime('%Y%m')
        out = filename_template.format(date_fmt)
        print("Saving to {}...".format(out))
        iris.save(c, out)


class cmor_var:
    def __init__(self, cmor_var_name, load, stream='', long_name=None, standard_name=None,
                 units=None, vertical_coord_type=None, scaling=1.0, comment=None, product=None):
        """
        A class for representing a single cmor valid variable. The attributes (standard_name, long_name etc are only
        needed if they aren't present in the model output (or are different, e.g. in the case of units).

        :param str cmor_var_name: The cmor-approved variable name
        :param callable or str load: Either CIS style variable name(s) or a callable function returning a single cube
        :param str stream: The filestream in which to find the above variable(s)

        The following are optional attributes
        :param str long_name: Variable long-name
        :param str standard_name: CF compliant standard name
        :param cf_units.Unit units: The required output units (which will be converted to)
        :param str vertical_coord_type: The vertical coordinate type (Surface, Column, etc). An educated guess will be
        made if not present
        :param float scaling: A constant scaling to apply
        :param str comment: Any comment attributes to add
        :param str product: CIS plugin name for reading the variable
        """
        self.product = product
        self.comment = comment
        self.scaling = scaling
        self.vertical_coord_type = vertical_coord_type
        self.units = units
        self.standard_name = standard_name
        self.long_name = long_name
        self.stream = stream
        self.load = load
        self.cmor_var_name = cmor_var_name

    def load_var(self, infile, product=None):
        from utils import filename_suffix
        # Add the stream name
        stream_file = filename_suffix(infile, self.stream)
        # Take the specific variable product over the general one if there is one
        product = self.product or product
        if callable(self.load):
            cube = self.load(stream_file, product=product)
        else:
            cube = cis.read_data(stream_file, self.load, product)
        return cube

    def write_var(self, cube, time=None, outbase=None, daily=False, monthly=False, experiment_info=None,
                  contact_info=None, overwrite=False, pdrmip_format=False, output_monthly=False):
        """
        Write the variable to a single file
        :param iris.cube.Cube cube: The data cube to be output
        :param str time: The time period to use for variables that don't have a unique one
        :param str outbase: The output filename base: aerocom3_<!ModelName>_<!ExperimentName>
        :param bool daily: Assume daily output
        :param bool monthly: Assume monthly output
        :param str experiment_info: e.g. "hindcast experiment (1980-2008); ACCMIP-MACCity emissions; nudged to ERAIA.";
        :param str contact_info: e.g. "Nick Schutgens (schutgens\@physics.ox.ac.uk)";
        :param: bool pdrmip_format: Construct filename in the pdrmip way?
        """
        from iris.std_names import STD_NAMES
        import os

        print("Process %s: '%s' / %s" % (self.cmor_var_name, self.long_name, self.standard_name))

        # Do scaling and conversions
        cube *= self.scaling
        if self.units is not None:
            if cube.units is not None and cube.units != '1':
                cube.convert_units(self.units)
            else:
                print("WARNING: Setting units to '{}'without conversion".format(self.units))
                cube.units = self.units

        # Update attributes
        cube.var_name = self.cmor_var_name
        if self.long_name is not None:
            cube.long_name = self.long_name
        if self.standard_name is not None:
            if self.standard_name not in STD_NAMES:
                STD_NAMES[self.standard_name] = {"canonical_units": cube.units}
            cube.standard_name = self.standard_name
        if self.comment is not None:
            cube.attributes['comment'] = self.comment
        if experiment_info is not None:
            cube.attributes['info_exp'] = experiment_info
        if contact_info is not None:
            cube.attributes['info_contact'] = contact_info

        # Figure out the vertical coordinate type
        if self.vertical_coord_type is not None:
            vert_coord = self.vertical_coord_type
        else:
            vert_coord = cube.coords(axis="Z", dim_coords=True)
            if len(vert_coord) == 0:
                vert_coord = cube.coords(axis="Z")
            print("vert_coord: %s" % (", ".join(c.name() for c in vert_coord)))
            if len(vert_coord) == 0:
                if cube.standard_name and cube.standard_name.startswith("atmosphere_boundary_layer_"):
                    vert_coord = "Surface"
                elif cube.standard_name and cube.standard_name.startswith("atmosphere_"):
                    vert_coord = "Column"
                elif cube.standard_name and cube.standard_name.startswith("surface_"):
                    vert_coord = "Surface"
                else:
                    raise ValueError("Unknown vertical coordinate type for %s" % self.cmor_var_name)
            elif len(vert_coord) == 1:
                vert_coord, = vert_coord
                if vert_coord.name() in ["model_level_number",
                                         "atmosphere_hybrid_height_coordinate",
                                         "atmosphere_hybrid_sigma_pressure_coordinate",
                                         "hybrid level at layer midpoints"]:
                    vert_coord = "ModelLevel"
                elif vert_coord.standard_name and vert_coord.standard_name == 'air_pressure':
                    vert_coord = "PressureLevel"
                else:
                    raise ValueError(
                        "Unknown vertical coordinate type (%s) for %s" % (vert_coord.name(), self.cmor_var_name))
            else:
                raise ValueError("Multiple vertical coordinates (%s) for %s" % (
                ", ".join(c.name() for c in vert_coord), self.cmor_var_name))
            print("vert_coord guessed: %s" % vert_coord)

        # Figure out the time period and frequency
        # TODO This doesn't actually shift the time reference, but Nick's script converts to 1850-01-01,00:00:00
        time_coord = cube.coords('time')
        print("time_coord: %s " % (", ".join(c.name() for c in time_coord)))
        if len(time_coord) == 1:
            time_coord, = time_coord
            time_period = np.unique([("%04d%02d%02d" % (d.year, d.month, d.day) if daily
                                         else "%04d%02d" % (d.year, d.month) if monthly
                                         else "%04d" % d.year)
                                        for d in time_coord.units.num2date(time_coord.points)])
            if len(time_period) == 0:
                time_period = time if time is not None else "9999"
            elif len(time_period) > 1:
                if time is not None:
                    time_period = time
                else:
                    raise ValueError("Multiple time periods (%s-%s) for %s" % (
                    min(time_period), max(time_period), self.cmor_var_name))
            else:
                time_period, = time_period
            time_freq = get_time_freq(time_coord, pdrmip_format)
        elif len(time_coord) == 0:
            time_period = time if time else "9999"
            time_freq = 'fx' if pdrmip_format else 'timeinvariant'
        else:
            raise ValueError("Multiple time coordinates (%s) for %s" % (", ".join(time_coord), self.cmor_var_name))
        print("time_period: %s" % time_period)
        print("time_coord: %s" % time_freq)

        if pdrmip_format:
            output_template = "{var}_{freq}_{model_exp}_{period}.nc"
        elif output_monthly:
            # Create an extra set of curly braces for the specific month. but escape them for the first format call
            output_template = "{model_exp}_{var}_{vert}_{{}}_{freq}.nc"
        else:
            output_template = "{model_exp}_{var}_{vert}_{period}_{freq}.nc"

        output_path, model_exp = os.path.split(outbase)
        outfile = os.path.join(output_path, output_template.format(model_exp=model_exp, var=self.cmor_var_name,
                                                                   vert=vert_coord, period=time_period,
                                                                   freq=time_freq))
        print("output filename: %s" % outfile)

        if os.path.isfile(outfile) and not overwrite:
            print("Skipping output as file already exists")
        else:
            if output_monthly:
                output_monthly_cubes(cube, outfile.replace("{{}}", "{}"))
            else:
                iris.save(cube, outfile)
