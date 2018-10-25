# -*- coding: utf-8 -*-

"""
Functions/classes dedicated to inputs/outputs using |pypago|
"""

from __future__ import print_function
import pickle
import re
import numpy as np
from netCDF4 import Dataset
from netcdftime import utime
import pypago.disp


def read_time(filename, time_varname):

    """
    Reads the time variable of a NetCDF time.
    If the time variable contains a 'units'
    attribute, then a list of
    py:class:`datetime.datetime` is returned.
    Else, a list of floats is returned.

    :param str filenaeme: NetCDF file name
    :param str time_varname: Name of the time variable

    """

    # open the netcdf time and extracts the time
    # netCDF variable
    with Dataset(filename) as fin:

        time = fin.variables[time_varname]

        if hasattr(time, 'units'):
            # if the variable has a units, we check the presence
            # of a calendar attribute. If None, we take gregorian
            # as default
            units = getattr(time, 'units')
            if hasattr(time, 'calendar'):
                calendar = getattr(time, 'calendar')
            else:
                message = 'No "calendar" attribute was found.\n'
                message += 'A Gregorian calendar is assumed.\n'
                print(message)
                calendar = 'gregorian'

            # conversion from numeric to datetime
            cdftime = utime(units, calendar)
            output = cdftime.num2date(time[:])

            # barrier.n: if the calendar is a non real one (360, all_leap, etc),
            # the output of the num2date function is not a datetime.datetime object
            # but a netcdftime._datetime.datetime, which can be converted back into
            # real datetime.
            # cf https://ocefpaf.github.io/python4oceanographers/blog/2015/08/10/cf_units_and_time/
            # if isinstance(output[0], 'netcdftime._datetime.datetime'):
            #       output = np.array([d._to_real_datetime() for d in output])

        else:
            # if no units attribute, we return a
            # numpy array of floats
            message = 'No "units" attribute was found.\n'
            message += 'Time is kept in his numeric form instead of dates.\n'
            output = time[:]

    return output


def count_ndim(filename, varname):

    """
    Counts the number of dimensions
    of a variable

    :param str filename: NetCDF file
    :param str varname: Name of the dimension
    :rtype: int
    :return: Number of dimensions

    """

    with Dataset(filename) as fin:
        output = fin.variables[varname].ndim

    return output


def count_dim(filename, varname):

    """
    Counts the number of elements of a
    netcdf file along a given dimension.

    :param str filename: NetCDF file
    :param str varname: Name of the dimension
    :rtype: int
    :return: Number of dimensions

    """

    with Dataset(filename) as fin:
        output = len(fin.dimensions[varname])

    return output


def correct_file_arch(filename):

    """
    This function corrects the architecture of `.pygo` files
    that contains any pago_obj instances.
    It reads the pickle file as a text file
    and virtually replace the module associated with
    the pago_obj by the pypago.io one. A new file
    is created, which has the same
    absolute path as the input, but with `_new` suffix added

    :param str filename: The name of the `.pygo` file to modify
    """

    # regular expression to match the pago_obj
    regexp = re.compile('^pago_obj$')

    # we open the input file (pickle file) and read the lines
    with open(filename, 'r') as fin:
        all_lines = fin.readlines()

    # index = line index of the "pago_obj" string
    index = 0
    # looping over the lines
    # when the regexp is matched we exit the loop
    for line in all_lines:
        if regexp.match(line):
            break
        index = index + 1

    # index is now the line where the
    # module associated with pago_obj is associated
    index = index - 1

    # initialisation of the output filename
    filenameout = filename.replace('.pygo', '_new.pygo')
    with open(filenameout, 'w') as fout:

        # we loop over all the lines of the INPUT file
        for indline in xrange(0, len(all_lines)):
            if indline == index:
                # if we are at the line of the module declaration
                # we replace by the right module
                fout.write('(cpypago.io\n')
            else:
                # else, we copy the input file
                fout.write(all_lines[indline])


def load(filename):

    """
    Loads a file and reads its content.
    The file must have been saved with the
    :py:func:`pypago.io.save` function.
    It uses the :py:func:`pickle.load` function.

    :param str filename: the name of the file (:file:`.pygo`)

    :return: The content of the file

    :rtype: :py:class:`list`, :py:class:`dict` or |pagoobj|

    """

    with open(filename, 'r') as pfile:
        data = pickle.load(pfile)

    return data


def modify_meshfile(filename):

    """
    Masks the vertical scale factors e3t, e3v and e3u on solid points.
    This is useful when using partial steps, since it allows
    to properly mask the solid points along the sections staircases,
    instead of masking where the current velocity is zero.

    :param str filename: Name of the mesh file which contains
       the scale factors to modify

    """

    with Dataset(filename, 'r+') as fin:

        e3t = fin.variables['e3t'][:]
        e3u = fin.variables['e3u'][:]
        e3v = fin.variables['e3v'][:]

        e3t = np.ma.masked_where(e3t != e3t, e3t)
        e3u = np.ma.masked_where(e3u != e3u, e3u)
        e3v = np.ma.masked_where(e3v != e3v, e3v)

        if e3t.mask.ndim == 0:

            if ('umask' in fin.variables.keys()) & \
                   ('vmask' in fin.variables.keys()) & \
                   ('tmask' in fin.variables.keys()):

                umask = fin.variables['umask'][:]
                tmask = fin.variables['tmask'][:]
                vmask = fin.variables['vmask'][:]

                e3t[tmask == 0] = np.nan
                e3u[umask == 0] = np.nan
                e3v[vmask == 0] = np.nan

            elif 'mbathy' in fin.variables.keys():

                mbathy = np.squeeze(fin.variables['mbathy'][:].astype(np.int))

                for i in xrange(0, e3t.shape[2]):
                    for j in xrange(0, e3t.shape[3]):
                        temp = mbathy[i, j]
                        if temp == 0:
                            e3u[:, :, i, j] = np.nan
                            e3v[:, :, i, j] = np.nan
                            e3t[:, :, i, j] = np.nan
                        else:
                            e3t[:, temp:, i, j] = np.nan
                            e3u[:, temp:, i, j] = np.nan
                            e3v[:, temp:, i, j] = np.nan

            else:

                print('No umask, tmask, vmask or mbathy variable ' +
                      'detected in the mesh_file -> Unchanged mesh file')

        else:

            print('The e3u, e3v and e3t arrays seem already masked -> ' +
                  'Unchanged mesh file')

        e3u = np.ma.masked_where(e3u != e3u, e3u)
        e3v = np.ma.masked_where(e3v != e3v, e3v)
        e3t = np.ma.masked_where(e3t != e3t, e3t)

        fin.variables['e3u'][:] = e3u
        fin.variables['e3v'][:] = e3v
        fin.variables['e3t'][:] = e3t


def read_bg_file(finname):

    """
    Reads a background netcdf file (bathy for instance).

    It takes as an input the name of a |netcdf| file, which must contain
    three variables:

    - A variable that has the same name as the first dimension (lat for
      instance, 1D or 2D)
    - A variable that has the same name as the second dimension (lon for
      instance, 1D or 2D)
    - And a third variable, whose name differs from the dimension names
      (bathy for instance, must be 2D)

    It is used in the
    :py:class:`pypago.guis.gui_sections_edition` module.

    .. versionchanged:: 20150813

    :param str finname: name of the background file

    """

    # barrier.n --- modified 2015-08-13
    # no more access to forbidden variable _name
    with Dataset(finname, 'r') as fin:
        dimnames = fin.dimensions.keys()

        for varname in fin.variables.keys():
            if varname not in dimnames:
                toplot = np.squeeze(fin.variables[varname][:])

        ycoord = fin.variables[dimnames[0]][:]
        xcoord = fin.variables[dimnames[1]][:]

        if xcoord.ndim == 1:
            xcoord, ycoord = np.meshgrid(xcoord, ycoord)

    return xcoord, ycoord, toplot


def readnc(filename, varname, start=None, end=None, stride=None):

    """
    Reads a variable from a |netcdf| file.
    Additional arguments allow to read only a part of the variable.
    If `start` and `end` are `None`, all the file is read.
    If they are defined, the `stride` variable is also read.
    If it is `None`, it is set to 1 along each dimension.

    :param str filename: name of the |netcdf| file
    :param str var: name of the variable to read
    :param list start: list that contains the first index to read for each
       dimension.
    :param list end: list that contains the last index to read for each
       dimension.
    :param list stride: list that contains the stride for each dimension
    :return: the array that contains the variable
    :rtype: numpy.array

    """

    with Dataset(str(filename)) as ncfile:

        if (start is not None) & (end is not None):

            start = np.array(start)
            end = np.array(end)

            # check if the end and start arrays have the same dimensions
            if np.any(end.shape != start.shape):
                message = "Start and end must have the same shape!"
                error = pypago.disp.PypagoErrors(message)
                raise error

            # if stride is deended, we check that it has the same size as the
            # start array
            if stride is not None:
                stride = np.array(stride)
                if np.any(start.shape != stride.shape):
                    message = "Start and stride must have the same shape!"
                    error = pypago.disp.PypagoErrors(message)
                    raise error

            # if stride is None, we set it to 1 for each dimensions
            else:
                stride = np.ones(start.shape, dtype=np.int)

            # we call the _readnc_span function
            output = _readnc_span(ncfile, varname, start, end, stride)

        else:
            output = ncfile.variables[varname][:]

    return output


def _readnc_span(ncfile, varname, start, end, stride):

    """ Function which is called if we are to read the variable `var` from \
    the `f` file, with `start`, `end` and `stride` which are not None

    :param netCDF4.Dataset ncfile: netcdf file
    :param str varname: name of the variable to read
    :param numpy.array start: array that contains the start indexes
    :param numpy.array end: array that contains the end indexes
    :param numpy.array stride: array that contains the stride indexes
    """

    vari = ncfile.variables[varname]
    varshape = vari.shape

    # checking that the sizes of the arguments are consistent with the var dimension
    if vari.ndim != len(start):
        message = "The size of the start argument is inconsistent with the variable dimension"
        error = pypago.disp.PypagoErrors(message)
        raise error

    # loop over the end array: where -1,
    # we read to the end of the file
    for idim in xrange(0, len(end)):
        if end[idim] == -1:
            end[idim] = varshape[idim]

    if len(start) == 1:
        output = ncfile.variables[varname][start[0]:end[0]:stride[0]]
    elif len(start) == 2:
        output = ncfile.variables[varname][start[0]:end[0]:stride[0],
                                           start[1]:end[1]:stride[1]]
    elif len(start) == 3:
        output = ncfile.variables[varname][start[0]:end[0]:stride[0],
                                           start[1]:end[1]:stride[1],
                                           start[2]:end[2]:stride[2]]
    elif len(start) == 4:
        output = ncfile.variables[varname][start[0]:end[0]:stride[0],
                                           start[1]:end[1]:stride[1],
                                           start[2]:end[2]:stride[2],
                                           start[3]:end[3]:stride[3]]

    else:
        raise IOError('Number of dimensions %d not implemented' % (len(start)))

    return output


def readnc_dtype(filename, variable):

    """
    Reads the precision of a variable within a file.
    Allows to determine whether the file is in 16 bits
    (hence with possible issues in the calculations).

    :param str filename: the filename which contains the variable
    :param str variable: the name of the variable

    """

    with Dataset(filename) as fin:
        fpres = fin.variables[variable].dtype

    return fpres


def readnc_wrap(filename, var, start, end, nlon):

    """
    Reads a variable from a |netcdf| file, with the possibility to
    do some kind of shiftgrid.
    It uses the :py:func:`pypago.io.readnc` function.

    :param str filename: name of the |netcdf| file
    :param str var: name of the variable to read
    :param list start: list that contains the first index to read for each
       dimension.
    :param list end: list that contains the last index to read for each
       dimension.
    :param list nlon: number of longitudes contained in the file (last
       dimension in the :file:`.nc` file)
    :return: the array that contains the variable
    :rtype: numpy.array

    """

    start = np.array(start).astype(np.int)
    end = np.array(end).astype(np.int)

    if start[-1] < 0:
        startint = np.zeros(start.shape, dtype=np.int) + start
        startint[-1] = nlon - np.abs(start[-1])
        endint = np.zeros(end.shape, dtype=np.int) + end
        endint[-1] = nlon
        res1 = readnc(filename, var, startint, endint)
        startint = np.zeros(start.shape, dtype=np.int) + start
        startint[-1] = 0
        res2 = readnc(filename, var, startint, end)
        res = np.concatenate((res1, res2), axis=-1)

    else:
        if end[-1] > nlon:
            endint = np.zeros(end.shape, dtype=np.int) + end
            endint[-1] = nlon
            res1 = readnc(filename, var, start, endint)
            startint = np.zeros(start.shape, dtype=np.int) + start
            startint[-1] = 0
            endint = np.zeros(end.shape, dtype=np.int) + end
            endint[-1] = end[-1] - nlon
            res2 = readnc(filename, var, startint, endint)
            res = np.concatenate((res1, res2), axis=-1)

        else:
            if end[-1] < start[-1]:
                endint = np.zeros(end.shape, dtype=np.int) + end
                endint[-1] = nlon
                res1 = readnc(filename, var, start, endint)
                startint = np.zeros(start.shape, dtype=np.int) + start
                startint[-1] = 0
                res2 = readnc(filename, var, startint, end)
                res = np.concatenate((res1, res2), axis=-1)

            else:
                res = readnc(filename, var, start, end)

    return res


def save(dictio, filename):

    """
    Function that saves a variable object into a file. It uses
    the :py:func:`pickle.dump` function

    :param dictio: the variable to save (:py:class:`dict`,
        :py:class:`list` or |pagoobj|)
    :param str filename: the name of the output file
        (use :file:`.pygo` for instance)
    :return: none

    """

    with open(filename, "w") as fileout:
        pickle.dump(dictio, fileout)


def check_ncvar_exist(filename, varname):

    """
    Determines whether a variable exists in a NetCDF file

    :param str filename: NetCDF filename
    :param str varname: NetCDF variable name
    :return: True if the variable is in the file, else False
    :rtype: bool

    """

    with Dataset(filename) as fin:
        output = varname in fin.variables

    return output


if __name__ == '__main__':

    FILENAME_TEST = "examples/data/Omon_CNRM-CM5_piControl_gridT.nc"
    VARNAME_TEST = 'time'
    DATES = read_time(FILENAME_TEST, VARNAME_TEST)
    print(DATES)
