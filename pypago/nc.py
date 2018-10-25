# -*- coding: utf-8 -*-

"""
Module that contains various functions dedicated to the \
conversion of |pypago| outputs into NetCDF files` \
to :file:`.nc` files
"""

from __future__ import print_function
import numpy as np
from netcdftime import utime
from netCDF4 import Dataset
from pypago.misc import PypagoErrors
import pypago.pyio
import pypago.disp
from pypago.sample_param import dictvname
try:
    from param import dictvname2
    dictvname.update(dictvname2)
except ImportError:
    pass


def gridsec_tonc(finname, varname, section_names=None, units='days since 1900-01-01 00:00:00', calendar='gregorian'):

    """
    Conversion of tracer and velocities along the
    sections into |netcdf|. There will
    be one file per section, in which the
    variables will be saved. The absolute path of
    the output files will be the same
    as for the input file name, except that
    the :file:`.pygo` will be replaced by
    :file:`_sec_{SECTIONNAME}.nc`

    :param str finname: output of the
     :py:func:`pypago.loaddata.loaddata` function (must contain
     a dictionary with the `MODEL_sections` and `MODEL_time` keys)

    :param list section_names: Default is None.
     If set, the list of the sections' names    from which to extract the data

    :param list varname: Default is None. If None, all the variables (`vect`,
     `vecs` and `vecv` are extracted.  If set, the list of the variables
     names to extract

    """

    modelsec = pypago.pyio.load(finname)

    print(r'Conversion of the %s file into netcdf' % finname)
    print(r'==============================================')

    if section_names is None:

        print(r'Processing of ALL the sections')

        for secint in modelsec:
            _write_gridsec_netcdf(finname, secint, varname, units, calendar)
            print('section %s: done' % secint.name)

    else:

        # Conversion of section_names into list
        if not isinstance(section_names, 'list'):
            section_names = [section_names]

        print(r'Processing of the section list %s'
              % section_names)

        for secnames in section_names:
            indice_sec = pypago.misc.findsecnum(modelsec, secnames)
            secint = modelsec[indice_sec]
            _write_gridsec_netcdf(finname, secint, varname, units, calendar)
            print('section %s: done' % secint.name)

    print('\n')


def gridvol_tonc(finname, varname, domain_names=None, units='days since 1900-01-01 00:00:00', calendar='gregorian'):

    """
    Conversion of tracer fields within a domain into |netcdf|.
    There will be one file per domain, in which the variables will be saved.
    The absolute path of the output files will be the same
    as for the input file name, except that the :file:`.pygo`
    will be replaced by
    :file:`_dom_{DOMAINNAME}.nc`

    :param str finname: output of the
        :py:func:`pypago.loaddata.loaddata` function (must contain
        a dictionary with the `MODEL_areas` and `MODEL_time` keys)

    :param list domain_names: Default is None. If set, the list of the domains'
        names from which to extract the data

    :param list varname: Default is None. If None, all the variables
        (`temperature`, `salinity`) are extracted.
        If set, the list of the variables names to extract

    """

    modeldom = pypago.pyio.load(finname)
    # modeldom = data['MODEL_areas']
    # modeltime = data['MODEL_time']

    print(r'Conversion of the %s file into netcdf' % finname)
    print(r'==============================================')

    if domain_names is None:

        print(r'Processing of ALL the domains')

        for domint in modeldom:
            _write_gridvol_netcdf(finname, domint, varname, units, calendar)
            print('domain %s: done' % domint.name)

    else:

        # Conversion of domain_names into list
        if not isinstance(domain_names, 'list'):
            domain_names = [domain_names]

        print(r'Processing of the domain list %s'
              % domain_names)

        for domnames in domain_names:

            indice_dom = pypago.misc.finddomnum(modeldom, domnames)
            domint = modeldom[indice_dom]
            _write_gridvol_netcdf(finname, domint, varname, units, calendar)
            print('domain %s: done' % domint.name)

    print('\n')


def secind_tonc(finname, section_names=None):

    """
    Conversion of section indices (issued from
    :py:func:`pypago.pypago.indices_MODEL` function) in |netcdf|.

    :param str finname: output of the
       :py:func:`pypago.pypago.indices_MODEL` function

    :param str section: Default is None. If set, the section to extract
    """

    data = pypago.pyio.load(finname)
    sectionind = data['MODEL_indices']
    modeltime = data['MODEL_time']

    print(r'Conversion of the %s file into netcdf' % finname)
    print(r'==========================================')

    if section_names is None:

        print(r'Processing of ALL the sections')

        for secint in sectionind:
            _write_secind_netcdf(finname, modeltime, secint)
            print('section %s: done' % secint.name)

    else:
        print(r'Processing of the section list %s' % section_names)

        for secnames in section_names:

            indice_section = pypago.misc.findsecnum(sectionind, secnames)
            secint = sectionind[indice_section]
            _write_secind_netcdf(finname, modeltime, secint)
            print('section %s: done' % secint.name)

    print('\n')


def sections_tonc(finname, section_names=None):

    """
    Extraction of sections' endpoints into |netcdf| files

    :param str finname: filename containing the list of sections' endpoints
       (output of the :py:mod:`pypago.guis.gui_sections_edition`)

    :param list section_names: Default is None. If set, the list of sections'
       names from which to extract the endpoints

    Convert all sections of :file:`sections_NA_Nico.pygo` to individual
    |netcdf| files:
    >>> _sections_tonc('sections_NA_Nico.pygo')

    Convert only ``ar7`` section
    >>> _sections_tonc('sections_NA_Nico.pygo', section_names=['ar7'])

    What append if missing file ?
    >>> _sections_tonc('badfile')

    What append if missing section name in the file ?
    >>> _sections_tonc('sections_NA_Nico.pygo', section_names=['badsection'])

    """

    try:
        pagosections = pypago.pyio.load(finname)
    except IOError:
        message = r'filename {0} not found'.format(finname)
        error = PypagoErrors(message)
        raise error

    print(r'Conversion of the %s file into netcdf' % finname)
    print(r'==============================================')

    if section_names is None:

        print(r'Processing of ALL the sections')

        for secint in pagosections:
            _write_sections_netcdf(finname, secint)
            print('section %s: done' % secint.name)

    else:

        print(r'Processing of the section list %s' % section_names)

        for secnames in section_names:

            try:
                indice_sec = pypago.misc.findsecnum(pagosections, secnames)
            except ValueError:
                message = r'Can not produce NetCDF file'
                error = PypagoErrors(message)
                raise error

            secint = pagosections[indice_sec]
            _write_sections_netcdf(finname, secint)
            print('section %s: done' % secint.name)

    print('\n')


def _volind_tonc(finname, domain_names=None):

    """
    Conversion of volume indices (issued from
    :py:func:`pypago.pypago.volumes_MODEL` function) in |netcdf|.

    :param str finname: output of
       the :py:func:`pypago.pypago.volumes_MODEL` function

    :param str domain: Default is None. If set, the domain to extract
    """

    data = pypago.pyio.load(finname)
    volumeind = data['MODEL_volumes']
    modeltime = data['MODEL_time']

    print(r'Conversion of the %s file into netcdf' % finname)
    print(r'==========================================')

    if domain_names is None:

        print(r'Processing of ALL the domains')

        for domint in volumeind:
            _write_volind_netcdf(finname, modeltime, domint)
            print('domain %s: done' % domint.name)

    else:

        print(r'Processing of the domain list %s' % domain_names)

        for domnames in domain_names:
            indice_dom = pypago.misc.finddomnum(volumeind, domnames)
            domint = volumeind[indice_dom]
            _write_volind_netcdf(finname, modeltime, domint)

    print('\n')


def _write_gridsec_netcdf(finname, secint, varname, units, calendar):

    """
    Function that handles the writing of individual |netcdf| section files
    (called by :py:func:`pypago.tonc.gridsec_tonc` function)

    :param str finname: output of the
       :py:func:`pypago.pypago.sections_MODEL` function
    :param float modeltime: time vector contained on the file
    :param pypago.sections.GridSection secint: |pypago| object that contains the
       index elements
    :param str varname: variable name to save in the file ('vecv', 'vect' or
       'vecv'. If None, all are saved)

    """

    timename = dictvname['time_varname']

    name = secint.name
    nz, nl = secint.areavect.shape
    nlvect = len(secint.lvect)

    foutname = finname.replace('.pygo', '_sec_'+name+'.nc')
    fout = Dataset(foutname, 'w')

    fout.createDimension('z', nz)
    fout.createDimension('l', nl)
    fout.createDimension('lvect', nlvect)
    fout.createDimension(timename, None)

    # If time exists in the section attribute:
    if hasattr(secint, timename):
        # create time variable
        fout.createVariable(timename, 'f', (timename, ))

        # recover time variable
        modeltime = secint.time

        # if the variable is not made of numbers: conversion into
        # numerical time
        if not isinstance(modeltime[0], (int, long, float)):
            cdftime = utime(units, calendar)
            modeltime = cdftime.date2num(modeltime)
            fout.variables[timename].calendar = calendar
            fout.variables[timename].units = units
        fout.variables[timename][:] = modeltime

    fout.createVariable('areavect', 'f', ('z', 'l'))
    fout.variables['areavect'][:] = secint.areavect

    fout.createVariable('depthvect', 'f', ('z', 'l'))
    fout.variables['depthvect'][:] = secint.depthvect

    fout.createVariable('veci', 'i', ('l', ))
    fout.variables['veci'][:] = secint.veci

    fout.createVariable('vecj', 'i', ('l', ))
    fout.variables['vecj'][:] = secint.vecj

    fout.createVariable('orient', 'i', ('l', ))
    fout.variables['orient'][:] = secint.orient

    facesout = np.empty(secint.faces.shape)
    facesout[secint.faces == 'N'] = 1
    facesout[secint.faces == 'W'] = 2

    fout.createVariable('faces', 'i', ('l', ))
    fout.variables['faces'][:] = facesout
    fout.variables['faces'].description = '1 for N, 2 for W'

    if not isinstance(varname, list):
        varname = [varname]

    for var in varname:
        fout.createVariable(var, 'f', (timename, 'z', 'l',))
        fout.variables[var][:] = secint.__dict__[var]

    fout.close()


def _write_gridvol_netcdf(finname, domint, varname, units, calendar):

    """
    Function that handles the writing of individual |netcdf| volume files
    (called by :py:func:`pypago.tonc.gridvol_tonc` function)

    :param str finname: output of the
       :py:func:`pypago.pypago.sections_MODEL` function
    :param float modeltime: time vector contained on the file
    :param pypago.areas.Areas domint: |pypago| object that contains the
       index elements
    :param str varname: variable to save in the file ('salinity' or
       'temperature'). If None, both are saved
    """

    timename = dictvname['time_varname']

    nz, npoints = domint.volume.shape

    foutname = finname.replace('.pygo', '_dom_' + domint.name + '.nc')
    fout = Dataset(foutname, 'w')

    fout.createDimension(timename, None)
    fout.createDimension('z', nz)
    fout.createDimension('points', npoints)

    # If time exists in the section attribute:
    if hasattr(domint, timename):
        # create time variable
        fout.createVariable(timename, 'f', (timename, ))

        # recover time variable
        modeltime = domint.time

        # if the variable is not made of numbers: conversion into
        # numerical time
        if not isinstance(modeltime[0], (int, long, float)):
            cdftime = utime(units, calendar)
            modeltime = cdftime.date2num(modeltime)
            fout.variables[timename].calendar = calendar
            fout.variables[timename].units = units
        fout.variables[timename][:] = modeltime

    fout.createVariable('indi', 'i', ('points',))
    fout.createVariable('indj', 'i', ('points', ))
    fout.createVariable('volume', 'f', ('z', 'points', ))
    fout.createVariable('surface', 'f', ('points', ))

    fout.variables['indi'][:] = domint.i
    fout.variables['indj'][:] = domint.j
    fout.variables['volume'][:] = domint.volume
    fout.variables['surface'][:] = domint.surface

    # If the domain has been defined from sections,
    # creation of section dimension and variables
    if domint.secnames is not None:
        nsec = len(domint.secnames)
        fout.createDimension('section', nsec)
        fout.createVariable('sectionsign', 'i', ('section', ))
        fout.variables['sectionsign'][:] = domint.signs
        fout.variables['sectionsign'].sections = ','.join(domint.secnames)

    # barrier.n --- modified 2015-08-13
    # this variable is not written in the file, since its
    # dimensions are (lat, lon). To keep it, we must add
    # these two dimensions in the file
    # fout.variables['mask'][:] = domint.mask
    # fout.createVariable('mask', 'i', ('points', ))

    if not isinstance(varname, list):
        varname = [varname]

    for var in varname:
        if domint.__dict__[var].ndim == 3:
            fout.createVariable(var, 'f', (timename, 'z', 'points'))
            fout.variables[var][:] = domint.__dict__[var]
        elif domint.__dict__[var].ndim == 2:
            fout.createVariable(var, 'f', (timename, 'points'))
            fout.variables[var][:] = domint.__dict__[var]
        else:
            print('Not implemented !')

    fout.close()


def _write_secind_netcdf(finname, modeltime, secint):

    """
    Function that handles the writing of individual |netcdf| index files
    (called by :py:func:`pypago.tonc.secind_tonc` function)

    :param str finname: output of the
       :py:func:`pypago.pypago.indices_MODEL` function
    :param float modeltime: time vector contained on the file
    :param pypago.sections.GridSection secint: |pypago| object that contains the
       volume elements
    """

    name = secint.name
    nl = secint.bpmt_vec.shape[1]

    foutname = finname.replace('.pygo', '_sec_' + name + '.nc')
    fout = Dataset(foutname, 'w')

    fout.createDimension('time', None)
    fout.createDimension('l', nl)

    fout.createVariable('time', 'f', ('time', ))
    fout.variables['time'] = modeltime

    fout.createVariable('lindex', 'i', ('l',))
    fout.variables['lindex'] = np.arange(0, nl)

    # barrier.n -- modified 2015-08-07
    # instead of using exceptions, we check for
    # the variable type and the variable dimension
    for key, val in zip(secint.__dict__.keys(),
                        secint.__dict__.values()):

        if isinstance(val, np.ndarray):
            if val.ndim == 2:  # for the _vec variables
                fout.createVariable(key, 'f', ('time', 'l'))
                fout.variables[key][:] = val
            else:
                fout.createVariable(key, 'f', ('time', ))
                fout.variables[key][:] = val


def _write_sections_netcdf(finname, secint):

    """
    Function that handles the writing of individual section endpoints files
    (called by :py:func:`pypago.tonc.sections_tonc` function)

    :param str finname: output of the
       :py:mod:`pypago.guis.ihm_editions_sections` module
    :param pypago.sections.Section secint: |pypago| object that contains the
       section element
    """

    name = secint.name
    lon = secint.lon
    lat = secint.lat
    dire = secint.dir

    dire = np.array(dire)
    dirout = np.empty(dire.shape)
    dirout[dire == 'NE'] = 1
    dirout[dire == 'NW'] = 2
    dirout[dire == 'SE'] = 3
    dirout[dire == 'SW'] = 4

    npoint = len(lon)
    nseg = len(dirout)

    foutname = finname.replace('.pygo', '_'+name+'.nc')
    fout = Dataset(foutname, 'w')
    fout.createDimension('npoint', npoint)
    fout.createDimension('nseg', nseg)
    fout.createVariable('lon', 'f', ('npoint', ))
    fout.createVariable('lat', 'f', ('npoint', ))
    fout.createVariable('dir', 'i', ('nseg', ))
    fout.variables['lon'][:] = lon
    fout.variables['lat'][:] = lat
    fout.variables['dir'][:] = dirout
    fout.variables['dir'].description = '1 for NE, 2 for NW, 3 \
    for SE and 4 for SW'

    fout.close()


def _write_volind_netcdf(finname, modeltime, domint):

    """
    Function that handles the writing of individual |netcdf| volume files
    (called by :py:func:`pypago.tonc.volind_tonc` function)

    :param str finname: output of the
       :py:func:`pypago.pypago.volumes_MODEL` function
    :param float modeltime: time vector contained on the file
    :param pypago.areas.Areas domint: |pypago| object that contains
        the volume index elements
    """

    name = domint.name

    foutname = finname.replace('.pygo', '_dom_'+name+'.nc')
    fout = Dataset(foutname, 'w')

    fout.createDimension('time', None)
    fout.createVariable('time', 'f', ('time', ))
    fout.variables['time'] = modeltime

    # barrier.n -- modified 2015-08-07
    # instead of using exceptions, we check for
    # the variable type and the variable dimension
    for key, val in zip(domint.__dict__.keys(),
                        domint.__dict__.values()):

        if isinstance(val, np.ndarray):
            fout.createVariable(key, 'f', ('time',))
            fout.variables[key][:] = val

    fout.close()
