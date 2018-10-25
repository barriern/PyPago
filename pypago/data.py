
""" Module for the reading of data """

from __future__ import print_function
import numpy as np
import pypago.pyio
import pypago.misc
from pypago.disp import PypagoErrors
from pypago.sample_param import dictvname
try:
    from param import dictvname2
    dictvname.update(dictvname2)
except ImportError:
    pass

def loadtime(model_sections, filename):

    """
    Load the time variable (using the :py:func:`pypago.pyio.read_time`
    function) and add it to the model sections attributes (time attribute)

    :param list model_sections: List of model sections
    :param str filename: NetCDF filename
    :return: List of model sections with the time attribute added
    :rtype: list

    """

    time_varname = dictvname['time_varname']   # output name of the time variable
    time = pypago.pyio.read_time(filename, time_varname)    # reads the time from the netcdf file

    for secint in model_sections:
        if not hasattr(secint, time_varname):
            setattr(secint, time_varname, time)   # if the 'time' attribute does not exist, it is added to data structure file
        else:
            setattr(secint, time_varname, np.concatenate((getattr(secint, time_varname), time), axis=0))    # else, it is concatenated to the existing one

    return model_sections


def loaddata_sec_t(model_sections, filename, dictvarname):

    """
    Extracts the model output from model T-grid points and interpolates
    them on section west and north faces. The extracted data are added to
    the input file::

        import pypago.data

        filenameT = 'g88_2000_00_gridT.nc'
        dictvarname = {'vect':'votemper, 'vecs':vosaline'}
        pypago.data.loaddata_sec_T('data.pygo', filenameT, dictvarname)

    :param str structfile: The name of the .pygo file
     where the list of :py:class:`pypago.sections.Gridsection` objects
     are stored.

    .. versionchanged:: 20120515

    JD PAGO WHOI
    """

    #   extraction of sub-grid corners and of the total number of longitudinal points
    nlon = model_sections[0].nlon
    jmin = model_sections[0].jmin
    jmax = model_sections[0].jmax
    imin = model_sections[0].imin
    imax = model_sections[0].imax
    nzlev = model_sections[0].areavect.shape[0]

    # count the total number of time steps in the file, over which to loop
    time_varname = dictvname['time_varname']
    ntime = pypago.pyio.count_dim(filename, time_varname)

    for timestep in xrange(0, ntime):

        # loop over the list of variables to extract (on T points)
        for outname, varname in dictvarname.items():

            if pypago.pyio.check_ncvar_exist(filename, varname):

                # extraction of the temperature on T points, but considering one more point left and top
                # to interpolate on W and E faces.
                tempc = pypago.pyio.readnc_wrap(filename, varname,
                                                [timestep, 0, jmin, imin-1],
                                                [timestep+1, nzlev+1, jmax+2, imax+1], nlon)
                tempc = np.ma.masked_where(tempc == 0, tempc)
                tempc = np.ma.masked_where(np.ma.getmaskarray(tempc), tempc)
                # interpolation on N/W faces.
                tempn = np.mean(np.ma.array((tempc[:, :, :-1, 1:], tempc[:, :, 1:, 1:])), axis=0)
                tempw = np.mean(np.ma.array((tempc[:, :, :-1, :-1], tempc[:, :, :-1, 1:])), axis=0)

                # extraction on the real T domain
                tempc = tempc[:, :, :-1, 1:]

                # function that modifies the MODEL_sections list
                _extract_sections(model_sections, tempn, tempw, outname, use_orient=False)

    return model_sections


def loaddata_sec_uv(model_sections, filename_u, filename_v, dictvarname):

    """

    Extracts the model output, interpolates if needed on west and north faces,
    and saves only the data along preselected sections and areas in structfile.

    .. warning::

       FOR NEMO MODEL ONLY

    :param str file_location: is where the model output can be
       found (full path), which name is
       :file:`{file_prefix}{SALT|TEMP|UVEL|VVEL}.{file_suffix}.nc`

    :param bool loadarea: the `loadarea` option determines whether to load
       the temperature and salinity on the volumes computed
       using the :py:func:`pypago.grid.areas_MODEL` function

    .. versionchanged:: 20120515

       JD PAGO WHOI
    """

    # model_sections = pypago.pyio.load(structfile)

    # extraction of sub-domain indexes.
    nlon = model_sections[0].nlon
    jmin = model_sections[0].jmin
    jmax = model_sections[0].jmax
    imin = model_sections[0].imin
    imax = model_sections[0].imax
    nzlev = model_sections[0].areavect.shape[0]

    # time steps in the file over which to loop
    time_varname = dictvname['time_varname']
    ntime = pypago.pyio.count_dim(filename_u, time_varname)

    for timestep in xrange(0, ntime):

        print ('time step = ', timestep)

        # loop over the variables which will be extracted
        for outname, (varname_u, varname_v) in dictvarname.items():

            if pypago.pyio.check_ncvar_exist(filename_u, varname_u) & pypago.pyio.check_ncvar_exist(filename_v, varname_v):

                # originally in centimeters/s >> convert to meters/s
                uvelw = pypago.pyio.readnc_wrap(filename_u, varname_u,
                                                [timestep, 0, jmin, imin-1],
                                                [timestep+1, nzlev+1, jmax+1, imax], nlon)

                # originally in centimeters/s >> convert to meters/s
                vveln = pypago.pyio.readnc_wrap(filename_v, varname_v,
                                                [timestep, 0, jmin, imin],
                                                [timestep+1, nzlev+1, jmax+1, imax+1], nlon)

                # function that modifies the MODEL_sections list
                _extract_sections(model_sections, vveln, uvelw, outname, use_orient=True)

            else:
                print("========================================")
                print("The %s or %s variable does not exist!!!" % (varname_u, varname_v))
                print("========================================")

    return model_sections


def _extract_sections(model_sections, tempn, tempw, outname, use_orient):

    """

    Extracts T, S, U and V along all the sections contained in
    the `model_sections` list, which is then updated.

    """

    for secint in model_sections:

        nzlev = secint.areavect.shape[0]

        veci = secint.veci.astype(np.int)
        vecj = secint.vecj.astype(np.int)
        faces = secint.faces

        vect = np.ma.zeros((1, nzlev, len(veci)))

        nl = len(veci)
        for l in xrange(0, nl):
            if faces[l] == 'N':
                vect[:, :, l] = tempn[:, :, vecj[l], veci[l]]
            elif faces[l] == 'W':
                vect[:, :, l] = tempw[:, :, vecj[l], veci[l]]
            if use_orient:
                vect[:, :, l] *= secint.orient[l]

        # mask data where areavect is NaN
        # first, need to broadcast arrays for them to have the same
        # same size.
        areavect, vect = np.broadcast_arrays(secint.areavect, vect)
        vect[np.isnan(areavect)] = np.nan

        # if len(indnan[0] > 0):
        #     i, j = np.nonzero(np.isnan(areavect))
        #     vect[:, i, j] = np.nan
        # else:
        #     if filepres == np.float32:
        #         #print 'Masking vecv where vecv==0'
        #         #vecv[vecv == 0] = np.nan
        #         #areavect[np.isnan(vecv[0, :, :])] = np.nan
        #         #secint.areavect = areavect
        #     else:
        #         #print r'data precision is unsufficient: vecv will be inappropriately masked.'
        #         #print r'expected errors on indices calculations!'
        #         #vecv[vecv == 0] = np.nan

        # vect[np.ma.getmaskarray(vect)] = np.nan

        if not hasattr(secint, outname):
            setattr(secint, outname, vect)
        else:
            setattr(secint, outname, np.concatenate((getattr(secint, outname), vect), axis=0))


def _trans_to_vel(model_sections):

    """
    Converts the transport arrays (in kg/s) into
    velocity arrays by diving the transport across each cell
    by the surface of the grid cell
    """

    for secint in model_sections:
        secint.vecv[-1, :, :] = secint.vecv[-1, :, :]/secint.areavect*1e-3


def loaddata_area_t(model_areas, filename, dictvarname):

    """
    Extracts the model output from model T-grid points and interpolates
    them on section west and north faces. The extracted data are added to
    the input file::

        import pypago.data

        filenameT = 'g88_2000_00_gridT.nc'
        dictvarname = {'vect':'votemper, 'vecs':vosaline'}
        pypago.data.loaddata_sec_T('data.pygo', filenameT, dictvarname)

    :param str structfile: The name of the .pygo file
     where the list of :py:class:`pypago.sections.Gridsection` objects
     are stored.

    .. versionchanged:: 20120515

       JD PAGO WHOI
    """

    nlon = model_areas[0].nlon
    jmin = model_areas[0].jmin
    jmax = model_areas[0].jmax
    imin = model_areas[0].imin
    imax = model_areas[0].imax
    nzlev = model_areas[0].volume.shape[0]

    time_varname = dictvname['time_varname']
    ntime = pypago.pyio.count_dim(filename, time_varname)

    for timestep in xrange(0, ntime):

        for outname, varname in dictvarname.items():

            if pypago.pyio.check_ncvar_exist(filename, varname):

                vardim = pypago.pyio.count_ndim(filename, varname)

                if vardim == 4:
                    tempc = pypago.pyio.readnc_wrap(filename, varname,
                                                    [timestep, 0, jmin, imin-1],
                                                    [timestep+1, nzlev+1, jmax+2, imax+1], nlon)
                    tempc = tempc[:, :, :-1, 1:]
                elif vardim == 3:
                    tempc = pypago.pyio.readnc_wrap(filename, varname,
                                                    [timestep, jmin, imin-1],
                                                    [timestep+1, jmax+2, imax+1], nlon)
                    tempc = tempc[:, :-1, 1:]
                else:
                    message = 'The number of dimensions within the variable is not '
                    message += 'right. Vardim = %d instead of 3 or 4.' % vardim
                    raise PypagoErrors(message)

                tempc = np.ma.masked_where(tempc == 0, tempc)
                tempc = np.ma.masked_where(np.ma.getmaskarray(tempc), tempc)

                _extract_areas(model_areas, tempc, outname)

                # function that modifies the MODEL_sections list
                # _extract_sections(model_sections, tempn, tempw, outname, use_orient=False)

            else:
                print("=================================")
                print("The %s variable does not exist!!!" % varname)
                print("=================================")

    return model_areas


def _extract_areas(model_areas, tempc, outname):

    """
    Extracts grid T variable on all the areas contained in
    the `MODEL_areas` list, which is then updated.
    """

    # loop over all the areas
    for areaint in model_areas:

        # set to NaN all the missing values of the variable
        mask = pypago.misc.extract_mask(tempc)
        tempc[mask] = np.nan

        # extract the data on the area indexes
        if tempc.ndim == 4:
            tempc = tempc[:, :, areaint.i, areaint.j]
        else:
            tempc = tempc[:, areaint.i, areaint.j]

        if not hasattr(areaint, outname):
            # if outname is not an attribute, we add it
            setattr(areaint, outname, tempc)
        else:
            # if outname is already an attribute, we update it
            setattr(areaint, outname, np.concatenate((getattr(areaint, outname), tempc), axis=0))


def loaddata_sec_uv_roms(model_sections, filename_u, filename_v, dictvarname):

    """

    Extracts the model output, interpolates if needed on west and north faces,
    and saves only the data along preselected sections and areas in structfile.

    .. warning::

       FOR ROMS MODEL ONLY

    :param str file_location: is where the model output can be
       found (full path), which name is
       :file:`{file_prefix}{SALT|TEMP|UVEL|VVEL}.{file_suffix}.nc`

    :param bool loadarea: the `loadarea` option determines whether to load
       the temperature and salinity on the volumes computed
       using the :py:func:`pypago.grid.areas_MODEL` function

    .. versionchanged:: 20120515

       JD PAGO WHOI
    """

    # extraction of sub-domain indexes.
    nlon = model_sections[0].nlon
    jmin = model_sections[0].jmin
    jmax = model_sections[0].jmax
    imin = model_sections[0].imin
    imax = model_sections[0].imax
    nzlev = model_sections[0].areavect.shape[0]

    # time steps in the file over which to loop
    time_varname = dictvname['time_varname']
    ntime = pypago.pyio.count_dim(filename_u, time_varname)

    for timestep in xrange(0, ntime):

        print("Time step %d/%d" % (timestep, ntime))

        # loop over the variables which will be extracted
        for outname, (varname_u, varname_v) in dictvarname.items():

            if pypago.pyio.check_ncvar_exist(filename_u, varname_u) & pypago.pyio.check_ncvar_exist(filename_v, varname_v):

                # originally in centimeters/s >> convert to meters/s
                uvelw = pypago.pyio.readnc_wrap(filename_u, varname_u,
                                                [timestep, 0, jmin, imin],
                                                [timestep+1, nzlev+1, jmax+1, imax+1], nlon)

                # originally in centimeters/s >> convert to meters/s
                # note that V variables are stored on southern faces -> needs to shift one cell up
                vveln = pypago.pyio.readnc_wrap(filename_v, varname_v,
                                                [timestep, 0, jmin+1, imin],
                                                [timestep+1, nzlev+1, jmax+2, imax+1], nlon)

                # function that modifies the MODEL_sections list
                _extract_sections(model_sections, vveln, uvelw, outname, use_orient=True)

            else:
                print("========================================")
                print("The %s or %s variable does not exist!!!" % (varname_u, varname_v))
                print("========================================")

    return model_sections


def loaddata_sec_t_space(model_sections, filename, dictvarname):

    """
    Extracts the model output at the surface from model T-grid points and interpolates
    them on section west and north faces. The extracted data are added to
    the input file::

        import pypago.data

        filenameT = 'g88_2000_00_gridT.nc'
        dictvarname = {'vect':'votemper, 'vecs':vosaline'}
        pypago.data.loaddata_sec_T('data.pygo', filenameT, dictvarname)

    :param str structfile: The name of the .pygo file
     where the list of :py:class:`pypago.sections.Gridsection` objects
     are stored.

    .. versionchanged:: 20120515

    JD PAGO WHOI
    """

    #   extraction of sub-grid corners and of the total number of longitudinal points
    nlon = model_sections[0].nlon
    jmin = model_sections[0].jmin
    jmax = model_sections[0].jmax
    imin = model_sections[0].imin
    imax = model_sections[0].imax

    # time steps in the file over which to loop
    time_varname = dictvname['time_varname']
    ntime = pypago.pyio.count_dim(filename, time_varname)

    for timestep in xrange(0, ntime):

        print("Time step %d/%d" % (timestep, ntime))

        # loop over the list of variables to extract (on T points)
        for outname, varname in dictvarname.items():

            if pypago.pyio.check_ncvar_exist(filename, varname):

                # extraction of the temperature on T points, but considering one more point left and top
                # to interpolate on W and E faces.
                tempc = pypago.pyio.readnc_wrap(filename, varname,
                                                [timestep, jmin, imin-1],
                                                [timestep+1, jmax+2, imax+1], nlon)
                tempc = np.ma.masked_where(tempc == 0, tempc)
                tempc = np.ma.masked_where(np.ma.getmaskarray(tempc), tempc)
                # interpolation on N/W faces.
                tempn = np.mean(np.ma.array((tempc[:, :-1, 1:], tempc[:, 1:, 1:])), axis=0)
                tempw = np.mean(np.ma.array((tempc[:, :-1, :-1], tempc[:, :-1, 1:])), axis=0)

                # extraction on the real T domain
                tempc = tempc[:, :-1, 1:]

                # function that modifies the MODEL_sections list
                _extract_sections_space(model_sections, tempn, tempw, outname, use_orient=False)

    return model_sections


def _extract_sections_space(model_sections, tempn, tempw, outname, use_orient):

    """

    Extracts T variables at the surface along all the sections contained in
    the `model_sections` list, which is then updated.

    """

    for secint in model_sections:

        veci = secint.veci.astype(np.int)
        vecj = secint.vecj.astype(np.int)
        faces = secint.faces

        vect = np.ma.zeros((1, len(veci)))

        nl = len(veci)
        for l in xrange(0, nl):
            if faces[l] == 'N':
                vect[:, l] = tempn[:, vecj[l], veci[l]]
            elif faces[l] == 'W':
                vect[:, l] = tempw[:, vecj[l], veci[l]]
            if use_orient:
                vect[:, l] *= secint.orient[l]

        # mask data where areavect is NaN
        # first, need to broadcast arrays for them to have the same
        # same size.
        # areavect, vect = np.broadcast_arrays(secint.areavect, vect)
        vect[:, np.isnan(secint.areavect[0, :])] = np.nan

        # if len(indnan[0] > 0):
        #     i, j = np.nonzero(np.isnan(areavect))
        #     vect[:, i, j] = np.nan
        # else:
        #     if filepres == np.float32:
        #         #print 'Masking vecv where vecv==0'
        #         #vecv[vecv == 0] = np.nan
        #         #areavect[np.isnan(vecv[0, :, :])] = np.nan
        #         #secint.areavect = areavect
        #     else:
        #         #print r'data precision is unsufficient: vecv will be inappropriately masked.'
        #         #print r'expected errors on indices calculations!'
        #         #vecv[vecv == 0] = np.nan

        # vect[np.ma.getmaskarray(vect)] = np.nan

        if not hasattr(secint, outname):
            setattr(secint, outname, vect)
        else:
            setattr(secint, outname, np.concatenate((getattr(secint, outname), vect), axis=0))
