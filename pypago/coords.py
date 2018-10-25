''' Module related to model coordinates '''

from __future__ import print_function
import numpy as np
import pylab as plt
from pypago.disp import PypagoErrors
import pypago.pyio
import pypago.misc
from pypago.sample_param import dictvname
try:
    from param import dictvname2
    dictvname.update(dictvname2)
except ImportError:
    pass

# we set the variables of the dictvname dict
# as global variables


class Coords(object):

    '''
    Initialisation of the :py:class:`Coords` class.
    It is initialised by providing a |netcdf| mesh file.
    The class reads and assigns the following attributes:
    - longitude
    - latitude
    - bathy (if exists, else None)
    - mask (if exists, else None)
    - dxt
    - dyt
    - dye (if exists, else None)
    - dxn (if exists, else None)
    - dzt (if exists, else None)

    :param str filename: Name of the |netcdf| mesh file

    .. note:: At this stage, all the file is read (no extraction
              over a subdomain)

    '''

    def __init__(self, filename):

        '''
        Initialisation of the Class. It initialises the filename attribute
        with the filename provided in argument. It first initialises all
        the arrays to None and then read them in the input NetCDF file.
        '''

        # Initialisation of the filename attribute
        self.filename = filename
        self.modelname = None

        # Initialisation of the attributes to None
        self.latt = None
        self.lont = None
        self.bathy = None
        self.mask = None
        self.dxt = None
        self.dyt = None
        self.dye = None
        self.dxn = None
        self.dyw = None
        self.dzt = None
        self.dzw = None
        self.dzn = None
        self.dze = None

        # reads the lon/lat variables from the file
        print('Reading longitude: variable %s' % dictvname['lon_varname'])
        self.lont = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['lon_varname']))

        print('Reading latitude: variable %s' % dictvname['lat_varname'])
        self.latt = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['lat_varname']))

        # Check whether the bathy variable exists.
        # If so, we define the bathy attribute. Else, remains None
        if 'bathy_varname' in dictvname:
            print('Reading bathymetry: variable %s' % dictvname['bathy_varname'])
            self.bathy = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['bathy_varname']))
        else:
            print('No bathymetry is read')

        # Check whether the tmask variable exists.
        # If so, we define the bathy attribute. Else, remains None
        if 'tmask_varname' in dictvname:
            print('Reading T-grid mask: variable %s' % dictvname['tmask_varname'])
            self.mask = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['tmask_varname']))

        # Cell meridional and zonal widths at center
        print('Reading T-grid zonal width: variable %s' % dictvname['dxt_varname'])
        self.dxt = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dxt_varname']))

        print('Reading T-grid meridional width: variable %s' % dictvname['dyt_varname'])
        self.dyt = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dyt_varname']))

        # Cell meridional width at eastern face of cell
        if 'dye_varname' in dictvname:
            print('Reading V-grid eastern meridional width: variable %s' % dictvname['dye_varname'])
            self.dye = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dye_varname']))

        # Cell zonal width at northern face of cell
        if 'dxn_varname' in dictvname:
            print('Reading U-grid northern meridional width: variable %s' % dictvname['dxn_varname'])
            self.dxn = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dxn_varname']))

        # Cell height at center
        if 'dzt_varname' in dictvname:
            print('Reading T-grid height: variable %s' % dictvname['dzt_varname'])
            self.dzt = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dzt_varname']))

    def __str__(self):

        ''' Redefinition of the string function '''

        output = 'Coords object:\n'
        output += '  -filename: %s\n' % self.filename
        output += '  -modelname: %s\n' % self.modelname
        attr_list = pypago.misc.extract_attrlist(self)
        attr_list = [attr for attr in attr_list if attr not in ['filename', 'modelname']]
        output += pypago.misc.extract_str(attr_list, self)

        return output

    def plot_mask(self, ax=None):

        ''' Contours the mask attribute '''

        if ax is None:
            ax = plt.gca()

        cs = ax.contour(self.mask, levels=[1 - np.spacing(1), 1], colors='k')
        ax.set_xlim(0, self.mask.shape[1]-1)
        ax.set_ylim(0, self.mask.shape[0]-1)

        return cs


class NemoCoords(Coords):

    """
    Coords class associated with the NEMO ocean model
    Inheritates from the :py:class:`pypago.coords.Coords`
    class::

        from pypago.coords import NemoCoords

        filename = 'nemo_mesh.nc'
        coords = NemoCoords(filename)

    :param str filename: Name of the |netcdf| mesh file

    """

    def __init__(self, filename):

        '''
        Initialisation of the NemoCoords class.
        '''

        # Initialisation through the mother class
        super(NemoCoords, self).__init__(filename)

        self.modelname = 'NEMO'

        # reading coordinates
        self.read_coord()

        # reading scale factors
        self.read_scalefactors()

    def read_coord(self):

        """
        Processes NEMO file coordinates (bathy, longitude, latitude, mask).
        - If the bathy array is None (i.e not read from input file)
        it is constructed from the mbathy and deptht arrays
        - Extracts the first level of the mask array
        """

        if self.bathy is None:
            # If bathymetry is not read directly from the file, it is
            # reconstructed from the 1D depth and mbathy arrays
            print('Reading 1D deptha array: variable %s' % dictvname['depth_varname'])
            deptht = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['depth_varname']))

            if deptht.ndim != 1:
                message = 'The deptht array must be 1D. Currently, ' + \
                          'ndim = %d. ' % deptht.ndim + 'This program will be stopped.'
                raise PypagoErrors(message)

            # extraction of the mbathy array
            print('Reading mbathy: variable %s' % dictvname['mbathy_varname'])
            mbathy = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['mbathy_varname']))
            mbathy = mbathy.astype(np.int)
            mbathy[mbathy < 0] = 0

            # reconstruction of the bathy from deptht and mbathy
            print('Reconstruction of bathy from mbathy and depth')
            self.bathy = deptht[mbathy]

        # Extraction of the surface (land sea) tmask array
        self.mask = np.squeeze(self.mask[0, :, :])

    def read_scalefactors(self):

        """
        Processes NEMO scale factors.
        - Initialises the dzw and dzn arrays. If dzt is 1D, then
        no partial step is assumed. If dzt is 2D, then 3D dzt, dzn
        and dze are reconstructed using the mbathy array. If dzt is
        3D then we assume that dze and dzn are in files
        - Masking the dzt, dzw, dzn arrays where the associated masks
        are 0 or masked
        - Reconstruction of dzw from dze

        """

        # If the e3t variable is 1D, no partial step.
        # e3t = e3u = e3v
        if self.dzt.ndim == 1:
            message = 'The dzt variable is 1D. '
            message += "Assumes no partial step. "
            message += "dzt = dzn = dzw = 1d array."
            nlat, nlon = self.lont.shape
            self.dzt = np.tile(self.dzt, (nlat, nlon, 1))
            self.dzt = np.transpose(self.dzt, (2, 0, 1))
            self.dzw = self.dzt.copy()  # no partial steps
            self.dzn = self.dzt.copy()  # no partial steps

        else:
            # if the dzt variable is 2D, then it gives the width of the last level
            # therefore, we reconstruct the e3t, e3u and e3v variables by using the
            # mbathy 2D variable
            if self.dzt.ndim == 2:
                print("The dzt variable is 2D")
                print("A 3D scale factor array at T, U and V " +
                      "points have been computed from the 2D scale factor")
                self.dzt = self.reconstruct_3d_e3t()

            else:
                print('Dzt is 3D. Model grid is in partial step')

            # putting dzt as NaN where 3D mask == 0
            self.dzt[np.squeeze(pypago.pyio.readnc(self.filename, dictvname['tmask_varname'])) == 0] = np.nan

            if ('dze_varname' in dictvname) and ('dzn_varname' in dictvname):
                print('Reading U-grid eastern height: variable %s' % dictvname['dze_varname'])
                self.dze = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dze_varname']))
                print('Reading V-grid northern height: variable %s' % dictvname['dzn_varname'])
                self.dzn = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dzn_varname']))
            else:
                message = 'The "dze_varname" and "dzn_varname" variables '
                message += 'are not defined. The dzn and dzw variables '
                message += 'are reconstructed from the dzt variable'
                print(message)
                self.reconstruct_3d_e3uv()

            if ('umask_varname' in dictvname) and ('vmask_varname' in dictvname):
                # if 3D umask and vmask variables are defined, then we mask the dzw and
                # dzn variables where masks are 0
                print('Reading U-grid mask: variable %s' % dictvname['umask_varname'])
                self.dze[np.squeeze(pypago.pyio.readnc(self.filename, dictvname['umask_varname'])) == 0] = np.nan
                print('Reading V-grid mask: variable %s' % dictvname['vmask_varname'])
                self.dzn[np.squeeze(pypago.pyio.readnc(self.filename, dictvname['vmask_varname'])) == 0] = np.nan

            # putting dzt as NaN where masked
            self.dzt[np.ma.getmaskarray(self.dzt) == 1] = np.nan
            self.dze[np.ma.getmaskarray(self.dze) == 1] = np.nan
            self.dzn[np.ma.getmaskarray(self.dzn) == 1] = np.nan

            # reconstructing the dzw array
            print('Reconstruction of U-grid western height from U-grid eastern height')
            self.dzw = np.concatenate((self.dze[:, :, -1:], self.dze[:, :, :-1]), axis=-1)

    def reconstruct_3d_e3t(self):

        """
        Creates, from 2D partial step value, a 3D cell width at T points.
        It is constructed by using the "mbathy" variable (index of the last
        non-ocean point) and the 1D constant scale factor

        :return: A numpy array containing the cell width at T, U and V points
        :rtype: numpy.array
        """

        dzt2d = self.dzt.copy()

        # extraction of the mbathy variable
        mbathy = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['mbathy_varname']))
        mbathy = mbathy.astype(np.int)

        # extraction of 1D scale factor
        e3t_0 = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dzt1d_varname']))

        # recovering the dimensions of input array
        nlat, nlon = mbathy.shape

        # initialisation of e3t, e3u and e3v arrays
        # by using the 1D constant scale factors
        self.dzt = np.tile(e3t_0, (nlat, nlon, 1))
        self.dzt = np.transpose(self.dzt, (2, 0, 1))

        for ilat in xrange(0, nlat):
            for ilon in xrange(0, nlon):

                # mbathy value at the ilat/ilon point
                mbathytemp = mbathy[ilat, ilon]

                # if 0, then land everywhere
                if mbathytemp > 0:
                    self.dzt[mbathytemp-1, ilat, ilon] = dzt2d[ilat, ilon]

    def reconstruct_3d_e3uv(self):

        """
        Creates, from the dzt variable, a 3D cell width at northern and
        southern grid faces.

        :return: A tuple containing the cell width at U and V points
        :rtype: numpy.array
        """

        self.dze = self.dzt.copy()
        self.dzn = self.dzt.copy()

        # recovering the dimensions of input array
        nlat, nlon = self.dzt.shape[1], self.dzt.shape[2]

        for ilat in xrange(0, nlat-1):
            for ilon in xrange(0, nlon):
                self.dzn[:, ilat, ilon] = np.min((self.dzt[:, ilat, ilon], self.dzt[:, ilat+1, ilon]), axis=0)

        for ilat in xrange(0, nlat):
            for ilon in xrange(0, nlon-1):
                self.dze[:, ilat, ilon] = np.min((self.dzt[:, ilat, ilon], self.dzt[:, ilat, ilon+1]), axis=0)


class GfdlCoords(Coords):

    """
    Coords class associated with the GFDL ocean model
    Inheritates from the :py:class:`pypago.coords.Coords`
    class::

        from pypago.coords import GfdlCoords

        filename = 'gfdl_mesh.nc'
        coords = GfdlCoords(filename)

    :param str filename: Name of the mesh file

    """

    def __init__(self, filename):

        ''' Initialisation of the GfdlCoords class '''

        super(GfdlCoords, self).__init__(filename)

        self.modelname = 'GFDL'
        self.dzc = None

        self.read_coord()
        self.read_scalefactors()

    def read_coord(self):

        """
        Processes GFDL file coordinates (bathy, longitude, latitude, mask).
        - If sets to 0 all the masked/NaN values within the tmask array
        """

        # Putting 0 when NaN or masked
        self.mask[self.mask != self.mask] = 0
        self.mask[np.ma.getmaskarray(self.mask) == 1] = 0

    def read_scalefactors(self):

        """
        Processes GFDL scale factors.
        - Extracts the dzc array
        - Set the dzt, dzc arrays to NaN where 0
        - Shift from west to east of the dzc array
        (depending on user input)

        """

        try:
            from param import gfdl_default
        except ImportError:
            message = 'Either the param.py file is missing \n'
            message += 'or the "gfdl_default" variable is missing'
            error = PypagoErrors(message)
            raise error

        # at north WEST corner of cell
        self.dzc = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dzc_varname']))

        # masking of the dzt and dzc variables
        self.dzt[self.dzt == 0] = np.nan
        self.dzc[self.dzc == 0] = np.nan

        if not gfdl_default:
            message = 'gfdl_default is True'
            message += 'hence we assume that corner '
            message += 'point is located at the north east corner of cell'
            print(message)
            self.dzc = np.concatenate((self.dzc[:, :, -1:], self.dzc[:, :, :-1]), axis=-1)


class CcsmCoords(Coords):

    """
    Coords class associated with the CCSM ocean model
    Inheritates from the :py:class:`pypago.coords.Coords`
    class::

        from pypago.coords import CcsmCoords

        filename = 'ccsm_mesh.nc'
        coords = CcsmCoords(filename)

    :param str filename: Name of the |netcdf| mesh file

    """

    def __init__(self, filename):

        ''' Initialisation of the CcsmCoords class '''

        super(CcsmCoords, self).__init__(filename)
        self.modelname = 'CCSM'

        self.read_coord()
        self.read_scalefactors()

    def read_coord(self):

        """
        Processes CCSM file coordinates (bathy, longitude, latitude, mask).
        - If converts the bathymetry array from centimeters to meters
        """

        # Conversion of bathymetry from centimeters to meters
        self.bathy = self.bathy * 1e-2

    def read_scalefactors(self):

        """
        Processes CCSM scale factors.
        - Converts the scale factors from cm into m
        - Copy the dzt into dzw and dzn (no partial steps)
        """

        # conversion from cm to m
        self.dxt = self.dxt * 1e-2
        self.dyt = self.dyt * 1e-2
        self.dxe = self.dxe * 1e-2
        self.dyn = self.dyn * 1e-2
        self.dzt = self.dzt * 1e-2

        self.dzw = self.dzt.copy()
        self.dzn = self.dzt.copy()


class MpioCoords(Coords):

    """
    Coords class associated with the MPIO ocean model
    Inheritates from the :py:class:`pypago.coords.Coords`
    class::

        from pypago.coords import MpioCoords

        filename = 'mpio_mesh.nc'
        coords = MpioCoords(filename)

    :param str filename: Name of the |netcdf| mesh file

    """

    def __init__(self, filename):

        ''' Initialisation of the MpioCoords class '''

        super(MpioCoords, self).__init__(filename)
        self.modelname = 'MPIO'
        self.read_coord()
        self.read_scalefactors()

    def read_coord(self):

        """
        Processes MPIO file coordinates (bathy, longitude, latitude, mask).
        - If flips all the arrays along their rightmost dimension (latitude)
        so that north is on the top of the figures (using flipud)
        - Extracts the first level of the tmask array
        - Setting the tmask to 0 where Nan/masked
        """

        # extracting the tmask at the first level
        self.mask = self.mask[0, :, :]

        # flipping up-down the arrays
        self.lont = np.flipud(self.lont)
        self.latt = np.flipud(self.latt)
        self.bathy = np.flipud(self.bathy)
        self.mask = np.flipud(self.mask)

        # where tmask is NaN/masked, we set it to 0
        self.mask[np.ma.getmaskarray(self.mask) == 1] = 0
        self.mask[self.mask != self.mask] = 0

    def read_scalefactors(self):

        """
        Processes MPIO scale factors.
        - Flipping the dxt, dyt, dxn and dye arrays along
        the first dimension (flipud)
        - Flipping the dzt array along the 2nd dimension
        - Creation of dzw by multiplication of dzt with the flipped
        amsuo variable
        - Creation of dzn by multiplication of dzt with the flipped
        amsue variable
        """

        # flipping the dimensions of the dxt/dyt arrays to
        # make north on top
        self.dxt = np.flipud(self.dxt)
        self.dyt = np.flipud(self.dyt)
        self.dxn = np.flipud(self.dxn)
        self.dye = np.flipud(self.dye)

        # flipping dzt along the 2nd dimension
        # flipud doesnt work here
        self.dzt = self.dzt[:, ::-1, :]
        self.dzw = self.dzt * np.squeeze(pypago.pyio.readnc(self.filename, dictvname['amsuo_varname']))[:, ::-1, :]
        self.dzn = self.dzt * np.squeeze(pypago.pyio.readnc(self.filename, dictvname['amsue_varname']))[:, ::-1, :]


class HycoCoords(Coords):

    """
    Coords class associated with the Hyco ocean model
    Inheritates from the :py:class:`pypago.coords.Coords`
    class::

        from pypago.coords import HycoCoords

        filename = 'Hyco_mesh.nc'
        coords = HycoCoords(filename)

    :param str filename: Name of the |netcdf| mesh file

    """

    def __init__(self, filename):

        ''' Initialisation of the HycoCoords class '''

        super(HycoCoords, self).__init__(filename)

        self.modelname = 'HYCO'

        self.read_coord()
        self.read_scalefactors()

    def read_coord(self):

        """
        Processes HYCO file coordinates (bathy, longitude, latitude, mask).
        - Reconstructs the mask array from the bathy (where bathy>0, mask=1)
        """

        # mask in Hyco does not exist. We reconstruct it
        self.mask = np.zeros(self.bathy)
        self.mask[self.bathy > 0] = 1

    def read_scalefactors(self):

        """
        Processes HYCO scale factors.
        - Reconstruction of dxn from dxt by adding a line of
        NaN at the bottom of dxt
        - Extraction of the dyw array
        - Extraction of the layer index in order to recover
        the number of vertical levels
        - Creation of dzt, dzn and dze as 3D arrays of ones.
        """

        # reconstruction of dxn using the scqx array
        # initially read as the dxn array
        # adding NaNs on the rightmost column
        self.dxn = pypago.pyio.readnc(self.filename, dictvname['dxn_varname'])
        temp = np.nan*np.ones(self.dxn[-1:, :].shape)
        self.dxn = np.concatenate((self.dxn[1:, :], temp), axis=0)

        self.dyw = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['dyw_varname']))

        layer_index = np.squeeze(pypago.pyio.readnc(self.filename, dictvname['layind_varname']))
        self.dzt = np.ones((len(layer_index), self.dxt.shape[0], self.dxt.shape[1]))
        self.dzn = self.dzt.copy()
        self.dzw = self.dzt.copy()


class RomsCoords(Coords):

    """
    Coords class associated with the Roms ocean model
    Inheritates from the :py:class:`pypago.coords.Coords`
    class::

        from pypago.coords import RomsCoords

        filename = 'Roms_mesh.nc'
        coords = RomsCoords(filename)

    :param str filename: Name of the |netcdf| mesh file

    """

    def __init__(self, filename):

        ''' Initialisation of the class '''

        super(RomsCoords, self).__init__(filename)
        self.modelname = 'ROMS'
        self.read_scalefactors()
        self.dys = None

    def read_scalefactors(self):

        """
        Processes ROMS scale factors.
        - Taking the inverse of dxt, dyt, dye and dxn
        - Adding a line of NaNs at the bottom of the dzt array
        - Extracting the dze and dzn arrays
        - Creating the dzw array from the dze array by adding a
        layer of NaNs on the easternmost face of the "cube"
        """

        # taking the inverse of scale factors
        self.dxt = 1./self.dxt
        self.dyt = 1./self.dyt

        # Reading dxs on southern faces, and reconstruct dxn by using dxn
        print('Reconstructing the V-grid northern width using width of T grid')
        self.dxn = np.ones(self.dxt.shape, dtype=np.float) * np.nan   # barrier.n: initialize dxn variable as NaNs
        self.dxn[:-1, :] = 0.5 * (self.dxt[1:, :] + self.dxt[:-1, :])    # computes the dxn values as the mean of dxt

        # reading dy on western faces
        print('Reconstructing the U-grid western width using width of T grid')
        self.dyw = np.ones(self.dyt.shape, dtype=np.float) * np.nan
        self.dyw[:, :-1] = 0.5 * (self.dyt[:, 1:] + self.dyt[:, :-1])

        if self.dzt is not None:
            # if the input mask contains the z_rho array, then it
            # is assumed that it also contains the z_u and z_v variables

            # extracting first time step for dzt
            self.dzt = self.dzt[0, :, :, :]

            # reads the dze 4d variables with start=[1,1,1,1], end=[1, -1, -1, -1]
            print('Reading dz variable on western faces: variable %s' % dictvname['dzw_varname'])
            self.dzw = pypago.pyio.readnc(self.filename, dictvname['dzw_varname'], 4*[1], [1]+3*[-1])

            print('Reading dz variable on southern faces: variable %s' % dictvname['dzs_varname'])
            dzs = pypago.pyio.readnc(self.filename, dictvname['dzs_varname'], 4*[1], [1]+3*[-1])
            print('Reconstruction of dz on northern faces from dz on southern faces')
            self.dzn = np.ones(self.dzt.shape, dtype=np.float) * np.nan
            self.dzn[:, :-1, :] = dzs[:, :, :]

        else:
            # if the z_rho variables, the same thing as in the GFDL
            # model is performed, i.e dzt = dzw = dzn = 1
            print('Initialisation of dzt, dzw and dzn as ')
            print('(roms_nsigma, nlat, lon) arrays of ones')
            self.dzt = np.ones([dictvname['roms_nsigma']] + list(self.lont.shape))
            self.dzw = self.dzt
            self.dzn = self.dzt


class OfamCoords(Coords):

    """
    Coords class associated with the Ofam ocean model
    Inheritates from the :py:class:`pypago.coords.Coords`
    class::

        from pypago.coords import OfamCoords

        filename = 'Ofam_mesh.nc'
        coords = OfamCoords(filename)

    :param str filename: Name of the |netcdf| mesh file

    """

    def __init__(self, filename):

        ''' Initialisation of the OfamCoords class '''

        super(OfamCoords, self).__init__(filename)
        self.modelname = 'OFAM'

        self.read_scalefactors()

    def read_scalefactors(self):

        """
        Processes OFAM scale factors.
        - Extraction of the dyw array
        - Extraction of the dzb array (width at z-levels)
        - Reconstruction of dzt array from dzb
        - Setting dzw and dzn as equal to dzt (no partial steps)
        """

        # extracting the dyw variable
        self.dyw = pypago.pyio.readnc(self.filename, dictvname['dyw_varname'])

        # extracting the dzb variable
        dzb = pypago.pyio.readnc(self.filename, dictvname['dzb_varname'])

        # calculating the dzt array as the difference between 2 cons levels
        self.dzt = np.diff(np.concatenate(([0], dzb)))

        # copying dzt into dzw and dzn (assumes no partial step)
        self.dzw = self.dzt.copy()
        self.dzn = self.dzt.copy()


class MicoCoords(Coords):

    """
    Coords class associated with the Mico ocean model
    Inheritates from the :py:class:`pypago.coords.Coords`
    class::

        from pypago.coords import MicoCoords

        filename = 'Mico_mesh.nc'
        coords = MicoCoords(filename)

    :param str filename: Name of the |netcdf| mesh file

    """

    def __init__(self, filename):

        ''' Initialisation of the MicoCoords class '''

        super(MicoCoords, self).__init__(filename)
        self.modelname = 'MICO'

    def read_scalefactors(self):

        """
        Processes MICO scale factors.
        - Extract dx at southern faces
        - Reconstruct dx at northern faces (dxn)
        - Extraction of dyw
        - Creation of vertical scale factors:
        a) In isopicnic coordinates: asking for
        the number of vertical levels, setting
        dzt, dzw and dzn as 1D array of ones
        b) In interpolated z-coordinates: opening of
        a file with the z-bounds, reconstruction
        of dzt. And copy of the dzt values to
        dzw and dzn (no scale factors)
        """

        try:
            from param import micom_isopycnic, micom_nlevels
        except ImportError:
            message = 'The "param.py" file is either missing \n'
            message += 'or does not contain the micom_isopycnic \n'
            message += 'and micom_nlevels variables.'
            raise PypagoErrors(message)

        # reading cell width at southern face
        dxs = pypago.pyio.readnc(self.filename, dictvname['dxs_varname'])

        # moving cell width at northern face
        self.dxn = np.zeros(dxs.shape)
        self.dxn[:-1, :] = dxs[1:, :]

        # reading the dyw variable
        self.dyw = pypago.pyio.readnc(self.filename, dictvname['dyw_varname'])

        if micom_isopycnic:
            self.dzt = np.ones(micom_nlevels)
            self.dzn = np.ones(micom_nlevels)
            self.dzw = np.ones(micom_nlevels)

        else:
            zbounds = pypago.pyio.readnc(dictvname['zfile'], dictvname['dzt_varname'])
            self.dzt = zbounds[:, 1] - zbounds[:, 0]
            self.dzn = self.dzt.copy()
            self.dzw = self.dzt.copy()


def create_coord(modelname, filename):

    """
    Returns a Coord object, depending on the modelname

    :param str modelname: Name of the model
    :param str filename: Name of the |netcdf| mesh file

    """

    # List of possible models
    modellist = ["GFDL", "CCSM", "NEMO", "MPIO", "HYCO", "ROMS", "OFAM", "MICO"]

    # checking that the model name exists.
    if modelname not in modellist:
        message = "The %s model name is unrecognized. " % modelname
        message += "Possible values are: " + "/".join(modellist) + ". "
        message += "The program will stop"
        raise PypagoErrors(message)

    # returning the appropriate Coord class
    # depending on the model name
    if modelname == "NEMO":
        output = NemoCoords(filename)
    elif modelname == "GFDL":
        output = GfdlCoords(filename)
    elif modelname == "CCSM":
        output = CcsmCoords(filename)
    elif modelname == "MPIO":
        output = MpioCoords(filename)
    elif modelname == "HYCO":
        output = HycoCoords(filename)
    elif modelname == "ROMS":
        output = RomsCoords(filename)
    elif modelname == "OFAM":
        output = OfamCoords(filename)
    elif modelname == "MICO":
        output = MicoCoords(filename)

    return output
