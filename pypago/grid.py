
''' Class that handles grid scale factors for a model configuration '''

from __future__ import print_function
import numpy as np
import pylab as plt
from pypago.disp import PypagoErrors
import pypago.disp
import pypago.coords


class Grid(object):

    """ Class that handles grid extraction on a specific domain

    :param int jmin: Index of the southernmost point of the subdomain

    :param int jmax: Index of the northernmost point of the subdomain

    :param int imin: Index of the westernmost point of the subdomain

    :param int imax: Index of the easternnmost point of the subdomain

    :param :py:class:`pypago.coords.Coords` coord: Coordinate object
     associated with the grid file (longitude, latitude, mask of the
     entire domain, used in the file extraction)

    """

    def __init__(self, coord, jmin, jmax,
                 imin, imax):

        ''' Initialisation of the grid class '''

        # here are the coords and scale factors
        # extracted on the domain grid.
        self.latt = None
        self.lont = None
        self.mask = None
        self.dyw = None
        self.bathy = None
        self.dxt = None
        self.dyt = None
        self.dyte = None
        self.dxn = None
        self.dzt = None
        self.dzw = None
        self.dzn = None
        self.areaw = None
        self.arean = None
        self.volume = None
        self.surface = None
        self.jmin = None
        self.jmax = None
        self.imin = None
        self.imax = None

        self.filename = None
        self.modelname = None
        self.nlon = None

        if coord is not None:

            for arg in [jmin, jmax, imin, imax]:
                if arg is None:
                    message = 'Either jmin, jmax, imin or imax '
                    message += 'is None.'
                    raise PypagoErrors(message)

            nlat, nlon = coord.lont.shape

            # Checking that the jmin/jmax/imin/imin are within the bounds
            # if not, correct their values
            if (imin < 0) | (imin > nlon - 1):
                print('The imin argument must be ' +
                      'between %d and %d. Currently, %d' % (0, nlon-1, imin))
                print('imin has been set to %d' % (0))
                imin = 0
            if (imax < 0) | (imax > nlon - 1):
                print('The imax argument must be ' +
                      'between %d and %d. Currently, %d' % (0, nlon-1, imax))
                print('imax has been set to %d' % (nlon-1))
                imax = nlon-1
            if (jmin < 1) | (jmin > nlat - 2):
                print('The jmin argument must be ' +
                      'between %d and %d. Currently, %d' % (1, nlat-2, jmin))
                print('jmin has been set to %d' % (1))
                jmin = 1
            if (jmax < 1) | (jmax > nlat - 2):
                print('The jmax argument must be ' +
                      'between %d and %d. Currently, %d' % (1, nlat-2, jmax))
                print('jmax has been set to %d' % (nlat-2))
                jmax = nlat-2

            if jmax < jmin:
                print("Currently, jmax<jmin. " +
                      "The variables have been switched")
                jmax, jmin = jmin, jmax

            self.filename = coord.filename
            self.modelname = coord.modelname
            self.nlon = nlon

            self.jmin = jmin
            self.jmax = jmax
            self.imin = imin
            self.imax = imax

            self.extract_all_var(coord)

    def plot_dom(self, ax=None):

        """
        Draws the domain defined by the jmax, jmin,
        imin, and imax attributes.

        :param matplotlib.axes.Axes ax: Axis on which to draw
         the figure. If None, draws on current axis.

        """

        if ax is None:
            ax = plt.gca()

        if self.imin < self.imax:
            lontp = [self.imin, self.imax, self.imax, self.imin, self.imin]
            lattp = [self.jmin, self.jmin, self.jmax, self.jmax, self.jmin]
            ax.plot(lontp, lattp, color='r')
        else:
            lontp = [self.nlon, self.imin, self.imin, self.nlon]
            lattp = [self.jmin, self.jmin, self.jmax, self.jmax]
            ax.plot(lontp, lattp, color='r')

            lontp = [0, self.imax, self.imax, 0]
            lattp = [self.jmin, self.jmin, self.jmax, self.jmax]
            ax.plot(lontp, lattp, color='r')

    def __str__(self):

        output = 'Model grid for the %s model:\n' % self.modelname
        output += '    - mesh file: %s\n' % self.filename
        output += '    - jmin: %d\n' % self.jmin
        output += '    - jmax: %d\n' % self.jmax
        output += '    - imin: %d\n' % self.imin
        output += '    - imax: %d\n' % self.imax
        output += '    - nlat: %d\n' % self.mask.shape[0]
        output += '    - nlon: %d\n' % self.mask.shape[1]
        if self.volume is not None:
            output += '    - nz: %d\n' % self.volume.shape[0]
        return output

    def compute_areas(self):

        """
        Computes the cell area at the northern and
        western faces, and the volume and surface
        of the T grid cells.
        """

        nz = self.dzt.shape[0]
        self.areaw = self.dzw * np.tile(self.dyw, (nz, 1, 1))
        self.arean = self.dzn * np.tile(self.dxn, (nz, 1, 1))

        self.areaw[self.areaw == 0] = np.nan
        self.areaw[np.ma.getmaskarray(self.areaw) == 1] = np.nan

        self.arean[self.arean == 0] = np.nan
        self.arean[np.ma.getmaskarray(self.arean) == 1] = np.nan

        self.surface = self.dxt * self.dyt
        self.volume = self.dzt * np.tile(self.surface, (nz, 1, 1))

    def extract_2d_var(self, varin):

        """
        Extraction of a 2D variable on the subdomain defined
        by the four domain limits. It provides the possibility
        to read a variable "discontinuously"
        (if `imin` > `imax`).

        :param numpy.array varin: Input 2D variable
        :return: The input array extracted on the specified subdomain
        :rtype: numpy.array

        """

        if self.imin > self.imax:
            out1 = varin[self.jmin:self.jmax+1, self.imin:]
            out2 = varin[self.jmin:self.jmax+1, :self.imax+1]
            varout = np.concatenate((out1, out2), axis=-1)
        else:
            varout = varin[self.jmin:self.jmax+1, self.imin:self.imax+1]

        return varout

    def extract_3d_var(self, varin):

        """
        Extraction of a 2D variable on the subdomain defined
        by the four domain limits. It provides the possibility
        to read a variable "discontinuously"
        (if `imin` > `imax`).

        :param numpy.array varin: Input 3D variable
        :return: The input array extracted on the specified subdomain
        :rtype: numpy.array

        """

        if self.imin > self.imax:
            out1 = varin[:, self.jmin:self.jmax+1, self.imin:]
            out2 = varin[:, self.jmin:self.jmax+1, :self.imax+1]
            varout = np.concatenate((out1, out2), axis=-1)

        else:
            varout = varin[:, self.jmin:self.jmax+1, self.imin:self.imax+1]

        return varout

    def create_3d_vert_scalefact(self):

        """
        If the vertical scale factor dzt is 1D,
        it reconstructs a 3D variable by repeating
        the 1D variable along the lon and lat dimensions
        """

        ny, nx = self.dxt.shape
        self.dzt = np.transpose(np.tile(self.dzt, [ny, nx, 1]), [2, 0, 1])
        self.dzw = self.dzt.copy()
        self.dzn = self.dzt.copy()

    def extract_all_var(self, coord):

        """
        Extracts all the variables (coordinates,
        scale factors, etc) on the sub-domain
        defined by the jmin, jmax, imin and imax variables.
        """

        print("Extraction of lont on the domain")
        self.lont = self.extract_2d_var(coord.lont)
        print("Extraction of latt on the domain")
        self.latt = self.extract_2d_var(coord.latt)
        print("Extraction of bathy on the domain")
        self.bathy = self.extract_2d_var(coord.bathy)
        print("Extraction of mask on the domain")
        self.mask = self.extract_2d_var(coord.mask)
        print("Extraction of dxt on the domain")
        self.dxt = self.extract_2d_var(coord.dxt)
        print("Extraction of dyt on the domain")
        self.dyt = self.extract_2d_var(coord.dyt)
        print("Extraction of dxn on the domain")
        self.dxn = self.extract_2d_var(coord.dxn)

        if coord.dyw is not None:
            # If the model contains the dyw variable (HYCO, OFAM, MICO models for instance)
            print("Extraction of dyw on the domain")
            self.dyw = self.extract_2d_var(coord.dyw)

        else:
            # If the model does not contain the dyw variable, we reconstruct it from the
            # dye variable
            print("Reconstruction of the dyw variable from dye")
            if self.imin < self.imax:
                if self.imin > 0:
                    self.dyw = coord.dye[self.jmin:self.jmax+1, self.imin-1:self.imax]
                else:
                    self.dyw = np.concatenate((coord.dye[self.jmin:self.jmax+1, -1:],
                                               coord.dye[self.jmin:self.jmax+1, self.imin:self.imax]),
                                              axis=-1)
            else:   # self.imax < self.imin
                if self.imax > 0:
                    self.dyw = np.concatenate((coord.dye[self.jmin:self.jmax+1, self.imin-1:],
                                               coord.dye[self.jmin:self.jmax+1, :self.imax]),
                                              axis=-1)
                else:
                    self.dyw = coord.dye[self.jmin:self.jmax, self.imin-1:]

    def plot_mask(self, ax=None):

        ''' Contours the mask attribute '''

        if ax is None:
            ax = plt.gca()

        cs = ax.contour(self.mask, levels=[1 - np.spacing(1), 1], colors='k')
        ax.set_xlim(0, self.mask.shape[1]-1)
        ax.set_ylim(0, self.mask.shape[0]-1)

        return cs


class CommonGrid(Grid):

    ''' Grid class associated with all models but GFDL '''

    def __init__(self, coord, jmin, jmax, imin, imax):

        ''' Initialisation of the common grid class '''

        # initialisation of the attributes to None
        # using the mother initialisation
        super(CommonGrid, self).__init__(coord, jmin, jmax,
                                         imin, imax)

        if coord.dzt.ndim > 2:
            self.dzt = self.extract_3d_var(coord.dzt)
            self.dzw = self.extract_3d_var(coord.dzw)
            self.dzn = self.extract_3d_var(coord.dzn)

        else:
            # if dzt is 1D, we create 3D scale factors
            # from 1D scale factors
            self.create_3d_vert_scalefact()

        self.compute_areas()


class GfdlGrid(Grid):

    ''' GFDL grid class '''

    def __init__(self, coord, jmin, jmax, imin, imax):

        ''' Initialisation of the GFDL grid class '''

        # initialisation of the attributes to None
        # using the mother initialisation
        super(GfdlGrid, self).__init__(coord, jmin, jmax,
                                       imin, imax)

        self.dzc = None

        self.extract_3d_var_gfdl(coord)
        self.compute_areas()

    def extract_3d_var_gfdl(self, coord):

        """
        Creates the vertical scale factors for the
        GFDL model.
        """

        if coord.dzt.ndims > 2:

            if self.imax < self.imin:

                dzc = np.concatenate((coord.dzc[:, self.jmin-1:self.jmax+1, self.imin:],
                                      coord.dzc[:, self.jmin-1:self.jmax+1, :self.imax+2]),
                                     axis=-1)

                temp = np.array([dzc[:, 1:, :-1],
                                 dzc[:, :-1, :-1]])
                temp = np.ma.masked_where(temp != temp, temp)
                self.dzw = np.mean(temp, axis=0)

                temp = np.array([dzc[:, 1:, 1:],
                                 dzc[:, 1:, :-1]])
                temp = np.ma.masked_where(temp != temp, temp)
                self.dzn = np.mean(temp, axis=0)
                self.dzt = self.extract_3d_var(coord.dzt)

            else:   # imin < imax

                if self.imax < self.mask.shape[-1]:
                    dzc = coord.dzc[:, self.jmin-1:self.jmax+1, self.imin:self.imax+2]
                else:
                    dzc = np.concatenate((coord.dzc[:, self.jmin-1:self.jmax+1, self.imin:self.imax+1],
                                          coord.dzc[:, self.jmin-1:self.jmax+1, 0:1]), axis=-1)

                temp = np.array([dzc[:, 1:, :-1],
                                 dzc[:, :-1, :-1]])
                temp = np.ma.masked_where(temp != temp, temp)
                self.dzw = np.mean(temp, axis=0)

                temp = np.array([dzc[:, 1:, 1:],
                                 dzc[:, 1:, :-1]])
                temp = np.ma.masked_where(temp != temp, temp)
                self.dzn = np.mean(temp, axis=0)

                self.dzt = self.extract_3d_var(coord.dzt)

        else:
            # if dzt is 1D, we create 3D scale factors
            # from 1D scale factors
            self.create_3d_vert_scalefact()


def create_grid(coord, jmin=None, jmax=None, imin=None, imax=None):

    """
    Returns a Grid instance associated with the input coord argument. If it
    is a `pypago.coord.GfdlCoords` object, then it returns a
    :py:class:`pypago.grid.GfdlGrid` object. Else, it returns a
    :py:class:`pypago.grid.CommonGrid` object.

    :param int jmin: Index of the southernmost point of the subdomain

    :param int jmax: Index of the northernmost point of the subdomain

    :param int imin: Index of the westernmost point of the subdomain

    :param int imax: Index of the easternnmost point of the subdomain

    :param :py:class:`pypago.coords.Coords` coord: Coordinate object
     associated with the grid file (longitude, latitude, mask of the
     entire domain, used in the file extraction)

    """

    if imin is None:
        print('The imin argument is None. Set to 0')
        imin = 0
    if imax is None:
        print('The imax argument is None. Set to %d' % (coord.lont.shape[1] - 1))
        imax = coord.lont.shape[1] - 1
    if jmin is None:
        jmin = 1
        print('The jmin argument is None. Set to 1')
    if jmax is None:
        print('The jmax argument is None. Set to %d' % (coord.lont.shape[0] - 2))
        jmax = coord.lont.shape[0] - 2

    if isinstance(coord, pypago.coords.GfdlCoords):
        output = GfdlGrid(coord, jmin, jmax, imin, imax)
    else:
        output = CommonGrid(coord, jmin, jmax, imin, imax)

    return output
