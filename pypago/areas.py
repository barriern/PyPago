# -*- coding: utf-8 -*-

"""
This function allows to extract spatial domains, by giving its section
boundaries.

"""

from __future__ import print_function
import numpy as np
import pypago.plot
import pypago.pyio
import pypago.misc
from pypago.disp import PypagoErrors


class Areas(object):

    """ Domain area class """

    def __init__(self, grid, name, i, j, secnames=None, signs=None):

        self.name = name
        self.i = i
        self.j = j

        self.mask = np.zeros(grid.mask.shape)
        self.mask[i, j] = 1

        self.volume = grid.volume[:, i, j]
        self.surface = grid.surface[i, j]

        self.secnames = secnames
        self.signs = signs

        self.jmin = grid.jmin
        self.jmax = grid.jmax
        self.imin = grid.imin
        self.imax = grid.imax
        self.nlon = grid.nlon
        self.modelname = grid.modelname

    def __str__(self):

        output = 'Domain area, %s model:\n' % self.modelname

        attrnames = pypago.misc.extract_attrlist(self)
        output += pypago.misc.extract_str(attrnames, self)

        return output

    def finalplot(self, grid, gridsec=None, ax=None):

        """
        Draws the areas mask.
        It draws the model mask and the model sections,
        and it fills the mask of the area

        :param pypago.grid.Grid grid: |pypago| object that contains the model grid
        :param pypago.sections.GridSection gridsec: |pypago| object that contains the model
           sections
        :param numpy.array mask: mask of the area
        """

        pypago.plot.plot_dom_mask(grid, gridsec, self.mask, ax=ax)


def extract_dom_from_pol(grid, lonpol, latpol):

    """
    Extracts the i, j indexes of the grid domain
    by provinding the longitudes and latitudes of
    a polygon.

    :param pypago.grid.Grid grid: Input grid
    :param numpy.array lonpol: Polygon longitude
    :param numpy.array lonpol: Polygon latitude

    :return: A tuple containing the i, j indexes
     of the domain point within the domain.

    :rtype: tuple

    """

    from matplotlib.path import Path

    # convert list inputs into arrays
    lonpol = np.array(lonpol)
    latpol = np.array(latpol)

    # add cyclic values (if not)
    if lonpol[-1] != lonpol[0]:
        lonpol = np.append(lonpol, lonpol[-1:])
    if latpol[-1] != latpol[0]:
        latpol = np.append(latpol, latpol[-1:])

    if len(lonpol) != len(latpol):
        message = 'The lonpol and latpol must have the same length.\n'
        message += 'Currently:\n'
        message += ' -len(lonpol) = %3d\n' % len(lonpol)
        message += ' -len(latpol) = %3d\n' % len(latpol)
        message += 'This program will be stopped.'
        raise PypagoErrors(message)

    nlat, nlon = grid.latt.shape

    # creation the input of the path.Path command:
    # [(x1, y1), (x2, y2), (x3, y3)]
    path_input = [(lontemp, lattemp) for lontemp, lattemp in zip(lonpol, latpol)]
    pathobj = Path(path_input)

    mask = np.array([pathobj.contains_point((lontemp, lattemp)) for lontemp, lattemp in zip(np.ravel(grid.lont), np.ravel(grid.latt))])

    mask = np.reshape(mask, (nlat, nlon))

    i, j = np.nonzero(mask == 1)

    return i, j
