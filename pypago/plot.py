
"""
Module that contains various functions dedicated to plotting using |pypago|
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from pypago.misc import PypagoErrors, extract_mask, make_percentile_cmap


def plot_dom_mask(grid, gridsec=None, mask=None, ax=None):

    """
    Draws the final plot.
    It draws the model mask and the model sections,
    and it fills the mask of the area

    :param pago_obj grid: |pypago| object that contains the model grid
    :param pago_obj gridsec: |pypago| object that contains the model
       sections
    :param numpy.array mask: mask of the area
    """

    if ax is None:
        ax = plt.gca()

    cmap2 = ListedColormap(['none', 'white'])
    cmap = ListedColormap(['black', 'gray', 'white'])

    if mask is not None:
        plt.imshow(mask, interpolation='none', cmap=cmap)
    else:
        plt.imshow(grid.mask, interpolation='none', cmap=cmap)
    if gridsec is not None:
        for sec in gridsec:
            sec.plotsecfaces(ax)
    plt.title('sections in grid')
    plt.xlabel('grid points')
    plt.ylabel('grid points')
    plt.xlim(0, grid.mask.shape[1])
    plt.ylim(0, grid.mask.shape[0])


def _data_prep(secstruct, vartoplot, output_form, istracer=None, itime=None):

    """ Data extraction, preparation and masking prior to plotting """

    if not hasattr(secstruct, vartoplot):
        message = 'The "%s" variables is not an attribute '
        message += 'of the input section.'
        raise PypagoErrors(message)

    secdata = getattr(secstruct, vartoplot)

    if secdata.ndim == 3:
        if itime is None:
            message = 'The time average of the %s ' % vartoplot
            message += 'is plotted.'
            print(message)
            secdata = np.mean(secdata, axis=0)
        else:
            try:
                secdata = secdata[itime, :, :]
            except:
                message = 'The data could not be extracted. '
                message += 'Check that the time index is ok. '
                message += 'This program will be stopped'
                raise PypagoErrors(message)

    atracer, zvect, lvect  = preplot(secstruct, secdata, istracer)

    return (atracer, zvect, lvect)


def pcolplot(secstruct, vartoplot, istracer, itime=None, ax=None):

    """ Pcolor of a variable contained in a gridded section.

    :param pypago.sections.GridSection secstruct: Gridded section
    :param str vartoplot: Variable to plot
    :param bool istracer: True if tracer field (i.e. temperature),
     False if velocity field (or transport field).
    :param int itime: Time index to extract. If None,
     temporal mean is computed
    :param ax ax: Axis. If None, gca()

    :return: cs, cb

    """

    atracer_pcol, zvect_pcol, lvect_pcol = _data_prep(secstruct, vartoplot, 'pcolor', istracer, itime)

    zmin, zmax = zvect_pcol.min(), zvect_pcol.max()
    lmin, lmax = lvect_pcol.min(), lvect_pcol.max()

    if ax is None:
        ax = plt.gca()

    cs = plt.tripcolor(lvect_pcol, zvect_pcol, atracer_pcol)
    cb = plt.colorbar(cs)

    cmin, cmax = make_percentile_cmap(atracer_pcol, 5)
    cs.set_clim(cmin, cmax)

    if not istracer:
        cmin, cmax = cs.get_clim()
        cmax = np.max(np.abs([cmin, cmax]))
        cs.set_clim(-cmax, cmax)

    ax.set_facecolor("gray")
    ax.set_xlim(lmin, lmax)
    ax.set_ylim(zmin, zmax)
    ax.set_xticks(np.linspace(lmin, lmax, 5))

    return cs, cb


def contourplot(secstruct, vartoplot, istracer, itime=None, ax=None, **kwargs):

    """ Contourplot of a variable contained in a gridded section.

    :param pypago.sections.GridSection secstruct: Gridded section
    :param str vartoplot: Variable to plot
    :param bool istracer: True if tracer field (i.e. temperature),
     False if velocity field (or transport field).
    :param int itime: Time index to extract. If None,
     temporal mean is computed
    :param ax ax: Axis. If None, gca()

    :return: cs, cb

    """

    atracer, zvect, lvect = _data_prep(secstruct, vartoplot, 'contour', istracer, itime)

    zmin, zmax = zvect.min(), zvect.max()
    lmin, lmax = lvect.min(), lvect.max()

    if ax is None:
        ax = plt.gca()

    cl = plt.tricontour(lvect, zvect, atracer, **kwargs)

    ax.set_facecolor("gray")
    ax.set_xlim(lmin, lmax)
    ax.set_ylim(zmin, zmax)
    ax.set_xticks(np.linspace(lmin, lmax, 5))

    return cl


def preplot(secstruct, secdata, istracer):

    """
    Function that prepares a section tracer (T or S array)
    variable to do both a contourf and a pcolor plot.

    :param pago_obj secstruct: The |pypago| object that contains the section
       variables (`veci`, `vecj`, etc.)
    :param numpy.array secdata: Array that contains the value to
       plot (vect or vecs, could be a mean, a snapshot, etc.)
    :param str istracer: True if a tracer field (vect, vecs) is plotted
    """

    nz = len(secstruct.depthvect[:, 0])
    nl = len(secstruct.lvect)

    # barrier.n: removed mask from zvect since it causes problem
    # when plotting. we set to 0 so that no effect on the cumsum
    zvect = secstruct.depthvect.copy()
    zvect = -np.cumsum(zvect, axis=0)
    zeros = np.zeros((nl, 1)).T
    zvect = np.concatenate((zeros, zvect), axis=0)
    zvect = 0.5 * (zvect[1:, :] + zvect[:-1, :])

    lvect = secstruct.lengthvect.copy()
    lvect = np.cumsum(np.tile(lvect, (nz, 1)), axis=1)
    zeros = np.zeros((nz, 1))
    lvect = np.concatenate((zeros, lvect), axis=1)
    lvect = 0.5 * (lvect[:, 1:] + lvect[:, :-1]) 

    secdata = np.ma.masked_where(secdata == 0, secdata)
    isok = ~np.ma.getmaskarray(secdata)
    atracer = np.ravel(secdata[isok])
    zvect = np.ravel(zvect[isok])
    lvect = np.ravel(lvect[isok])
    

    return atracer, zvect, lvect