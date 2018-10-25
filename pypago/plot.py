
"""
Module that contains various functions dedicated to plotting using |pypago|
"""

from __future__ import print_function
import numpy as np
import pylab as plt
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
    cmap = ListedColormap(['black', 'white', 'gray'])

    plt.imshow(grid.mask, interpolation='none', cmap=cmap)
    if mask is not None:
        plt.imshow(mask, interpolation='none', cmap=cmap2)
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

    atracer, zvect, lvect, atracer_pcol, zvect_pcol, lvect_pcol = preplot(secstruct, secdata, istracer)

    atracer_pcol = np.ma.masked_where(extract_mask(atracer_pcol), atracer_pcol)
    atracer = np.ma.masked_where(extract_mask(atracer), atracer)

    lvect_pcol = np.ma.masked_where(extract_mask(lvect_pcol), lvect_pcol)
    lvect = np.ma.masked_where(extract_mask(lvect), lvect)

    zvect_pcol = np.ma.masked_where(zvect_pcol == 0, zvect_pcol)
    lvect_pcol = np.ma.masked_where(np.ma.getmaskarray(zvect_pcol), lvect_pcol)

    if output_form == 'pcolor':
        output = (atracer_pcol, zvect_pcol, lvect_pcol)
    else:
        output = (atracer, zvect, lvect)

    return output


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

    cs = plt.pcolormesh(lvect_pcol, zvect_pcol, atracer_pcol)
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

    cl = plt.contour(lvect, zvect, atracer, **kwargs)

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

    atracer = np.zeros((nz, nl))
    zvect = np.zeros((nz, nl))

    for indz in xrange(0, nz):
        atracer[indz, :] = _nodouble_tsr(secstruct.veci,
                                         secstruct.vecj, secdata[indz, :])
        if istracer:
            zvect[indz, :] = _nodouble_tsr(secstruct.veci, secstruct.vecj,
                                           secstruct.depthvect[indz, :])
        else:
            zvect[indz, :] = _nodouble_v(secstruct.veci, secstruct.vecj,
                                         secstruct.depthvect[indz, :])

    # barrier.n: removed mask from zvect since it causes problem
    # when plotting. we set to 0 so that no effect on the cumsum
    zvect = np.array(zvect)
    zvect = np.ma.masked_where(np.ma.getmaskarray(zvect), zvect)
    zvect = np.ma.masked_where(np.isnan(zvect), zvect)
    zvect = -np.cumsum(zvect, axis=0)

    lvect = np.tile(secstruct.lvect, (nz, 1))

    lvect_pcol = _make_lvect_pcol(lvect)
    zvect_pcol = _make_zvect_pcol(zvect)
    atracer_pcol = _make_atracer_pcol(atracer)

    return atracer, zvect, lvect, \
        atracer_pcol, zvect_pcol, lvect_pcol


def _nodouble_v(veci, vecj, vecv):

    """
    Returns the velocity field after removing
    the duplicated values of (`veci`, `vecj`)

    .. warning::
        Relies on the assumption that all sections go
        from West to East

    :param numpy.array veci: i indexes of the section
    :param numpy.array vecj: j indexes of the section
    :param numpy.array vecv: velocity field of the section

    :return: the new velocity

    :rtype: numpy.array

    """

    veci = np.squeeze(veci)
    vecj = np.squeeze(vecj)
    vecv = np.squeeze(vecv)

    newvecv = np.ma.masked_array([])
    newvecv = np.ma.append(newvecv, vecv[0])

    for indl in xrange(1, veci.shape[0]):
        if not (veci[indl] == veci[indl-1]) & (vecj[indl] == vecj[indl-1]):
            newvecv = np.ma.append(newvecv, vecv[indl])
        else:
            temp = np.arctan(vecv[indl]/(-vecv[indl-1]))
            if -vecv[indl-1] < 0:
                if vecv[indl] >= 0:
                    temp += np.pi
                else:
                    temp -= np.pi

            if (temp > np.pi/4.) | (temp < -3*np.pi/4.):
                signe = 1
            else:
                signe = -1

            # barrier.n: correction according to Matlab's version
            # if vecv[indl-1]*vecv[indl-1] + vecv[indl]*vecv[indl] >= 0:
            #     newvecv[-1] = signe*np.sqrt(vecv[indl-1]*vecv[indl-1] +
            #                                 vecv[indl]*vecv[indl])
            # else:
            #     newvecv[-1] = np.nan
            newvecv[-1] = signe * np.sqrt(vecv[indl - 1]**2 + vecv[indl]**2)

    return newvecv


def _nodouble_tsr(veci, vecj, vectsr):

    """
    Returns the tracer field after removing
    the duplicated values of (`veci`, `vecj`)

    :param numpy.array veci: i indexes of the section
    :param numpy.array vecj: j indexes of the section
    :param numpy.array vectsr: tracer field of the section
    :return: the tracer field without the duplications
    :rtype: list
    """

    veci = np.squeeze(veci)
    vecj = np.squeeze(vecj)
    vectsr = np.squeeze(vectsr)

    newvectsr = []
    newvectsr.append(vectsr[0])

    for indl in xrange(1, veci.shape[0]):

        if not (veci[indl] == veci[indl-1]) & (vecj[indl] == vecj[indl-1]):
            newvectsr.append(vectsr[indl])

    newvectsr = np.array(newvectsr)
    return newvectsr


def _make_lvect_pcol(lvect):

    """
    Returns the length vector for a pcolor plot

    :param numpy.array lvect: length vector
    :return: array that contains the `lvect` array formatted for a pcolor plot
    :rtype: numpy.array
    """

    lvect = np.ma.masked_where(extract_mask(lvect), lvect)

    lvect_pcol = 1.5 * lvect[:, 0:1] - 0.5 * lvect[:, 1:2]
    lvect_pcol = np.concatenate((lvect_pcol, .5 * (lvect[:, :-1] + lvect[:, 1:])), axis=1)
    lvect_pcol = np.concatenate((lvect_pcol, 1.5 * lvect[:, -1:] - 0.5 * lvect[:, -2:-1]), axis=1)
    lvect_pcol = np.concatenate((lvect_pcol[0:1, :], lvect_pcol), axis=0)

    return lvect_pcol


def _make_zvect_pcol(zvect):

    """
    Returns the depth vector for a pcolor plot.

    :param numpy.array zvect: vector that contains the
      `zvect` array without the duplications (obtained with the
      :py:func:`pypago.plot._nodouble_tsr` module)
    :return: array that contains the `zvect` array formatted for a pcolor plot
    :rtype: numpy.array
    """

    nz, nl = zvect.shape
    zvect = np.ma.masked_where(extract_mask(zvect), zvect)

    # initialisation of a temporary zvect_pcol array
    zvect_pcol_temp = []
    for indl in xrange(0, nl):
        zl = len(np.nonzero(np.ma.getmaskarray(zvect[:, indl]) == 0)[0])
        if zl > 1:
            temp = np.concatenate(([0],
                                   .5*(zvect[:zl-1, indl] + zvect[1:zl, indl]),
                                   1.5*zvect[zl-1:zl, indl] - 0.5*zvect[zl-2:zl-1, indl],
                                   zvect[zl:nz+1, indl]))
        else:
            temp = np.concatenate(([0], zvect[0:1, indl],
                                   zvect[1:nz, indl]), axis=0)
        zvect_pcol_temp.append(temp)

    zvect_pcol_temp = np.array(zvect_pcol_temp)
    zvect_pcol_temp = np.ma.masked_where(zvect_pcol_temp != zvect_pcol_temp,
                                         zvect_pcol_temp).T

    # creation of the final zvect_pcol array
    zvect_pcol = np.zeros((zvect_pcol_temp.shape[0],
                           zvect_pcol_temp.shape[1]+1))
    zvect_pcol[:, 0] = zvect_pcol_temp[:, 0]

    for indl in xrange(1, len(zvect_pcol_temp[0, :])):

        ind1 = np.nonzero(np.ma.getmaskarray(zvect_pcol_temp[:, indl]) == 0)[0]
        ind2 = np.nonzero(np.ma.getmaskarray(zvect_pcol_temp[:, indl-1]) == 0)[0]
        zmin = np.min([len(ind1), len(ind2)])
        zmax = np.max([len(ind1), len(ind2)])

        zvect_pcol[:zmin, indl] = 0.5*(zvect_pcol_temp[:zmin, indl-1] + zvect_pcol_temp[:zmin, indl])

        for indz in xrange(zmin, zmax):
            if np.ma.getmaskarray(zvect_pcol_temp[indz, indl]):
                zvect_pcol[indz, indl] = zvect_pcol_temp[indz, indl-1]
            else:
                zvect_pcol[indz, indl] = zvect_pcol_temp[indz, indl]
        if zmax <= nz:
            zvect_pcol[zmax:, indl] = np.nan

    zvect_pcol[:, -1] = zvect_pcol[:, -2]
    zvect_pcol = np.ma.array(zvect_pcol, mask=np.isnan(zvect_pcol))

    return zvect_pcol


def _make_atracer_pcol(atracer):

    """
    Returns the field vector
    for a pcolor plot

    :param numpy.array atracer: vector that contains the
     field array without the duplications (obtained with the
     :py:func:`pypago.plot._nodouble_tsr` or
     :py:func:`pypago.plot._nodouble_v`)
    :return: array that contains the field array formatted for a pcolor plot
    :rtype: numpy.array

    """

    nz, nl = atracer.shape

    atracer_pcol = np.concatenate((atracer, np.nan*np.ones((nz, 1))), axis=1)
    atracer_pcol = np.concatenate((atracer_pcol, np.nan*np.ones((1, nl+1))), axis=0)

    return atracer_pcol
