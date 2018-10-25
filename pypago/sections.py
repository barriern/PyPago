""" Module for the manipulation of section endpoints and gridded section """

import numpy as np
import pylab as plt
try:
    from param import ee
except ImportError:
    from pypago.sample_param import ee
import pypago.toolsec as toolsec
from pypago.disp import PypagoErrors
import pypago.misc


class Section(object):

    """ Section endpoint object. """

    def __init__(self, name, lon, lat, dire):

        """

        :param str name: Name of the section
        :param numpy.array lon: Section endpoints longitudes
        :param numpy.array lat: Section endpoints latitudes
        :param numpy.array dire: Section segments directions

        """
        self.name = name
        self.lat = np.array(lat)
        self.lon = np.array(lon)
        self.dire = np.array(dire)

    def __str__(self):

        """ Redefinition of the string function """

        output = 'Section %s\n' % self.name
        output += '    -lon: %s\n' % self.lon
        output += '    -lat: %s\n' % self.lat
        output += '    -dire: %s\n' % self.dire
        return output


#    @property
#    def name(self):
#        """ Section name """
#        return self.__name
#
#    @property
#    def lon(self):
#        """ Section endpoints longitude """
#        return self.__lon
#
#    @property
#    def lat(self):
#        """ Section endpoints latitude """
#        return self.__lat
#
#    @property
#    def dire(self):
#        """ Section segments direction """
#        return self.__dire
#
#    @name.setter
#    def name(self, name):
#        '''
#        name setter
#        :param str name: Section name
#        '''
#        self.__name = np.array(name)
#
#    @lon.setter
#    def lon(self, lon):
#        '''
#        lon setter
#        :param str lon: Section lon
#        '''
#        self.__lon = np.array(lon)
#
#    @lat.setter
#    def lat(self, lat):
#        '''
#        lat setter
#        :param str lat: Section lat
#        '''
#        self.__lat = np.array(lat)
#
#    @dire.setter
#    def dire(self, dire):
#        '''
#        dire setter
#        :param str dire: Section dire
#        '''
#        self.__dire = np.array(dire)


class GridSection(object):

    """ Gridded section object """

    def __init__(self, grid, section):

        """

        :param pypago.grid.Grid grid: Model grid
        :param list section: List of :py:class:`pypago.sections.Section` objects

        """

        self.name = None
        self.modelname = None
        self.i = None
        self.j = None
        self.dire = None
        self.faces = None
        self.orient = None
        self.lengthvect = None
        self.areavect = None
        self.depthvect = None
        self.lvect = None
        self.veci = None
        self.vecj = None

        if grid is None:
            self.jmin = None
            self.jmax = None
            self.imin = None
            self.imax = None
            self.nlon = None

        elif (grid is not None) and (section is not None):
            self.jmin = grid.jmin
            self.jmax = grid.jmax
            self.imin = grid.imin
            self.imax = grid.imax
            self.nlon = grid.nlon
            self.modelname = grid.modelname

            self._initsection(grid, section)
            if self.name is not None:
                # It the initsection worked (name is not None) then we
                # finalise the section
                self._finalisesection(grid)
        else:
            message = 'The grid and section arguments must BOTH be'
            message += 'None or not None. This is not '
            message += 'the case.'
            raise PypagoErrors(message)

    def __str__(self):

        """ Redefinition of the str function """

        output = 'Gridded section, %s model:\n' % self.modelname

        attrnames = pypago.misc.extract_attrlist(self)
        output += pypago.misc.extract_str(attrnames, self)

        return output

    def _initsection(self, grid, section):

        """
        Function which allows to initialise the model sections from the
        sections' endpoints, the model longitude and latitudes.
        The index of the sections which are out of the domain
        (should not be further processed)
        are stored in a list. If the section is in the domain, the i,j
        indexes of its
        endpoints and the directions of the segments are stored in a
        |pypago| object.

        :param list sections: list of |pypago| objects
           that contains the sections' endpoints
        :param numpy.array lont: longitude of the model subdomain
        :param numpy.array latt: latitude of the model subdomain
        :return: a list (`index_badsec`) that contains the indexes of
           the bad sections and a list (`goodsec`) that contains
           the |pypago| objects of the good sections, containing the section
           endpoints locations in the model grid

        """

        # We try to use the locingrid function.
        # If failed, we append the index of the section to the bad section index list
        try:
            [vecj, veci] = toolsec.locingrid(section.lon,
                                             section.lat,
                                             grid.lont,
                                             grid.latt)
            dire = section.dire

            if veci[-1] < veci[0]:
                veci = veci[::-1]
                vecj = vecj[::-1]
                dire = dire[::-1]

            self.i = veci
            self.j = vecj
            self.dire = dire
            self.name = section.name

        except IOError:
            pass

    def _finalisesection(self, grid):

        """
        Function that performs the computation of section staircases.
        The `sections` arguments are the good sections which have been
        returned by the :py:func:`pypago.sec.initsections`.
        It creates a list of |pypago| objects, each containing the sections'
        properties (`areavect`, `veci`, `vecj`, `depthvect`, `lvect`, etc)

        :param list sections: list of the good sections, obtained by the
            :py:func:`pypago.sec.initsections` function

        :param pago_obj grid: Pago
            object that contains all the model attributes (area, depthvect,
            scale factor, etc...)

        """

        [vecj, veci] = toolsec.secingrid(self.i[0], self.j[0],
                                         self.i[1], self.j[1])

        [faces, newveci, newvecj, orientation] = toolsec.facesinsec(veci, vecj,
                                                                    self.dire[0])

        self.veci = newveci
        self.vecj = newvecj
        self.faces = faces
        self.orient = orientation

        nl1 = np.floor(len(newveci)/2.)
        if len(self.i) > 2:
            for l in xrange(2, len(self.i)):
                [vecj, veci] = toolsec.secingrid(self.i[l-1], self.j[l-1],
                                                 self.i[l], self.j[l])

                [faces, newveci, newvecj, orientation] = toolsec.facesinsec(veci, vecj, self.dire[l-1])

                [finalfaces, finalveci, finalvecj, finalorient] = toolsec.consec(self.veci, self.vecj, self.faces,
                                                                                 self.orient, newveci, newvecj,
                                                                                 faces, orientation)

                self.veci = finalveci
                self.vecj = finalvecj
                self.faces = finalfaces
                self.orient = finalorient
                nl2 = np.floor(len(newveci)/2.)

        # in case section is closed
        if (self.veci[0] == self.veci[-1]) & (self.vecj[0] == self.vecj[-1]):

            [finalfaces, finalveci, finalvecj, finalorient] = \
                    toolsec.consec(self.veci[-nl2-1:], self.vecj[-nl2-1:],
                                   self.faces[-nl2-1:], self.orient[-nl2:-1],
                                   self.veci[:nl1], self.vecj[:nl1],
                                   self.faces[:nl1], self.orient[:nl1])

            self.veci = np.append(self.veci[nl1:-nl2-1], finalveci)
            self.vecj = np.append(self.vecj[nl1:-nl2-1], finalvecj)
            self.faces = np.append(self.faces[nl1:-nl2-1], finalfaces)
            self.orient = np.append(self.orient[nl1:-nl2-1], finalorient)

        self.lengthvect = toolsec.lengthinsec(self.veci, self.vecj, self.faces,
                                              grid.dyw, grid.dxn)

        # if 'e3t_ps' in grid.__dict__.keys():
        # self.depthvect = _depthinsec_e3tps(self.veci, self.vecj,
        # self.faces, grid)
        # self.areavect = self.depthvect*np.tile(np.transpose(self.lengthvect[:, np.newaxis]), (grid.areaW.shape[0], 1))
        # else:
        self.areavect = toolsec.areainsec(self.veci, self.vecj, self.faces,
                                          grid.areaw, grid.arean)
        self.depthvect = toolsec.areainsec(self.veci, self.vecj, self.faces,
                                           grid.dzw, grid.dzn)

        nl = len(np.nonzero(self.veci == self.veci)[0])
        [vecj, veci] = toolsec.nodouble(self.veci[:nl], self.vecj[:nl])
        self.lvect = np.atleast_1d([0])
        for l in xrange(1, len(veci)):
            self.lvect = np.append(self.lvect, self.lvect[l-1] +
                                   toolsec.distance(grid.latt[vecj[l], veci[l]], grid.lont[vecj[l], veci[l]],
                                                    grid.latt[vecj[l-1], veci[l-1]], grid.lont[vecj[l-1], veci[l-1]], ee))

    def plotsecfaces(self, axes=None, **kwargs):

        """

        Draws a gridded section as 'staircase',
        with grid points as the x and y coordinates.

        Should be plotted on a map background that shows
        the `mask` variables. Must be used to check whether
        the all the dots associated with one section are on the same side
        of the line

        :param matplotlib.axes.Axes axes: The
           :py:class:`matplotlib axes <matplotlib:matplotlib.axes.Axes>` instance
           where to draw the lines

        :param dict kwargs: Additional arguments to the plot function
        :return: The color of the line
        :rtype: list|str|unicode

        """

        if axes is None:
            axes = plt.gca()

        lw = 1
        cpt = 1

        for indl in xrange(0, len(self.veci)):

            if self.faces[indl] == 'N':
                if cpt:
                    lll = axes.plot([self.veci[indl]-.5, self.veci[indl]+.5],
                                    [self.vecj[indl]+.5, self.vecj[indl]+.5], **kwargs)
                    col = lll[0].get_color()
                    cpt = 0
                else:
                    lll = axes.plot([self.veci[indl]-.5, self.veci[indl]+.5],
                                    [self.vecj[indl]+.5, self.vecj[indl]+.5],
                                    color=col)

                if self.orient[indl] == 1:
                    axes.scatter(self.veci[indl], self.vecj[indl]+.75,
                                 lw, marker='.', color=col)
                else:
                    axes.scatter(self.veci[indl], self.vecj[indl]+.25,
                                 lw, marker='.', color=col)

            else:
                if cpt:
                    lll = axes.plot([self.veci[indl]-.5, self.veci[indl]-.5],
                                    [self.vecj[indl]-.5, self.vecj[indl]+.5], **kwargs)
                    col = lll[0].get_color()
                    cpt = 0
                else:
                    lll = axes.plot([self.veci[indl]-.5, self.veci[indl]-.5],
                                    [self.vecj[indl]-.5, self.vecj[indl]+.5],
                                    color=col)

                if self.orient[indl] == 1:
                    axes.scatter(self.veci[indl]-.25, self.vecj[indl],
                                 lw, marker='.', color=col)
                else:
                    axes.scatter(self.veci[indl]-.75, self.vecj[indl],
                                 lw, marker='.', color=col)

        axes.text(np.mean(self.veci), np.mean(self.vecj), self.name, color=col)

        return col


def extract_grid_sections(grid, sectionslist):

    """
    Extract a list of GridSection objects from a grid
    object (containing coordinates and scale factors) and
    a list of section endpoints.

    :param pypago.grid.Grid grid: Input
     grid containing the coordinates, mask and scale factors.

    :param list sectionslist: List of pypago.sections.Section
     objects, containing section endpoints definitions.

    :return: A tuple containing the section in the model world
     (pypago.sections.GridSection objects) and the indexes of
     the discarded sections (i.e. sections out of the domain)

    """

    # if the input is a section, converion into a list
    if isinstance(sectionslist, Section):
        sectionslist = [sectionslist]

    if not isinstance(sectionslist, list):
        message = 'The sectionslist argument must either be a single or a list of '
        message += 'pypago.sections.Section objects. This program will '
        message += 'be stopped.'
        raise PypagoErrors(message)

    output = []
    badsection = []
    for secind in xrange(0, len(sectionslist)):
        gridsec = GridSection(grid, sectionslist[secind])
        if gridsec.name is not None:
            output.append(gridsec)
        else:
            badsection.append(secind)

    return output, badsection


def correct_gridsec(gridsec, secname, offset, position):

    """
    Function that allows to correct section staircases.
    Especially useful in order to correct sections' junctions. It extracts
    the section variables (`veci`, `vecj`, etc) along
    the length coordinates, from `offset` to `end` or from
    `start` to `offset`, depending on the value of `position` argument.

    :param list sec: section list that contains the |pypago|
       objects, containing the veci, vecj, vect, ... obtained
       from the :py:func:`pypago.pypago_sec.finalisesections`

    :param type secname: section points to correct

    :param int offset: number of points to remove

    :param str position: {'start','end'}
       whether we remove the first (position='start')
       or last (position='end') points.

    .. warning::

        The input list is modified

    """

    nw = pypago.misc.findsecnum(gridsec, secname)

    if position == 'end':

        gridsec[nw].veci = gridsec[nw].veci[:-offset]
        gridsec[nw].vecj = gridsec[nw].vecj[:-offset]
        gridsec[nw].faces = gridsec[nw].faces[:-offset]
        gridsec[nw].orient = gridsec[nw].orient[:-offset]
        gridsec[nw].areavect = gridsec[nw].areavect[:, :-offset]
        gridsec[nw].depthvect = gridsec[nw].depthvect[:, :-offset]
        gridsec[nw].lvect = gridsec[nw].lvect[:-offset]
        gridsec[nw].lengthvect = gridsec[nw].lengthvect[:-offset]

    else:
        gridsec[nw].veci = gridsec[nw].veci[offset:]
        gridsec[nw].vecj = gridsec[nw].vecj[offset:]
        gridsec[nw].faces = gridsec[nw].faces[offset:]
        gridsec[nw].orient = gridsec[nw].orient[offset:]
        gridsec[nw].areavect = gridsec[nw].areavect[:, offset:]
        gridsec[nw].depthvect = gridsec[nw].depthvect[:, offset:]
        gridsec[nw].lvect = gridsec[nw].lvect[offset:]
        gridsec[nw].lengthvect = gridsec[nw].lengthvect[offset:]

    return gridsec
