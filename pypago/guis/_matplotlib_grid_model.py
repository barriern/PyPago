# -*- coding: utf-8 -*-

"""
Module that handles the Matplotlib figure associated with
the :py:class:`pypago.pypago_gui.gui_grid_model.GridModel`
class
"""

import tkMessageBox
import pypago.pyio
import pypago.grid
import pypago.plot
import pypago.sections
import pypago.coords
import pypago.misc
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvas
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from matplotlib import rcParams
rcParams['lines.linewidth'] = 1


class MatplotlibGridModel(object):

    """
    Class that defines the GUIs associated with the program
    for the grid model definition. Inherits from the
    :py:class:`Tkinter.Tk` class.

    .. versionadded:: 20150825
       barrier.n

    """

    def __init__(self, gui):

        # GUI associated with the module
        self.gui = gui

        # domain coordinates
        self.imin = None
        self.imax = None
        self.jmin = None
        self.jmax = None

        # dictionary that contains the coordinates of
        # the global domain and the domain limits
        self.dic_coord = None
        self.latlim = None
        self.lonlim = None

        # dictionary that contains the coordinates
        # and scale factors of the subdomain
        self.gridmodel = None

        # initialisation of the matplotlib events
        self.dompic_event = None
        self.editsec_event = None
        self.dom_point_move_event = None
        self.dom_point_release_event = None
        self.sec_point_move_event = None
        self.sec_point_release_event = None

        # initialisation of some stuff related to the matplotlib
        # events
        self.xdomain = None
        self.ydomain = None
        self.xdomain2 = None
        self.xdomain1 = None
        self.xclick = None
        self.yclick = None
        self.deltax = None
        self.deltay = None

        # matplotlib object containing the scatter plot
        self.scatter = None

        # attributes related to the section edition
        # list of sections read on the files (endpoints)
        self.sections = None
        # list of sections processed ('bad' sections are discarded)
        self.model_sections = None
        # section name
        self.secname = None
        # indice of the section
        self.indsec = None
        # indice of the segment
        self.segment = None
        # segment direction
        self.segdir = None

        # indice of the point selected on the scatter plots
        # used in domain/section edition
        self.indpoint = None

        # creation of the Matplotlib frame
        self.figure = matplotlib.figure.Figure(figsize=(2, 2))
        self.plotax = self.figure.add_axes([0.1, 0.15, 0.8, 0.8])
        self.canvas = FigureCanvas(self.figure, master=gui)

    def init_domain(self):

        """
        Function that is called when a |netcdf| file is open.

        It loads the model grid file, draw the masks
        and initialise the MODEL_grid |pypago| object
        """

        # we first disconnect the dompic and editsec events (if possible)
        # prevents from conflicting events
        self.disconnect_events()

        # we load the global coordinates of the file
        self.dic_coord = pypago.coords.create_coord(self.gui.modelproc, self.gui.gridfilename)

        # We define the maximum lon/lat limits of the domain
        # Longitude: 0:END
        # Latitude: 1:END-1
        self.latlim = [1, self.dic_coord.mask.shape[0] - 2]
        self.lonlim = [0, self.dic_coord.mask.shape[1] - 1]

        # initialisation of the domain attributes
        self.imin = 0
        self.imax = self.dic_coord.mask.shape[1] - 1
        self.jmin = 1
        self.jmax = self.dic_coord.mask.shape[0] - 2

        # we activate the widgets for the domain edition
        self.gui.entry_imin.configure(state='normal')
        self.gui.entry_imax.configure(state='normal')
        self.gui.entry_jmin.configure(state='normal')
        self.gui.entry_jmax.configure(state='normal')

        # We display the default domain values
        self.gui.entry_imin_lab.set(str(self.imin))
        self.gui.entry_imax_lab.set(str(self.imax))
        self.gui.entry_jmin_lab.set(str(self.jmin))
        self.gui.entry_jmax_lab.set(str(self.jmax))

        # we draw the global mask
        self.plot_map_global()

        # we draw the domain
        self.plot_dom()

        # we activate the mpl domain event
        self.dompic_event = self.figure.canvas.mpl_connect('pick_event', self.on_dompic)

    def plot_map_global(self):

        """
        Function that plots the 'global' mask'
        """

        self.plotax.cla()
        print self.dic_coord.mask.min(), self.dic_coord.mask.max()
        self.plotax.contour(self.dic_coord.mask, levels=[1 - np.spacing(1), 1], colors='k')
        self.plotax.set_xlim(self.lonlim)
        self.plotax.set_ylim(self.latlim)
        self.plotax.set_xlabel('grid points')
        self.plotax.set_ylabel('grid points')

    def plot_dom(self):

        """
        Function that handles the plotting of the domain
        """

        self.plotax.lines = []  # We remove the lines

        try:
            self.scatter.remove()
        except:
            pass

        # we define the domain plot as a closed box
        self.ydomain = [self.jmin, self.jmax, self.jmax, self.jmin, self.jmin]
        self.xdomain = [self.imin, self.imin, self.imax, self.imax, self.imin]

        if self.imin < self.imax:  # If western lon < eastern lon (classical domain)

            self.plotax.plot(self.xdomain, self.ydomain, linewidth=3, color='black')
            # We draw the scatter plot
            self.scatter = self.plotax.scatter(self.xdomain[:4], self.ydomain[:4],
                                               100, color='DarkOrange', zorder=1000, picker=10)

        else:  # If western lon > eastern lon (discontinuous domain)

            # we define two boxes
            self.xdomain1 = [self.imin, self.imin, self.lonlim[-1],
                             self.lonlim[-1], self.imin]  # left box
            self.xdomain2 = [self.lonlim[0], self.lonlim[0],
                             self.imax, self.imax, self.lonlim[0]]  # right box

            # array that contains the longitude corners of the box
            xxx = [self.imin, self.imin, self.imax, self.imax, self.imin]

            # We draw separately the 2 boxes
            # right box
            self.plotax.plot(self.xdomain1, self.ydomain, linewidth=3, color='black')
            # left box
            self.plotax.plot(self.xdomain2, self.ydomain, linewidth=3, color='black')

            # We draw the four corners as scatter plot
            self.scatter = self.plotax.scatter(xxx[:4], self.ydomain[:4],
                                               100, color='DarkOrange', zorder=1000, picker=5)

        # we draw the figure
        self.canvas.draw()

    def on_dompic(self, event):    # pylint: disable=unused-argument

        """
        Matplotlib event that handles the displacement of the points (only way to edit the domain)
        """

        # self.indpoint = value of the point that has been selected
        # (0=bottom left, 1=top left, 2=top right, 3=bottom right)
        self.indpoint = event.ind[0]

        # x0,y0 contains the point of the domain on which which have clicked
        self.xclick = self.xdomain[self.indpoint]
        self.yclick = self.ydomain[self.indpoint]

        # we activate the point move and point release event
        self.dom_point_move_event = self.figure.canvas.mpl_connect('motion_notify_event',
                                                                   self.on_dompointmove)
        self.dom_point_release_event = self.figure.canvas.mpl_connect('button_release_event',
                                                                      self.on_dompointrelease)

    def on_dompointmove(self, event):    # pylint: disable=unused-argument

        """
        Matplotlib event that handles the scatter points displacements for the domain
        """

        if event.inaxes:  # if the click is on the axes

            # deltax,deltay is the mouse displacement
            deltax = event.xdata - self.xclick
            deltay = event.ydata - self.yclick

            if self.imin < self.imax:  # if we have a classical domain (i.e one single box)

                # We create two temporary arrays that contain the domain coordinates
                xxx = np.empty(len(self.xdomain[:]))
                xxx[:] = self.xdomain
                yyy = np.empty(len(self.ydomain[:]))
                yyy[:] = self.ydomain

                # if bottom left point, the top left and bottom right points are modified
                if self.indpoint == 0:
                    xxx[0] = xxx[0] + deltax
                    yyy[0] = yyy[0] + deltay
                    xxx[1] = xxx[1] + deltax
                    yyy[3] = yyy[3] + deltay

                # if top left point, the bottom left and top right points are modified
                if self.indpoint == 1:
                    xxx[1] = xxx[1] + deltax
                    yyy[1] = yyy[1] + deltay
                    xxx[0] = xxx[0] + deltax
                    yyy[2] = yyy[2] + deltay

                # if top right point, the bottom right and top left points are modified
                if self.indpoint == 2:
                    xxx[2] = xxx[2] + deltax
                    yyy[2] = yyy[2] + deltay
                    xxx[3] = xxx[3] + deltax
                    yyy[1] = yyy[1] + deltay

                # if bottom right point, the top right and bottom left points are modified
                if self.indpoint == 3:
                    xxx[3] = xxx[3] + deltax
                    yyy[3] = yyy[3] + deltay
                    xxx[2] = xxx[2] + deltax
                    yyy[0] = yyy[0] + deltay

                # we close the domain
                xxx[-1] = xxx[0]
                yyy[-1] = yyy[0]

                # we modify the line coordinates
                self.plotax.lines[0].set_xdata(xxx)
                self.plotax.lines[0].set_ydata(yyy)

                # we remove and redraw the scatter plot
                self.scatter.remove()
                self.scatter = self.plotax.scatter(xxx[:4], yyy[:4], 100, marker='o',
                                                   picker=5, color='DarkOrange', zorder=10000)

            else:

                # same thing as before but with two boxes, hence two temporary
                # x coordinates (xxx1 and xxx2).
                # we define only one yyy array

                xxx1 = np.empty(len(self.xdomain1[:]))
                xxx1[:] = self.xdomain1

                xxx2 = np.empty(len(self.xdomain2[:]))
                xxx2[:] = self.xdomain2

                yyy = np.empty(len(self.ydomain[:]))
                yyy[:] = self.ydomain

                # bottom left point
                if self.indpoint == 0:
                    xxx1[0] = xxx1[0] + deltax
                    yyy[0] = yyy[0] + deltay
                    xxx1[1] = xxx1[1] + deltax
                    yyy[3] = yyy[3] + deltay

                # top left point
                if self.indpoint == 1:
                    xxx1[1] = xxx1[1] + deltax
                    yyy[1] = yyy[1] + deltay
                    xxx1[0] = xxx1[0] + deltax
                    yyy[2] = yyy[2] + deltay

                # top right point
                if self.indpoint == 2:
                    xxx2[2] = xxx2[2] + deltax
                    yyy[2] = yyy[2] + deltay
                    xxx2[3] = xxx2[3] + deltax
                    yyy[1] = yyy[1] + deltay

                # bottom right point
                if self.indpoint == 3:
                    xxx2[3] = xxx2[3] + deltax
                    yyy[3] = yyy[3] + deltay
                    xxx2[2] = xxx2[2] + deltax
                    yyy[0] = yyy[0] + deltay

                # we close the domain
                xxx1[-1] = xxx1[0]
                xxx2[-1] = xxx2[0]
                yyy[-1] = yyy[0]

                # we change the coordinates of the two boxes
                self.plotax.lines[0].set_xdata(xxx1)
                self.plotax.lines[1].set_xdata(xxx2)
                self.plotax.lines[0].set_ydata(yyy)
                self.plotax.lines[1].set_ydata(yyy)

                # xxx=valeur of the points on the scatter plot
                xxx = [xxx1[0], xxx1[1], xxx2[2], xxx2[3]]

                # we remove and redo the scatter plot
                self.scatter.remove()
                self.scatter = self.plotax.scatter(xxx, yyy[:4], 100, marker='o',
                                                   picker=5, color='DarkOrange', zorder=10000)

            self.canvas.draw()  # on update le plot

    def on_dompointrelease(self, event):    # pylint: disable=unused-argument

        """
        Function that is called when we release a point of the model
        """

        if self.imin < self.imax:  # if classical domain

            # we change the imin, imax, jmin and jmax attributes of the MODEL_grid
            self.imin = np.round(np.min(np.array(self.plotax.lines[0].get_xdata()))).astype(np.int)
            self.imax = np.round(np.max(np.array(self.plotax.lines[0].get_xdata()))).astype(np.int)
            self.jmin = np.round(np.min(np.array(self.plotax.lines[0].get_ydata()))).astype(np.int)
            self.jmax = np.round(np.max(np.array(self.plotax.lines[0].get_ydata()))).astype(np.int)

        else:  # if we have a discontinuous domain

            self.imin = np.round(np.min(np.array(self.plotax.lines[0].get_xdata()))).astype(np.int)
            self.imax = np.round(np.max(np.array(self.plotax.lines[1].get_xdata()))).astype(np.int)
            self.jmin = np.round(np.min(np.array(self.plotax.lines[0].get_ydata()))).astype(np.int)
            self.jmax = np.round(np.max(np.array(self.plotax.lines[0].get_ydata()))).astype(np.int)

        # we change the display of the domain txt ctr boxes
        self.gui.entry_imin_lab.set(str(self.imin))
        self.gui.entry_imax_lab.set(str(self.imax))
        self.gui.entry_jmin_lab.set(str(self.jmin))
        self.gui.entry_jmax_lab.set(str(self.jmax))
        self.plot_dom()  # we redraw the domain (because we have rounded the imin,imax, and so on)

        # we disconnect the point move and point release event
        self.figure.canvas.mpl_disconnect(self.dom_point_move_event)
        self.figure.canvas.mpl_disconnect(self.dom_point_release_event)

    def plot_sections(self):

        """
        Function that draws the section of the section edition
        """

        self.disconnect_events()

        # nsec = number of section that are on the domain
        nsec = len(self.model_sections)

        # we draw the subdomain mask
        self.plot_map_subdomain()

        # we plot the sections that are on the domain
        for sec in range(0, nsec):
            self.plotax.plot(self.model_sections[sec].i, self.model_sections[sec].j,
                             picker=1, label=self.model_sections[sec].name)

        # we connect the section edition mpl event
        self.editsec_event = self.figure.canvas.mpl_connect('pick_event', self.on_editsec)
        self.canvas.draw()

    def on_editsec(self, event):    # pylint: disable=unused-argument

        """
        This is the mpl event for the section edition
        """

        # we activate the widget that prompts the i,j coordinates and the delete button
        self.gui.label_i.configure(state='normal')
        self.gui.label_j.configure(state='normal')
        self.gui.button_del.configure(state='normal')

        # we recover the position of the click
        self.xclick = event.mouseevent.xdata
        self.yclick = event.mouseevent.ydata

        # if click on a point
        if isinstance(event.artist, matplotlib.collections.PathCollection):

            # we recover the section of the point which has been clicked on
            self.indpoint = event.ind

            # we activate the point move and point release event
            self.sec_point_move_event = self.figure.canvas.mpl_connect('motion_notify_event',
                                                                       self.on_secpointmove)
            self.sec_point_release_event = self.figure.canvas.mpl_connect('button_release_event',
                                                                          self.on_secpointrelease)

        # if click on a text
        elif isinstance(event.artist, matplotlib.text.Text):

            # if a segment was selected, we draw it in black
            if self.segdir is not None:
                self.plotax.texts[self.segment].set_color('k')

            # we recover the number of the segment and the direction of the segment
            self.segment = int(event.artist.get_text()) - 1
            self.segdir = self.model_sections[self.indsec].dire[self.segment]

            # we draw the selected segment in red
            self.plotax.texts[self.segment].set_color('r')

            # we enable the widgets for the section edition
            self.gui.combobox_segment.configure(state='normal')

            # we change the value of the segdir widget
            self.gui.combobox_segment.set(self.segdir)

            # we change the display of the segment
            self.gui.label_segment_lab.set('Segment #' + str(self.segment + 1) + ':')

            # we update the plot
            self.canvas.draw()

        # if click on a line
        elif isinstance(event.artist, matplotlib.lines.Line2D):

            # If we click on a line, we change its display (scatter plot appears)
            # if a section was selected, we remove the text and the scatter plot
            if self.secname is not None:
                self.plotax.texts = []  # remove text
                self.scatter.remove()  # remove scatter plot

            # we recover the name of the section and find its number
            self.secname = event.artist.get_label()
            self.indsec = pypago.misc.findsecnum(self.model_sections, self.secname)

            # we reinitialize the selected segment
            self.segdir = None
            self.gui.combobox_segment.configure(state='disabled')

            # we activate the widget of the section name and update it
            self.gui.label_secname_lab.set('Section name: ' +
                                           str(self.model_sections[self.indsec].name))

            secint = self.model_sections[self.indsec]

            # we draw the scatter plot
            self.scatter = self.plotax.scatter(self.plotax.lines[self.indsec].get_xdata(),
                                               self.plotax.lines[self.indsec].get_ydata(), 50,
                                               marker='o', picker=5, color='DarkOrange',
                                               zorder=1000)

            # initialisation of the string section coordinates
            xstr = []
            ystr = []
            for xpoint, ypoint in zip(secint.i, secint.j):
                xstr.append(str(np.round(xpoint).astype(np.int)))
                ystr.append(str(np.round(ypoint).astype(np.int)))
            x_string = ', '.join(xstr)
            y_string = ', '.join(ystr)

            # we change the prompt of the coordinates
            self.gui.label_i_lab.set('i coordinates: ' + x_string)
            self.gui.label_j_lab.set('j coordinates: ' + y_string)

            # We add the text for the segment directions
            for indice_seg in range(0, len(secint.i) - 1):
                x_text = 0.5 * (secint.i[indice_seg] + secint.i[indice_seg + 1])
                y_text = 0.5 * (secint.j[indice_seg] + secint.j[indice_seg + 1])
                seg_text = str(indice_seg + 1)
                self.plotax.text(x_text, y_text, seg_text, bbox=dict(boxstyle="round", fc="1"),
                                 va='center', ha='center', fontsize=10, picker=5)

            # we force the drawing
            self.canvas.draw()

    def on_secpointrelease(self, event):    # pylint: disable=unused-argument

        """
        Function that is called when the button is released
        """

        # we first disconnect the point move and release mpl events
        self.figure.canvas.mpl_disconnect(self.sec_point_move_event)
        self.figure.canvas.mpl_disconnect(self.sec_point_release_event)

        # we recover the i,j coordinates of the section by rounding the line coordinates
        i = np.round(self.plotax.lines[self.indsec].get_xdata()).astype(np.int)
        j = np.round(self.plotax.lines[self.indsec].get_ydata()).astype(np.int)

        self.model_sections[self.indsec].i = i
        self.model_sections[self.indsec].j = j

        # we recover the lon,lat coordinates of the section
        self.sections[self.indsec].lon = self.gridmodel.lont[j, i]
        self.sections[self.indsec].lat = self.gridmodel.latt[j, i]

        # we update the i,j coordinates on the widgets
        x_string = []
        y_string = []
        for x_point, y_point in zip(i, j):
            x_string.append(str(np.round(x_point).astype(np.int)))
            y_string.append(str(np.round(y_point).astype(np.int)))
        x_string = ', '.join(x_string)
        y_string = ', '.join(y_string)
        self.gui.label_i_lab.set('i coordinates: ' + x_string)
        self.gui.label_j_lab.set('j coordinates: ' + y_string)

    def on_secpointmove(self, event):    # pylint: disable=unused-argument

        """
        Function that handles the maintained click on a point
        """

        # We recover the section that was selected
        secint = self.model_sections[self.indsec]

        # this is the mouse displacement
        self.deltax = event.xdata - self.xclick
        self.deltay = event.ydata - self.yclick

        # this is the coordinates of the section
        xint, yint = (secint.i, secint.j)

        # we create empty coords arrays that have the same coordinates as the section
        xxx = np.empty(secint.i.shape)
        xxx[:] = xint
        xxx[self.indpoint] = xxx[self.indpoint] + self.deltax

        yyy = np.empty(secint.j.shape)  # same as above for y
        yyy[:] = yint
        yyy[self.indpoint] = yyy[self.indpoint] + self.deltay

        # We update the line coordinates
        self.plotax.lines[self.indsec].set_xdata(xxx)
        self.plotax.lines[self.indsec].set_ydata(yyy)

        # we remove and redraw the scatter plot and the text
        self.scatter.remove()
        self.scatter = self.plotax.scatter(xxx, yyy, 50, marker='o', picker=5,
                                           color='DarkOrange', zorder=10000)

        for indice_seg in range(0, len(secint.i) - 1):
            x_text = 0.5 * (xxx[indice_seg] + xxx[indice_seg + 1])
            y_text = 0.5 * (yyy[indice_seg] + yyy[indice_seg + 1])
            self.plotax.texts[indice_seg].set_position([x_text, y_text])

        # we update the plot
        self.canvas.draw()

    def finalize_domain(self):

        """
        Function that computes the grid properties (scale factors,
        etc) on the subdomain
        """
                                    
        self.gridmodel = pypago.grid.create_grid(self.dic_coord,
                                                 self.jmin,
                                                 self.jmax,
                                                 self.imin,
                                                 self.imax)

        #if (gridmodel.dzt.ndim) == 1:
        #ndepth = len(gridmodel.dzt)
        #[nlat, nlon] = gridmodel.dxt.shape
        #gridmodel.dzt = np.tile(gridmodel.dzt, ([1, nlat, nlon]))
        #gridmodel.dzw = np.tile(gridmodel.dzw, ([1, nlat, nlon]))
        #gridmodel.dzn = np.tile(gridmodel.dzn, ([1, nlat, nlon]))
        #else:
        #[ndepth, nlat, nlon] = gridmodel.dzt.shape
        #
        #gridmodel.areaW = gridmodel.dzw * np.tile(gridmodel.dytw, [ndepth, 1, 1])
        #gridmodel.areaN = gridmodel.dzn * np.tile(gridmodel.dxtn, [ndepth, 1, 1])
        #
        #volume = gridmodel.dzt * np.tile(gridmodel.dxt * gridmodel.dyt, [ndepth, 1, 1])
        #volume[volume != volume] = 0
        #surface = gridmodel.dxt * gridmodel.dyt
        #gridmodel.volume = volume
        #gridmodel.surface = surface
        #
        #self.gridmodel = gridmodel

    def generate_sections(self):

        """
        Function that computes the section staircases associated with
        the model grid
        """

        # These sections will not be considered in the output MODEL_sections
        self.model_sections, index_badsec = pypago.sections.extract_grid_sections(self.gridmodel,
                                                                      self.sections)  
        #goodsec, index_badsec = pypago.sec.initsections(self.sections,
        #                                                       self.gridmodel.lont,
        #                                                       self.gridmodel.latt)
        #self.model_sections = goodsec

        #pypago.sections.finalisesections(self.model_sections, self.gridmodel)

        nsec = len(self.model_sections)
        vecname = []
        for indice_section in range(0, len(self.model_sections)):
            vecname.append(self.model_sections[indice_section].name)
        vecname = np.array(vecname)

        # we clean the axes, redraw the mask
        self.plot_map_subdomain()

        # we loop over the good sections and we draw them as steps
        for sec in range(0, nsec):
            col = self.model_sections[sec].plotsecfaces(self.plotax)

        self.canvas.draw()

        # if bad section exists, we prompt a message box with the list of discarded section
        if len(index_badsec):
            badsectionstring = []
            for ibad in index_badsec:
                badsectionstring.append('-' + str(self.sections[ibad].name + '\n'))
                self.sections.pop(ibad)
            stout = ''.join(badsectionstring)
            tkMessageBox.showinfo('Oups!', 'Sections:\n' + stout + ' are out of the domain ' + \
                                  'and are not considered!\n\n' + \
                                  'Consider enlarging your domain and ' + \
                                  'reopening the section file')

    def plot_map_subdomain(self):

        """
        Function that plots the subdomain masks
        """

        self.plotax.cla()
        self.plotax.contour(self.gridmodel.mask, levels=[1 - np.spacing(1), 1], colors='k')
        self.plotax.set_xlim([0, self.gridmodel.mask.shape[1]-1])
        self.plotax.set_ylim([0, self.gridmodel.mask.shape[0]-1])
        self.plotax.set_title('please check carefully that the direction of transport ' + \
                              'across each section is correct using the dots\nif not, edit the sections!')
        self.plotax.set_xlabel('grid points')
        self.plotax.set_ylabel('grid points')

    def disconnect_events(self):

        """
        Function that disconnects the Matplotlib events
        (prevents from conflicts)
        """

        try:
            self.figure.canvas.mpl_disconnect(self.editsec_event)
        except:
            pass

        try:
            self.figure.canvas.mpl_disconnect(self.dompic_event)
        except:
            pass
