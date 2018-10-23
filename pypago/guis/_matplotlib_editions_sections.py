# -*- coding: utf-8 -*-

"""
Module that handles the Matplotlib figure associated with
the :py:class:`pypago.guis.gui_edition_sections.EditionSections`
class
"""

import pypago.pyio
import pypago.misc
import pypago.sections
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvas
import matplotlib
matplotlib.use('TkAgg')
import pylab as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import rcParams
import tkMessageBox
rcParams['lines.linewidth'] = 2
rcParams['text.usetex'] = False


class MatplotlibEditionsSections(object):

    """
    Class that defines the GUIs associated with the program
    for the edition of section endpoints. Inherits from the
    :py:class:`Tkinter.Tk` class.

    .. versionadded:: 2015-07-31
       barrier.n

    """

    def __init__(self, gui):

        """

        Initialisation of the program window.

        :param parent:

        """

        # GUI with which the program interacts
        self.gui = gui

        # List of sections endpoints
        self.sections = []

        # coordinates of the map (ax1) and of the colorbar (ax2)
        self.plotaxcoords = [0.1, 0.15, 0.8, 0.8]
        self.cbaraxcoords = [0, 0.07, 1, 0.2]

        # creation of the Matplotlib frame
        self.figure = matplotlib.figure.Figure(figsize=(2, 2))
        self.plotax = self.figure.add_axes(self.plotaxcoords)
        self.cbarax = self.figure.add_axes(self.cbaraxcoords)
        self.cbarax.axis('off')
        self.canvas = FigureCanvas(self.figure, master=gui)

        # list of colormaps that can be used
        self.cmapnames = [m for m in plt.cm.datad]
        self.cmapnames.sort()
        # colorbar limits
        self.clim = None

        # settings of the basemap object
        self.bmap = None
        self.latn = 90
        self.lats = -90
        self.lonw = -180
        self.lone = 180
        self.res = 'l'
        self.lon0 = 270
        self.blat = 30
        self.proj = 'cyl'
        self.lonbg = None
        self.latbg = None
        self.mapbg = None

        # a the section edition
        self.secname = None
        self.indsec = None
        self.segment = None
        self.segdir = None
        self.indpoint = None

        # attributes used in the creation of new sections
        self.newoutx = []  # x coordinates
        self.newouty = []  # y coordinates
        self.neworient = []  # orientation
        self.numbernew = 1  # number of new sections (used in the name)
        self.cptnew = 0  # Number of points added

        # initialisation the matplotlib events
        self.secrelease = None
        self.secmove = None
        self.newsections = None
        self.editsections = None
        self.pointmove = None
        self.pointrelease = None

        # initialisation of the stuff used in the
        # matplotlib events
        self.deltax = None
        self.deltay = None
        self.xclick = None
        self.yclick = None
        self.scatter = None

    def init_bmap(self):

        """
        This function initialised the basemap object, when the coordinates,
        resolution and projections are changed.
        """

        if str(self.proj) == 'cyl':  # cylindrical projection: parameters are domain limits
            self.bmap = Basemap(llcrnrlon=self.lonw, urcrnrlon=self.lone, llcrnrlat=self.lats,
                                urcrnrlat=self.latn, projection=self.proj, ax=self.plotax,
                                resolution=self.gui.combobox_res.get(), suppress_ticks=False)

        elif str(self.proj) == 'npstere':  # northern hemisphere projection
            self.bmap = Basemap(lon_0=self.lon0, boundinglat=90-self.blat,
                                projection=self.proj, ax=self.plotax,
                                resolution=self.gui.combobox_res.get())

        elif str(self.proj) == 'spstere':  # southern hemisphere projection
            self.bmap = Basemap(lon_0=self.lon0, boundinglat=-90+self.blat,
                                projection=self.proj, ax=self.plotax,
                                resolution=self.gui.combobox_res.get())

        else:  # this is the Lambert Conformal Projection (not masked)
            self.make_lambert()

        self.init_plot()
        self.canvas.draw()

    def init_plot(self):

        """
        This function handles the plot of the background map
        (filled continents, etopo or background file).

        The sections are drawn here by using the :py:func:`draw_sections` function
        """

        self.figure.clf()  # we first clean the figure
        self.plotax = self.figure.add_axes(self.plotaxcoords)  # we redefine the map axes
        self.cbarax = self.figure.add_axes(self.cbaraxcoords)  # and the colormap axes
        self.cbarax.axis('off')  # we "hide" the axis of the colormap

        if self.gui.combobox_mode.get() == 'ETOPO':
            self.bmap.etopo(zorder=-1000, ax=self.plotax)
            self.bmap.drawcoastlines(color='k', linewidth=1, ax=self.plotax)

        elif self.gui.combobox_mode.get() == 'Filled continents':
            self.bmap.drawcoastlines(color='k', linewidth=1, zorder=-1000, ax=self.plotax)
            self.bmap.fillcontinents(color='DarkGray', zorder=-1001, ax=self.plotax)

        if self.gui.combobox_mode.get() == 'Map Background':

            if 1:#try:

                xmap, ymap = self.bmap(self.lonbg, self.latbg)

                if self.clim is None:
                    cmin = self.mapbg.min()
                    cmax = self.mapbg.max()
                    self.clim = str(cmin) + ',' + str(cmax)
                    self.gui.entry_clim_lab.set(self.clim)
                    levels = np.linspace(cmin, cmax, 21)

                else:
                    stout = self.gui.entry_clim_lab.get().split(',')
                    cmin = float(stout[0])
                    cmax = float(stout[1])
                    levels = np.linspace(cmin, cmax, 21)

                cscont = self.bmap.contourf(xmap, ymap, self.mapbg, levels=levels,
                                            ax=self.plotax, cmap=self.gui.combobox_cmap.get(),
                                            extend='both')

                try:
                    plt.colorbar(cscont, ax=self.cbarax, orientation='horizontal')
                except:
                    pass

                self.bmap.drawcoastlines(color='k', linewidth=1, ax=self.plotax)

            else:#except:
                tkMessageBox.showinfo('Oups!', 'Problem with the map background. You might need to check the NetCDF file')  # pylint: disable=line-too-long
                self.gui.combobox_mode.set('Filled continents')
                self.init_plot()

        self.draw_sections()

    def draw_sections(self):

        """
        Draw the sections that has been validated (the black lines).

        If we are in edition section mode, we activate the right mpl_event.

        If we create new sections, we keep the event deactivated.

        At the end, the matplotlib selection event is connected.
        At the beginning, if connected, we disconnect it.
        """

        try:
            self.figure.canvas.mpl_disconnect(self.editsections)
        except:
            pass

        try:
            self.figure.canvas.mpl_disconnect(self.newsections)
        except:
            pass

        self.plotax.lines = []  # lines of the ax instance
        self.plotax.texts = []  # texts of the ax instance

        # we run all the sections
        for indice_sec in range(0, len(self.sections)):
            section_int = self.sections[indice_sec]
            # x,y coords of the sections
            x_section, y_section = self.bmap(section_int.lon, section_int.lat)
            self.bmap.plot(x_section, y_section, picker=1, label=section_int.name,
                           color='Black', ax=self.plotax)  # we plot the sections

        if self.secname is not None:  # if one section is selected
            section_int = self.sections[self.indsec]
            # we note the position of the line nodes
            x_section_points = self.plotax.lines[self.indsec].get_xdata()
            y_section_points = self.plotax.lines[self.indsec].get_ydata()
            
            # we add a scatter plot for the edges, we set a picker of 5
            self.scatter = self.bmap.scatter(x_section_points, y_section_points,
                                             50, marker='o', picker=5,
                                             color='DarkOrange',
                                             zorder=10000, ax=self.plotax)

            # we add the text of the segment
            for indice_segment in range(0, len(section_int.lon)-1):
                x_text = 0.5*(x_section_points[indice_segment] + x_section_points[indice_segment+1])
                y_text = 0.5*(y_section_points[indice_segment] + y_section_points[indice_segment+1])
                seg_text = str(indice_segment+1)
                self.plotax.text(x_text, y_text, seg_text,
                                 bbox=dict(boxstyle="round", fc="1"),
                                 va='center', ha='center', fontsize=10, picker=5)

            if self.segdir is not None:  # if a segment is selected, we draw the segment in red
                self.plotax.texts[self.segment].set_color('r')

        # if section edition, we activate the secedit event
        if not self.gui.radiobox_editmode_lab.get():
            # we connect the event
            self.editsections = self.figure.canvas.mpl_connect('pick_event',
                                                               self.on_secedit_mplevent)
        else:
            self.newsections = self.figure.canvas.mpl_connect('button_press_event',
                                                              self.on_newsections_mplevent)

        self.draw_par_med()  # we draw the lines of the meridians/parallels
        self.canvas.draw()  # we draw the canvas

    def draw_par_med(self):

        """ Function that handles the drawing of the parallels/meridians """

        if self.proj == 'cyl':
            self.bmap.drawparallels(np.linspace(self.lats, self.latn, 11), ax=self.plotax)
            self.bmap.drawmeridians(np.linspace(self.lonw, self.lone, 11), ax=self.plotax)

        elif self.proj in ['npstere', 'spstere']:
            if self.proj == 'npstere':
                self.bmap.drawparallels(np.arange(90-self.blat, 90, 2), ax=self.plotax)
            else:
                self.bmap.drawparallels(-np.arange(90-self.blat, 90, 2), ax=self.plotax)
            self.bmap.drawmeridians(np.linspace(-180, 180, 11), ax=self.plotax)
        else:
            self.bmap.drawparallels(np.linspace(self.lats, self.latn, 11), ax=self.plotax)
            self.bmap.drawmeridians(np.linspace(self.lonw, self.lone, 11), ax=self.plotax)

        self.canvas.draw()

    def clean_selected_section(self):

        """
        Function that is called when the selected section is deselected.

        The section edition widgets are masked, the self.secname attribute is
        set to None
        """

        self.gui.entry_loncord.configure(state='disabled')
        self.gui.entry_latcord.configure(state='disabled')
        self.gui.entry_loncord_lab.set('')
        self.gui.entry_latcord_lab.set('')

        if self.secname is not None:
            self.secname = None  # we set the secname to None
            self.scatter.remove()  # we remove the scatter
            self.draw_sections()  # we redraw the section
            self.gui.combobox_seg.configure(state='disabled')  # we hide the "seg dir" txtctr widget
            self.gui.combobox_seg.set('')
            self.gui.label_seg_lab.set('Segment:')  # we change the "segdir" status
            self.gui.entry_secname_lab.set('')  # we clean the "secname" txtcrt
            self.gui.entry_secname.configure(state='disabled')  # we disable it
            self.gui.button_del.configure(state='disabled')

    def on_newsections_mplevent(self, event):   # pylint: disable=unused-argument

        """
        This is the matplotlib event that is used in the section creation.

        It is an event that is driven by mouse clicks.
        """

        if event.button == 1 and event.inaxes:  # if left click and in the axes
            self.newoutx.append(event.xdata)  # append the x data coordinate to the list
            self.newouty.append(event.ydata)  # append the y data coordinate to the list
            self.cptnew = self.cptnew+1  # number of points of the section
            # Plot the clicking coordinate
            self.plotax.plot([event.xdata], [event.ydata], color='Magenta', marker='o')
            self.canvas.draw()  # force drawing

            if self.cptnew >= 2:  # if more than two points (i.e. a section)
                # Draw the line section
                self.plotax.plot(self.newoutx[-2:], self.newouty[-2:], color='Magenta')
                self.canvas.draw()

        elif (event.button == 2) | (event.button == 3): 
            # if left click and two points have been defined, creation of a new section
            if len(self.newoutx) >= 2:
                lonnew, latnew = self.bmap(self.newoutx, self.newouty, inverse=True)
                self.sections.append(pypago.sections.Section('section'+str(self.numbernew),
                                                             lonnew, latnew, ['NE']*(len(self.newoutx)-1)))

            # reinitialisation of the section variables
            self.numbernew = self.numbernew + 1
            self.clean_edit_sections()

    def on_secedit_mplevent(self, event):    # pylint: disable=unused-argument

        """
        Matplotlib event for the section edition

        Based on mouse events
        """

        # Coordinates of the click event
        self.xclick = event.mouseevent.xdata
        self.yclick = event.mouseevent.ydata

        if isinstance(event.artist, matplotlib.collections.PathCollection):
        
            # If click on the scatter plot (orange points)
            # we disconnect the mpl event of section move and section release
            self.figure.canvas.mpl_disconnect(self.secmove)
            self.figure.canvas.mpl_disconnect(self.secrelease)

            self.indpoint = event.ind  # we recover the index of the point that was selected
            # we recover the coordinates of the point
            self.xclick, self.yclick = self.bmap(self.sections[self.indsec].lon[self.indpoint],
                                                 self.sections[self.indsec].lat[self.indpoint])

            # we activate the point move and point release events
            self.pointmove = self.figure.canvas.mpl_connect('motion_notify_event',
                                                            self.on_pointmove)
            self.pointrelease = self.figure.canvas.mpl_connect('button_release_event',
                                                               self.on_pointrelease)

        elif isinstance(event.artist, matplotlib.text.Text):

            # If click on text, we disconnect the mpl event of section move and section release

            self.figure.canvas.mpl_disconnect(self.secmove)
            self.figure.canvas.mpl_disconnect(self.secrelease)

            # here, we change the color of the previously
            if self.segdir is not None:
                self.plotax.texts[self.segment].set_color('k')

            # we retrieve the number of the segment
            self.segment = int(event.artist.get_text())-1
            # we retrieve the orientation of the segment
            self.segdir = self.sections[self.indsec].dire[self.segment]
            # we draw the text in red
            self.plotax.texts[self.segment].set_color('r')
            # we activate the wid_seg_dir widget
            self.gui.combobox_seg.config(state='normal')
            # we force its value to segdir
            self.gui.combobox_seg.set(self.segdir)
            # we modify the label of the segment status
            self.gui.label_seg_lab.set('Segment #' + str(self.segment+1)+':')
            self.canvas.draw()

        elif isinstance(event.artist, matplotlib.lines.Line2D):

            # If we click on a line

            # we first clear the segment dir
            self.segdir = None
            self.gui.combobox_seg.set('')
            self.gui.combobox_seg.config(state='disabled')

            # if a section was selected first
            # we remove the texts
            # we remove the scatter plot
            if self.secname is not None:
                self.plotax.texts = []
                self.scatter.remove()

            # we get the name of the line (the name of the section)
            self.secname = event.artist.get_label()

            # we enable the secname widget
            self.gui.entry_secname.config(state='normal')

            # we set the value of the text to the name of the section
            self.gui.entry_secname_lab.set(self.secname)

            self.gui.entry_loncord.configure(state='normal')
            self.gui.entry_latcord.configure(state='normal')

            # we retrieve the index of the section
            self.indsec = pypago.misc.findsecnum(self.sections, self.secname)
            secint = self.sections[self.indsec]

            # we retrieve the coordinates
            lonlon = ','.join(np.round(np.array(secint.lon), 3).astype(np.str))
            latlat = ','.join(np.round(np.array(secint.lat), 3).astype(np.str))

            # we add the coordinates to the txtctrl widget
            self.gui.entry_loncord_lab.set(lonlon)
            self.gui.entry_latcord_lab.set(latlat)
            self.gui.button_del.configure(state='normal')

            # we draw the scatter plot on all the points of the selected section
            self.scatter = self.bmap.scatter(self.plotax.lines[self.indsec].get_xdata(),
                                             self.plotax.lines[self.indsec].get_ydata(),
                                             50, marker='o', picker=5,
                                             color='DarkOrange', zorder=10000, ax=self.plotax)

            x_section, y_section = self.bmap(secint.lon, secint.lat)  # coordinates of the section

            # we draw the text of segments
            for indice_seg in range(0, len(secint.lon)-1):
                x_text = 0.5*(x_section[indice_seg] + x_section[indice_seg+1])
                y_text = 0.5*(y_section[indice_seg] + y_section[indice_seg+1])
                seg_text = str(indice_seg+1)
                self.plotax.text(x_text, y_text, seg_text, bbox=dict(boxstyle="round", fc="1"),
                                 va='center', ha='center', fontsize=10, picker=5)

            self.canvas.draw()

            # we activate the section move and section release events
            #self.secmove = self.figure.canvas.mpl_connect('motion_notify_event',
            #                                              self.on_secmove)
            #self.secrelease = self.figure.canvas.mpl_connect('button_release_event',
            #                                                 self.on_secrelease)

    #def on_secmove(self, event):    # pylint: disable=unused-argument

    #    """
    #    Event that handles when keep holding on a section
    #    """

    #    # this is the section that has been selected
    #    section_int = self.sections[self.indsec]
    #    # this is the map coordinates of the section
    #    x_section, y_section = self.bmap(section_int.lon, section_int.lat)

    #    if event.inaxes:  # if we move on the axis
    #        self.deltax = event.xdata - self.xclick  # this is the x displacement
    #        self.deltay = event.ydata - self.yclick  # this is the y displacement
    #        # we change the position of the line
    #        self.plotax.lines[self.indsec].set_xdata(x_section + self.deltax)
    #        self.plotax.lines[self.indsec].set_ydata(y_section + self.deltay)

    #        # here we change the position of the text
    #        for indice_seg in range(0, len(section_int.lon)-1):
    #            x_text = 0.5*(x_section[indice_seg] + x_section[indice_seg+1]) + self.deltax
    #            y_text = 0.5*(y_section[indice_seg] + y_section[indice_seg+1]) + self.deltay
    #            self.plotax.texts[indice_seg].set_position([x_text, y_text])

    #        # we remove and redraw the scatter plot
    #        self.scatter.remove()
    #        self.scatter = self.bmap.scatter(x_section + self.deltax, y_section + self.deltay,
    #                                         50, marker='o', picker=5,
    #                                         color='DarkOrange', zorder=10000, ax=self.plotax)

    #        # we force the drawing
    #        self.canvas.draw()

    #def on_secrelease(self, event):    # pylint: disable=unused-argument

    #    """
    #    Event that handles when we release a click on a line
    #    """

    #    # First we disconnect the two events relative to the section
    #    self.figure.canvas.mpl_disconnect(self.secrelease)
    #    self.figure.canvas.mpl_disconnect(self.secmove)

    #    # we recover the lon/lat coordinates of the section
    #    lonout, latout = self.bmap(self.plotax.lines[self.indsec].get_xdata(),
    #                               self.plotax.lines[self.indsec].get_ydata(),
    #                               inverse=True)

    #    # we change the lon/lat coordinates of the section object
    #    self.sections[self.indsec].lon = lonout
    #    self.sections[self.indsec].lat = latout

    #    # We change the lon/lat coordinates of the txtctr widgets
    #    lonlon = ','.join(np.round(np.array(lonout), 3).astype(np.str))
    #    latlat = ','.join(np.round(np.array(latout), 3).astype(np.str))
    #    self.gui.entry_loncord_lab.set(lonlon)
    #    self.gui.entry_latcord_lab.set(latlat)

    def on_pointmove(self, event):    # pylint: disable=unused-argument

        """
        Event that handles when we keep the mouse on a scatter point
        """

        secint = self.sections[self.indsec]  # we recover the section that is selected

        self.deltax = event.xdata - self.xclick  # this is the displacement
        self.deltay = event.ydata - self.yclick

        xint, yint = self.bmap(secint.lon, secint.lat)  # we recover the section coordinates

        # this is the new X coordinate of the point
        xxx = np.empty(secint.lon.shape)
        xxx[:] = xint
        xxx[self.indpoint] = xxx[self.indpoint] + self.deltax

        # this is the new Y coordinate of the point
        yyy = np.empty(secint.lat.shape)
        yyy[:] = yint
        yyy[self.indpoint] = yyy[self.indpoint] + self.deltay

        # we change the point coordinate on the line
        self.plotax.lines[self.indsec].set_xdata(xxx)
        self.plotax.lines[self.indsec].set_ydata(yyy)

        # we remove and redo the scatter plot
        self.scatter.remove()
        self.scatter = self.bmap.scatter(xxx, yyy, 50, marker='o',
                                         picker=5, color='DarkOrange',
                                         zorder=10000, ax=self.plotax)

        # we add the segment text
        for indice_seg in range(0, len(secint.lon)-1):
            x_text = 0.5*(xxx[indice_seg] + xxx[indice_seg+1])
            y_text = 0.5*(yyy[indice_seg] + yyy[indice_seg+1])
            self.plotax.texts[indice_seg].set_position([x_text, y_text])

        self.canvas.draw()

    def on_pointrelease(self, event):    # pylint: disable=unused-argument

        """
        Event that handles when we release the click on a point
        """

        # We disconnect the point events
        self.figure.canvas.mpl_disconnect(self.pointmove)
        self.figure.canvas.mpl_disconnect(self.pointrelease)

        # we recover the lon/lat coordinates of the new line section
        lonout, latout = self.bmap(self.plotax.lines[self.indsec].get_xdata(),
                                   self.plotax.lines[self.indsec].get_ydata(),
                                   inverse=True)

        # we change the lon/lat of the section
        self.sections[self.indsec].lon = lonout
        self.sections[self.indsec].lat = latout

        # We change the lon/lat coordinates of the txtctr widgets
        lonlon = ','.join(np.round(np.array(lonout), 3).astype(np.str))
        latlat = ','.join(np.round(np.array(latout), 3).astype(np.str))
        self.gui.entry_loncord_lab.set(lonlon)
        self.gui.entry_latcord_lab.set(latlat)

    def plot_section_newcoord(self):

        """
        Module that handle the plotting of the sections when
        new coordinates have been enterred manually by the user
        """

        if self.gui.radiobox_editmode_lab.get():  # creation of a new section

            try:
                # we recover the variables of lon,lat
                lonint = np.array(self.gui.entry_loncord_lab.get().split(',')).astype(np.float)
                latint = np.array(self.gui.entry_latcord_lab.get().split(',')).astype(np.float)
                if len(lonint) == len(latint):  # if both have the same sizes
                    # we add our new section
                    self.sections.append(pypago.sections.Section('section'+str(self.numbernew),
                                                                 lonint, latint, ['NE']*(len(lonint)-1)))
                    self.numbernew = self.numbernew+1  # we iterate our number of new sections
                    x_section, y_section = self.bmap(lonint, latint)  # we plot this section
                    # Plot the clicking coordinate
                    self.plotax.plot(x_section, y_section, color='Magenta', marker='o')
                    self.figure.canvas.draw()  # force drawing
                else:
                    tkMessageBox.showinfo('Oups!', 'Section coordinates must have the same size')
            except:
                pass

        else:  # we edit a section

            try:

                lonint = np.array(self.gui.entry_loncord_lab.get().split(',')).astype(np.float)
                latint = np.array(self.gui.entry_latcord_lab.get().split(',')).astype(np.float)

                # if the new coords has the same length that old ones: keep the same
                # directions
                if (len(self.sections[self.indsec].lon) == len(lonint)) & \
                   (len(self.sections[self.indsec].lat) == len(latint)):
                    self.sections[self.indsec].lon = lonint
                    self.sections[self.indsec].lat = latint
                    self.scatter.remove()
                    self.draw_sections()
                else:

                    if len(lonint) == len(latint):
                        tkMessageBox.showinfo('Oups!',
                                              "Old and new coordinates have not the same size.\n" +
                                              "The segments's directions have been set to 'NE'")
                        self.sections[self.indsec].lon = lonint
                        self.sections[self.indsec].lat = latint
                        self.sections[self.indsec].dire = ['NE']*(len(lonint)-1)
                        self.scatter.remove()
                        self.segdir = None
                        self.gui.combobox_seg.configure(state='disabled')
                        self.gui.combobox_seg.set('')
                        self.draw_sections()

                    else:
                        tkMessageBox.showinfo('Oups!',
                                              'Section coordinates must have the same size')

            except:
                pass

    def clean_edit_sections(self):

        """
        This function resets the arrays which are used in the section creation.
        """

        self.newoutx = []
        self.newouty = []
        self.neworient = []
        self.cptnew = 0

    def make_lambert(self):

        """
        Function that creates a NCL/Matlab like
        lambert conformal projection (but without masking)

        """

        # We first recover the hemisphere
        if (self.lats >= 0) & (self.latn > 0):
            hem = 'NH'
        elif (self.latn <= 0) & (self.lats < 0):
            hem = 'SH'
        else:
            hem = None

        if hem is not None:  # if hemisphere is not none, the lcc projection can be used is

            # Here, we convert the Python lcc to a NCL-like lcc (cf. my script lambert.py)
            if hem == 'NH':
                lat2 = 89.999
                lat1 = 0.001
            else:
                lat2 = -89.999
                lat1 = -0.001
            lon0 = 0.5*(self.lone + self.lonw)
            self.bmap = Basemap(llcrnrlon=self.lonw, llcrnrlat=self.lats,
                                urcrnrlon=self.lone, urcrnrlat=self.latn,
                                lon_0=lon0, lat_1=lat1, lat_2=lat2,
                                projection=self.proj, resolution='c')

            if hem == 'NH':
                xmin = self.bmap(self.lonw, self.lats)[0]
                xmax = self.bmap(self.lone, self.lats)[0]
                ymin = self.bmap(lon0, self.lats)[1]
                ymax = self.bmap(self.lonw, self.latn)[1]
            else:
                xmin = self.bmap(self.lonw, self.latn)[0]
                xmax = self.bmap(self.lone, self.latn)[0]
                ymax = self.bmap(lon0, self.latn)[1]
                ymin = self.bmap(self.lonw, self.lats)[1]

            lonmin, latmin = self.bmap(xmin, ymin, inverse=True)
            lonmax, latmax = self.bmap(xmax, ymax, inverse=True)

            self.bmap = Basemap(llcrnrlon=lonmin, llcrnrlat=latmin,
                                urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon0,
                                lat_1=lat1, lat_2=lat2, projection=self.proj,
                                resolution=self.gui.combobox_res.get(), ax=self.plotax)

        else:
            # If our domain crosses the equator, lcc cannot be used
            # -> we change the projection to cyl
            self.proj = 'cyl'
            self.gui.combobox_proj.set(self.proj)
            tkMessageBox.showinfo('Oups!', "'To use lcc projection, your domain must not cross \
            the equator\nBack to the cyl projection.'")
            self.init_bmap()
