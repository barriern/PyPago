#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module that handles the GUI for the definition of section endpoints
"""

from __future__ import print_function

import sys
import os
import Tkinter as tk
import tkFileDialog
import ttk

import pypago.pyio
import pypago.misc
from pypago.guis.NewToolbar import NewToolbar
from pypago.guis._matplotlib_editions_sections import MatplotlibEditionsSections

class EditionSections(tk.Tk):

    """
    Class that defines the GUIs associated with the program for the edition of section endpoints
    Inherits from the :py:class:`Tkinter.Tk` class.

    .. versionadded:: 2015-07-31
       barrier.n

    """

    def __init__(self, parent):

        """

        Initialisation of the program window.

        :param parent:

        """

        tk.Tk.__init__(self, parent)
        self.minsize(height=700, width=900)  # set the minimum size of the window
        self.parent = parent

        self.bmapplot = MatplotlibEditionsSections(self)

        self.backgroud_file_name = None
        self.section_file_name = None  # name of the section file

        self.make_layout()  # add the widgets to the window
        self.bind_functions()  # bind the widgets to their functions
        self.make_menu()  # create the menus and their functions

        self.bmapplot.init_bmap()  # initialize the drawing of the map background

    def make_layout(self):

        """
        In this function, all the widgets are created.

        The layout is managed by the :py:func:`Tkinter.grid` function
        """

        self.grid()

        self.maphandling = tk.Label(self, text='Map handling', bg='black', fg='white')
        self.maphandling.grid(row=0, column=0, columnspan=4, sticky='EW')
        self.sechandling = tk.Label(self, text='Section handling', bg='navyblue', fg='white')
        self.sechandling.grid(row=0, column=4, columnspan=4, sticky='EW')

        self.label_lonw = tk.Label(self, text='lonw')
        self.label_lonw.grid(row=1, column=0)
        self.entry_lonw_lab = tk.StringVar()
        self.entry_lonw = tk.Entry(self, textvariable=self.entry_lonw_lab)
        self.entry_lonw.grid(row=1, column=1)
        self.entry_lonw_lab.set(self.bmapplot.lonw)

        self.label_lone = tk.Label(self, text='lone', anchor="w")
        self.label_lone.grid(row=1, column=2)
        self.entry_lone_lab = tk.StringVar()
        self.entry_lone = tk.Entry(self, textvariable=self.entry_lone_lab)
        self.entry_lone.grid(row=1, column=3)
        self.entry_lone_lab.set(self.bmapplot.lone)

        self.label_lats = tk.Label(self, text='lats', anchor="w")
        self.label_lats.grid(row=2, column=0)
        self.entry_lats_lab = tk.StringVar()
        self.entry_lats = tk.Entry(self, textvariable=self.entry_lats_lab)
        self.entry_lats.grid(row=2, column=1)
        self.entry_lats_lab.set(self.bmapplot.lats)

        self.label_latn = tk.Label(self, text='latn', anchor="w")
        self.label_latn.grid(row=2, column=2)
        self.entry_latn_lab = tk.StringVar()
        self.entry_latn = tk.Entry(self, textvariable=self.entry_latn_lab)
        self.entry_latn.grid(row=2, column=3)
        self.entry_latn_lab.set(self.bmapplot.latn)

        self.label_lon0 = tk.Label(self, text='lon0', anchor="w")
        self.label_lon0.grid(row=3, column=0)
        self.entry_lon0_lab = tk.StringVar()
        self.entry_lon0 = tk.Entry(self, textvariable=self.entry_lon0_lab)
        self.entry_lon0.grid(row=3, column=1)
        self.entry_lon0_lab.set(self.bmapplot.lon0)

        self.label_blat = tk.Label(self, text='Bounding lat.', anchor="w")
        self.label_blat.grid(row=3, column=2)
        self.entry_blat_lab = tk.StringVar()
        self.entry_blat = tk.Entry(self, textvariable=self.entry_blat_lab)
        self.entry_blat.grid(row=3, column=3)
        self.entry_blat_lab.set(self.bmapplot.blat)

        self.label_proj = tk.Label(self, text='Projection', anchor="w")
        self.label_proj.grid(row=4, column=0)
        self.combobox_proj = ttk.Combobox(self, values=['cyl', 'npstere', 'spstere', 'lcc'])
        self.combobox_proj.grid(row=4, column=1)
        self.combobox_proj.set('cyl')

        self.label_res = tk.Label(self, text='Resolution', anchor="w")
        self.label_res.grid(row=4, column=2)
        self.combobox_res = ttk.Combobox(self, values=['c', 'l', 'i', 'h', 'f'])
        self.combobox_res.grid(row=4, column=3)
        self.combobox_res.set('c')

        self.label_mode = tk.Label(self, text='Plot mode', anchor="w")
        self.label_mode.grid(row=5, column=0)
        self.combobox_mode = ttk.Combobox(self,
                                          values=['Filled continents',
                                                  'ETOPO',
                                                  'Map Background'])
        self.combobox_mode.grid(row=5, column=1)
        self.combobox_mode.set('Filled continents')

        self.label_cmap = tk.Label(self, text='Colormap', anchor="w")
        self.label_cmap.grid(row=5, column=2)
        self.combobox_cmap = ttk.Combobox(self, values=self.bmapplot.cmapnames)
        self.combobox_cmap.grid(row=5, column=3)
        self.combobox_cmap.set('jet')

        self.label_clim = tk.Label(self, text='Colormap lim.', anchor="w")
        self.label_clim.grid(row=6, column=0)
        self.entry_clim_lab = tk.StringVar()
        self.entry_clim = tk.Entry(self, textvariable=self.entry_clim_lab)
        self.entry_clim.grid(row=6, column=1)
        self.entry_clim_lab.set('')

        self.label_secname = tk.Label(self, text='Section name', anchor="w")
        self.label_secname.grid(row=1, column=4)

        self.entry_secname_lab = tk.StringVar()
        self.entry_secname = tk.Entry(self, textvariable=self.entry_secname_lab)
        self.entry_secname.grid(row=1, column=5, columnspan=2, sticky='EW')
        self.entry_secname_lab.set('')
        self.entry_secname.configure(state='disabled')

        self.label_seg_lab = tk.StringVar()
        self.label_seg = tk.Label(self, textvariable=self.label_seg_lab, anchor="w")
        self.label_seg.grid(row=2, column=4)
        self.label_seg_lab.set('Segment:')
        self.combobox_seg = ttk.Combobox(self, values=['NE', 'NW', 'SE', 'SW'])
        self.combobox_seg.grid(row=2, column=5, columnspan=2, sticky='EW')
        self.combobox_seg.configure(state='disabled')

        self.label_loncord = tk.Label(self, text='Lon. coord.', anchor="w")
        self.label_loncord.grid(row=3, column=4)
        self.entry_loncord_lab = tk.StringVar()
        self.entry_loncord = tk.Entry(self, textvariable=self.entry_loncord_lab)
        self.entry_loncord.grid(row=3, column=5, columnspan=2, sticky='EW')
        self.entry_loncord_lab.set('')
        self.entry_loncord.configure(state='disabled')

        self.label_latcord = tk.Label(self, text='Lat. coord.', anchor="w")
        self.label_latcord.grid(row=4, column=4)
        self.entry_latcord_lab = tk.StringVar()
        self.entry_latcord = tk.Entry(self, textvariable=self.entry_latcord_lab)
        self.entry_latcord.grid(row=4, column=5, columnspan=2, sticky='EW')
        self.entry_latcord_lab.set('')
        self.entry_latcord.configure(state='disabled')

        self.button_del = tk.Button(self, text='Del. section', command=self.on_delete_button)
        self.button_del.grid(row=5, column=4)
        self.button_del.configure(state='disabled')

        self.label_secmode = tk.Label(self, text='Section mode', anchor="w")
        self.label_secmode.grid(row=6, column=4)

        self.radiobox_editmode_lab = tk.IntVar()
        self.radiobox_editmode = tk.Radiobutton(self, text='Edit',
                                                variable=self.radiobox_editmode_lab, value=0,
                                                command=self.on_change_working_mode)
        self.radiobox_editmode.grid(row=6, column=5)
        self.radiobox_editmode = tk.Radiobutton(self, text='Add',
                                                variable=self.radiobox_editmode_lab, value=1,
                                                command=self.on_change_working_mode)
        self.radiobox_editmode.grid(row=6, column=6)
        self.radiobox_editmode_lab.set(0)

        self.bmapplot.canvas.get_tk_widget().grid(row=7, column=0, columnspan=8, sticky='NSEW')

        # Creation of the Matplotlib toolbar
        toolbar_frame = tk.Frame(self)
        self.toolbar = NewToolbar(self.bmapplot.canvas, toolbar_frame)
        toolbar_frame.grid(row=8, column=0, columnspan=8)

        # Management of the resizable parts of the grid
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        self.grid_columnconfigure(3, weight=1)
        self.grid_columnconfigure(4, weight=1)
        self.grid_columnconfigure(5, weight=1)
        self.grid_columnconfigure(6, weight=1)
        self.grid_rowconfigure(7, weight=1)

        self.entry_clim.configure(state='disabled')
        self.combobox_cmap.configure(state='disabled')
        self.label_clim.configure(state='disabled')
        self.label_cmap.configure(state='disabled')

    def bind_functions(self):

        """
        Function that binds the widgets to their respective functions
        """

        self.entry_lonw.bind("<Return>", self.on_change_domain)
        self.entry_lone.bind("<Return>", self.on_change_domain)
        self.entry_lats.bind("<Return>", self.on_change_domain)
        self.entry_latn.bind("<Return>", self.on_change_domain)
        self.entry_lon0.bind("<Return>", self.on_change_domain)
        self.entry_blat.bind("<Return>", self.on_change_domain)

        self.combobox_res.bind("<<ComboboxSelected>>", self.on_change_basemap)
        self.combobox_proj.bind("<<ComboboxSelected>>", self.on_change_basemap)
        self.combobox_mode.bind("<<ComboboxSelected>>", self.on_change_plot)

        # Change in clim -> change the clim attribute, redraw
        self.entry_clim.bind("<Return>", self.on_change_clim)

        # Change in colormap -> only need to redraw
        self.combobox_cmap.bind("<<ComboboxSelected>>", self.on_change_plot)

        # Function associated with change in direction
        self.combobox_seg.bind("<<ComboboxSelected>>", self.on_change_dir)

        # Function associated with a change in the section name
        self.entry_secname.bind('<Return>', self.on_sec_name_change)

        # Functions associated with change in section coordinates
        self.entry_loncord.bind('<Return>', self.on_change_sec_coords)
        self.entry_latcord.bind('<Return>', self.on_change_sec_coords)

    def make_menu(self):

        """
        Function that creates the Menu and that makes the link
        with their respective functions
        """

        menubar = tk.Menu(self)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open", command=self.on_open)
        file_menu.add_command(label="Save", command=self.on_save)
        file_menu.add_command(label="Save As", command=self.on_save_as)
        file_menu.add_command(label="Quit", command=self.on_quit)
        menubar.add_cascade(label="File", menu=file_menu)

        background_menu = tk.Menu(menubar, tearoff=0)  # create a background menu
        # Add a load background item
        background_menu.add_command(label='Load background',
                                    command=self.on_background)
        menubar.add_cascade(label="Background", menu=background_menu)

        self.config(menu=menubar)

    # =============================================== Menu functions
    def on_open(self):

        """
        Function that is called when the Open menu is activated.
        Allows to open a :file:`.pygo` file that contains sections
        """

        options = {}
        options['defaultextension'] = '.pygo'
        options['filetypes'] = [('PyPago Files', '.pygo')]
        options['initialdir'] = '{0}/.pago'.format(os.environ['HOME'])
        options['title'] = u'Choose the .pygo file to be opened'
        filename = tkFileDialog.askopenfilename(**options)

        if filename:
            self.section_file_name = filename
            try:
                self.bmapplot.sections = pypago.pyio.load(self.section_file_name)
                self.bmapplot.numbernew = len(self.bmapplot.sections) + 1
            except:
                print('perdu')
                print("Unexpected error:", sys.exc_info()[0])
                raise

            self.radiobox_editmode_lab.set(0)
            self.bmapplot.clean_selected_section()
            self.on_change_working_mode()

    def on_save(self):

        """
        Function that is called when the Save menu is activated.

        If a filename is defined, it saves into the file.

        If no filename is defined, it prompts a file window.
        """

        if self.section_file_name is None:  # if the section_file_name is None -> SaveAs
            options = {}
            options['defaultextension'] = '.pygo'
            options['filetypes'] = [('PyPago Files', '.pygo')]
            options['initialdir'] = '{0}/.pago'.format(os.environ['HOME'])
            options['title'] = u'Choose the .pygo file to be saved'
            filename = tkFileDialog.asksaveasfilename(**options)
            if filename:
                self.section_file_name = filename
                pypago.pyio.save(self.bmapplot.sections, self.section_file_name)
        else:
            # if section file exists -> we save
            pypago.pyio.save(self.bmapplot.sections, self.section_file_name)

    def on_save_as(self):

        """
        Function that is called when the Save As menu is activated.

        It prompts a file window, in which the sections will be saved.
        """

        options = {}
        options['defaultextension'] = '.pygo'
        options['filetypes'] = [('PyPago Files', '.pygo')]
        options['initialdir'] = '{0}/.pago'.format(os.environ['HOME'])
        options['title'] = u'Choose the .pygo file to be saved'
        filename = tkFileDialog.asksaveasfilename(**options)
        if filename:
            self.section_file_name = filename
            pypago.pyio.save(self.bmapplot.sections, self.section_file_name)

    def on_quit(self):

        """
        Function that is called when the Quit menu is activated.

        It destroys the Tkinter frame
        """

        self.destroy()

    def on_background(self):

        """
        Function that is called when the Background menu is activated.

        It opens a file window and the user chooses a |netcdf| file.

        The :py:func:`pypago.pyio.read_bg_file` function is then called.

        The plot mode is switched to the background mode
        and the :py:func:`init_plot` function is called
        """

        options = {}
        options['defaultextension'] = '.nc'
        options['filetypes'] = [('NetCDF files', '.nc')]
        filename = tkFileDialog.askopenfilename(**options)

        if filename:
            self.backgroud_file_name = filename
            # we read the background file
            self.bmapplot.lonbg, self.bmapplot.latbg, self.bmapplot.mapbg = \
                pypago.pyio.read_bg_file(self.backgroud_file_name)
            self.combobox_mode.set('Map Background')  # we change the plotmode to "bg" mode
            self.bmapplot.clim = None  # we first reset the colorbar limits
            self.on_change_plot(None)
            self.bmapplot.lonw = self.bmapplot.lonbg.min() 
            self.bmapplot.lone = self.bmapplot.lonbg.max()
            self.bmapplot.lats = self.bmapplot.latbg.min()
            self.bmapplot.latn = self.bmapplot.latbg.max()
            self.entry_lonw_lab.set(self.bmapplot.lonw)
            self.entry_lone_lab.set(self.bmapplot.lone)
            self.entry_lats_lab.set(self.bmapplot.lats)
            self.entry_latn_lab.set(self.bmapplot.latn)
            self.bmapplot.init_bmap()
            self.bmapplot.init_plot()  # we redo the plot

    def on_change_working_mode(self):

        """
        Function that handles the change in the working mode (i.e. edition or
        creation of sections)

        The only way to edit the domain is in the edit section mode.
        Therefore, all the widgets for the map handling are deactivated if
        new sections are created
        """

        self.bmapplot.clean_selected_section()

        if self.radiobox_editmode_lab.get():  # if new sections
            self.entry_lonw.configure(state='disabled')
            self.entry_lone.configure(state='disabled')
            self.entry_latn.configure(state='disabled')
            self.entry_lats.configure(state='disabled')
            self.entry_lon0.configure(state='disabled')
            self.entry_blat.configure(state='disabled')
            self.combobox_proj.configure(state='disabled')
            self.combobox_res.configure(state='disabled')
            self.combobox_mode.configure(state='disabled')
            self.entry_loncord.configure(state='normal')
            self.entry_latcord.configure(state='normal')

        else:

            # barrier.n, 2016-09-19. Needed so that if the user edits a section
            # without a right-click (i.e without validation
            self.bmapplot.clean_edit_sections() 

            self.combobox_proj.configure(state='normal')
            self.combobox_res.configure(state='normal')
            self.combobox_mode.configure(state='normal')
            if self.bmapplot.proj == 'lcc' or self.bmapplot.proj == 'cyl':
                self.entry_lonw.configure(state='normal')
                self.entry_lone.configure(state='normal')
                self.entry_latn.configure(state='normal')
                self.entry_lats.configure(state='normal')
            else:
                self.entry_blat.configure(state='normal')
                self.entry_lon0.configure(state='normal')

        self.bmapplot.draw_sections()

    def on_change_dir(self, event):    # pylint: disable=unused-argument

        """
        Function that handles a change in the segment direction

        We change the value of the section self.indsec direction of segment
        self.segment
        """

        self.bmapplot.sections[self.bmapplot.indsec].dire[self.bmapplot.segment] = self.combobox_seg.get()

    def on_change_plot(self, event):    # pylint: disable=unused-argument

        """
        Function that is called when the map needs to be redrawn,
        but without having defined a new Basemap object
        """

        if self.combobox_mode.get() in ['Filled Continents', 'ETOPO']:
            self.entry_clim.configure(state='disabled')
            self.combobox_cmap.configure(state='disabled')
            self.label_clim.configure(state='disabled')
            self.label_cmap.configure(state='disabled')
        else:
            self.entry_clim.configure(state='normal')
            self.combobox_cmap.configure(state='normal')
            self.label_clim.configure(state='normal')
            self.label_cmap.configure(state='normal')
        self.bmapplot.init_plot()

    def on_change_basemap(self, event):    # pylint: disable=unused-argument

        """
        Function that is called when a new Basemap has been defined.

        The new basemap is defined and the map is redrawn.
        """

        self.bmapplot.proj = self.combobox_proj.get()  # we recover the projection

        if self.bmapplot.proj in ['cyl', 'lcc']:
            self.entry_lonw.configure(state='normal')
            self.entry_lone.configure(state='normal')
            self.entry_latn.configure(state='normal')
            self.entry_lats.configure(state='normal')
            self.entry_lon0.configure(state='disabled')
            self.entry_blat.configure(state='disabled')

        else:
            self.entry_lonw.configure(state='disabled')
            self.entry_lone.configure(state='disabled')
            self.entry_latn.configure(state='disabled')
            self.entry_lats.configure(state='disabled')
            self.entry_lon0.configure(state='normal')
            self.entry_blat.configure(state='normal')

        self.bmapplot.init_bmap()

    def on_delete_button(self):

        """
        Function that is called when a section is deleted
        """

        # if a section has been selected (i.e. the self.indsec is not None)
        if self.bmapplot.secname is not None:
            self.bmapplot.sections.pop(self.bmapplot.indsec)  # we remove the section
            self.bmapplot.clean_selected_section()

    def on_sec_name_change(self, event):    # pylint: disable=unused-argument

        """
        Function that is called when the section name is changed
        """

        # we change the name of the section
        self.bmapplot.sections[self.bmapplot.indsec].name = self.entry_secname_lab.get()
        # we change the label of the line that draws the plot
        self.bmapplot.plotax.lines[self.bmapplot.indsec].set_label(self.bmapplot.sections[self.bmapplot.indsec].name)

    def on_change_domain(self, event):    # pylint: disable=unused-argument

        """
        Function that is called when the domain limits (lonw, etc.) are changed

        If the entries can be converted into float, a new basemap is created
        and the map is redrawn
        """

        try:
            self.bmapplot.lonw = float(self.entry_lonw_lab.get())
            self.bmapplot.lone = float(self.entry_lone_lab.get())
            self.bmapplot.lats = float(self.entry_lats_lab.get())
            self.bmapplot.latn = float(self.entry_latn_lab.get())
            self.bmapplot.lon0 = float(self.entry_lon0_lab.get())
            self.bmapplot.blat = float(self.entry_blat_lab.get())
            self.bmapplot.init_bmap()
        except:
            self.entry_lonw_lab.set(self.bmapplot.lonw)
            self.entry_lone_lab.set(self.bmapplot.lone)
            self.entry_lats_lab.set(self.bmapplot.lats)
            self.entry_latn_lab.set(self.bmapplot.latn)
            self.entry_lon0_lab.set(self.bmapplot.lon0)
            self.entry_blat_lab.set(self.blat)

    def on_change_sec_coords(self, event):    # pylint: disable=unused-argument

        """
        Function that handles a change in the loncoords, latcoords entries.

        Two cases:

        - we edit a section
        - we create new section
        """

        self.bmapplot.plot_section_newcoord()

    def on_change_clim(self, event):    # pylint: disable=unused-argument

        """
        Function that is called when the clim has been changed
        """

        try:
            stout = self.entry_clim_lab.get().split(',')
            cmin = float(stout[0])
            cmax = float(stout[1])
            self.bmapplot.clim = str(cmin)+','+str(cmax)
            self.bmapplot.init_plot()
        except:
            self.entry_clim_lab.set(self.bmapplot.clim)

if __name__ == "__main__":
    PAGOFRAME = EditionSections(None)
    PAGOFRAME.title("PyPago Sections' editions")
    PAGOFRAME.mainloop()
    sys.exit(0)
