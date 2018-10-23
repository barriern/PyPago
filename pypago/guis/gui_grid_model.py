#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Class that defines the GUIs associated with the program for the extraction of
grid model and section staircases
"""

import pypago.pyio
import pypago.grid
import pypago.misc
import pypago.sections
import Tkinter as tk
import tkFileDialog
import ttk
import numpy as np
from pypago.guis.NewToolbar import NewToolbar
from pypago.guis._matplotlib_grid_model import MatplotlibGridModel
import sys
import os


class GridModel(tk.Tk):

    """
    Class that defines the GUIs associated with the program for the extraction of
    grid model and section staircases

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
        # set the minimum size of the window
        self.minsize(height=700, width=900)
        self.parent = parent

        # initialize the Matplotlib object associated with the GUI
        self.bmapplot = MatplotlibGridModel(self)

        # name of the meshgrid file
        self.gridfilename = None

        # name of the output file
        self.outputfilename = None

        # this is the name of the model that we process
        self.modelproc = 'NEMO'

        # creates the frame layout
        self.make_layout()

        self.reprocess_grid = False

        # creates the menu and its associated functions
        self.make_menu()

        # bind the widgets to their functions
        self.bind_functions()

    def make_layout(self):

        """
        In this function, all the widgets are created.

        The layout is managed by the :py:func:`Tkinter.grid` function
        """

        self.grid()

        # Description of the widgets (replace wx StaticBox)
        self.domainhandling = tk.Label(self, text='Domain handling', bg='black', fg='white')
        self.domainhandling.grid(row=0, column=0, columnspan=4, sticky='EW')
        self.sechandling = tk.Label(self, text='Section handling', bg='navyblue', fg='white')
        self.sechandling.grid(row=0, column=4, columnspan=4, sticky='EW')

        # Widgets for domain handling
        label_imin = tk.Label(self, text='imin')
        label_imin.grid(row=1, column=0)
        self.entry_imin_lab = tk.StringVar()
        self.entry_imin = tk.Entry(self, textvariable=self.entry_imin_lab)
        self.entry_imin.grid(row=1, column=1)

        label_imax = tk.Label(self, text='imax')
        label_imax.grid(row=1, column=2)
        self.entry_imax_lab = tk.StringVar()
        self.entry_imax = tk.Entry(self, textvariable=self.entry_imax_lab)
        self.entry_imax.grid(row=1, column=3)

        label_jmin = tk.Label(self, text='jmin')
        label_jmin.grid(row=2, column=0)
        self.entry_jmin_lab = tk.StringVar()
        self.entry_jmin = tk.Entry(self, textvariable=self.entry_jmin_lab)
        self.entry_jmin.grid(row=2, column=1)

        label_jmax = tk.Label(self, text='jmax')
        label_jmax.grid(row=2, column=2)
        self.entry_jmax_lab = tk.StringVar()
        self.entry_jmax = tk.Entry(self, textvariable=self.entry_jmax_lab)
        self.entry_jmax.grid(row=2, column=3)

        # Radiobox for working mode (edit domain, edit sections, check sections)
        self.radiobox_editmode_lab = tk.IntVar()

        self.radiobox_editmode0 = tk.Radiobutton(self, text='Edit domain',
                                                 variable=self.radiobox_editmode_lab, value=0,
                                                 command=self.on_change_status)
        self.radiobox_editmode0.grid(row=3, column=3)

        self.radiobox_editmode1 = tk.Radiobutton(self, text='Edit sections',
                                                 variable=self.radiobox_editmode_lab, value=1,
                                                 command=self.on_change_status)
        self.radiobox_editmode1.grid(row=4, column=3)

        self.radiobox_editmode2 = tk.Radiobutton(self, text='Check sections',
                                                 variable=self.radiobox_editmode_lab, value=2,
                                                 command=self.on_change_status)
        self.radiobox_editmode2.grid(row=5, column=3)

        # Widgets that handle the model information and selection
        label_model = tk.Label(self, text='Model:')
        label_model.grid(row=4, column=0)
        self.combobox_model = ttk.Combobox(self, values=['NEMO', 'GFDL'])
        self.combobox_model.set(self.modelproc)
        self.combobox_model.grid(row=4, column=1)

        # static text for secname
        self.label_secname_lab = tk.StringVar()
        self.label_secname = tk.Label(self, textvariable=self.label_secname_lab)
        self.label_secname.grid(row=1, column=4)
        self.label_secname_lab.set('Section name:')

        # static text for i coordinates
        self.label_i_lab = tk.StringVar()
        self.label_i = tk.Label(self, textvariable=self.label_i_lab)
        self.label_i.grid(row=2, column=4)
        self.label_i_lab.set('i coordinates:')

        # static text for j coordinates
        self.label_j_lab = tk.StringVar()
        self.label_j = tk.Label(self, textvariable=self.label_j_lab)
        self.label_j.grid(row=3, column=4)
        self.label_j_lab.set('j coordinates:')

        # static text for segment direction
        self.label_segment_lab = tk.StringVar()
        self.label_segment = tk.Label(self, textvariable=self.label_segment_lab)
        self.label_segment.grid(row=4, column=4)
        self.label_segment_lab.set('Segment:')

        # combobox for segment direction
        self.combobox_segment = ttk.Combobox(self, values=['NE', 'NW', 'SE', 'SW'])
        self.combobox_segment.grid(row=4, column=5)

        # button for section deletion
        self.button_del = tk.Button(self, text='Del. section', command=self.on_del_button)
        self.button_del.grid(row=5, column=4)

        # a matplotlib canvas is added
        self.bmapplot.canvas.get_tk_widget().grid(row=6, column=0, columnspan=8, sticky='NSEW')

        # Creation of the Matplotlib toolbar
        toolbar_frame = tk.Frame(self)
        self.toolbar = NewToolbar(self.bmapplot.canvas, toolbar_frame)
        toolbar_frame.grid(row=7, column=0, columnspan=8)

        # disabling of the widgets when the frame is started
        self.button_del.configure(state='disabled')
        self.combobox_segment.configure(state='disabled')
        self.entry_imin.configure(state='disabled')
        self.entry_imax.configure(state='disabled')
        self.entry_jmin.configure(state='disabled')
        self.entry_jmax.configure(state='disabled')
        self.radiobox_editmode0.configure(state='disabled')
        self.radiobox_editmode1.configure(state='disabled')
        self.radiobox_editmode2.configure(state='disabled')

        # frame reshape properties
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        self.grid_columnconfigure(3, weight=1)
        self.grid_columnconfigure(4, weight=1)
        self.grid_columnconfigure(5, weight=1)
        self.grid_rowconfigure(6, weight=1)

    def make_menu(self):

        """
        Function that creates the menu
        """

        menubar = tk.Menu(self)

        # creation of the file menu
        self.file_menu = tk.Menu(menubar, tearoff=0)
        self.file_menu.add_command(label="Open", command=self.on_open)
        self.file_menu.add_command(label="Save", command=self.on_save)
        self.file_menu.add_command(label="Save As", command=self.on_save_as)
        self.file_menu.add_command(label="Quit", command=self.on_quit)
        menubar.add_cascade(label="File", menu=self.file_menu)

        # creation of the section menu
        self.section_menu = tk.Menu(menubar, tearoff=0)
        self.section_menu.add_command(label='Load section', command=self.on_open_sections)
        menubar.add_cascade(label='Section', menu=self.section_menu)

        # deactivation of the 'Save' and 'Save as' menu items
        self.file_menu.entryconfig(1, state='disabled')
        self.file_menu.entryconfig(2, state='disabled')

        # deactivation of the 'Load Section' menu items
        self.section_menu.entryconfig(0, state='disabled')

        self.config(menu=menubar)

    def bind_functions(self):

        """
        Function that binds the widgets to their respective functions
        """

        # functions associated with a change in the domain coordinates
        self.entry_imin.bind("<Return>", self.on_change_coordinates)
        self.entry_imax.bind("<Return>", self.on_change_coordinates)
        self.entry_jmin.bind("<Return>", self.on_change_coordinates)
        self.entry_jmax.bind("<Return>", self.on_change_coordinates)

        # functions associated with a change in the model
        self.combobox_model.bind("<<ComboboxSelected>>", self.on_change_model)

        # functions associated with a change in the segment
        self.combobox_segment.bind("<<ComboboxSelected>>", self.on_change_segment)

    def on_quit(self):

        """
        Function that is called when the Quit menu is activated.

        It destroys the Tkinter frame
        """

        self.destroy()

    def on_open(self):

        """
        Function that is called when the Open menu is activated.
        Allows to open a .nc meshfile
        """

        options = {}
        options['defaultextension'] = '.nc'
        options['filetypes'] = [('Meshgrid files', '.nc')]
        options['title'] = u'Choose the .nc file containing meshgrid to be opened'

        try:
            filename = tkFileDialog.askopenfilename(**options)
        except:
            print('perdu')
            print "Unexpected error:", sys.exc_info()[0]
            raise

        if filename:

            # we switch the working mode to 'edit domain'
            self.radiobox_editmode_lab.set(0)

            # we define the gridfilename attribute, in which
            # scale factors will be read
            self.gridfilename = filename

            # we extracts the lon/lat/mask/bathy arrays in the file
            # we draw the mask and the domain
            self.bmapplot.init_domain()

            # we activate the widgets associated with the
            # domain edition
            self.entry_imin.configure(state='normal')
            self.entry_imax.configure(state='normal')
            self.entry_jmin.configure(state='normal')
            self.entry_jmax.configure(state='normal')

            # we activate the 'edit domain' widget
            self.radiobox_editmode0.configure(state='normal')
            self.radiobox_editmode1.configure(state='disabled')
            self.radiobox_editmode2.configure(state='disabled')

            # we activate the 'load section' menu
            self.file_menu.entryconfig(1, state='disabled')
            self.file_menu.entryconfig(2, state='disabled')
            self.section_menu.entryconfig(0, state='normal')

    def on_save(self):

        """
        Function that is called when the Save menu is activated.

        If a filename is defined, it saves into the file.
        Else, it prompts a file window.
        """

        # if the section_file_name is None -> SaveAs
        if self.outputfilename is None:
            options = {}
            options['defaultextension'] = '.pygo'
            options['filetypes'] = [('PyPago Grid Files', '.pygo')]
            filename = tkFileDialog.asksaveasfilename(**options)
            if filename:
                self.outputfilename = filename
                pypago.pyio.save(self.bmapplot.gridmodel, self.outputfilename.replace(".pygo", "_grid.pygo"))
                pypago.pyio.save(self.bmapplot.model_sections, self.outputfilename.replace(".pygo", "_gridsec.pygo"))
        else:
            # if section file exists -> we save
            pypago.pyio.save(self.bmapplot.gridmodel, self.outputfilename.replace(".pygo", "_grid.pygo"))
            pypago.pyio.save(self.bmapplot.model_sections, self.outputfilename.replace(".pygo", "_gridsec.pygo"))

    def on_save_as(self):

        """
        Function that is called when the Save As menu is activated.

        It prompts a file window, in which the sections will be saved.
        """

        options = {}
        options['defaultextension'] = '.pygo'
        options['filetypes'] = [('PyPago Grid Files', '.pygo')]
        filename = tkFileDialog.asksaveasfilename(**options)
        if filename:
            self.outputfilename = filename
            pypago.pyio.save(self.bmapplot.gridmodel, self.outputfilename.replace(".pygo", "_grid.pygo"))
            pypago.pyio.save(self.bmapplot.model_sections, self.outputfilename.replace(".pygo", "_gridsec.pygo"))

    def on_open_sections(self):

        """
        Function that is called when the Load Section menu is activated.
        It extracts the scale factors on the domain
        defined by `imin`, `imax`, `jmin`, `jmax` variables. And it generates
        the section staircases
        """

        options = {}
        options['defaultextension'] = '.pygo'
        options['filetypes'] = [('PyPago Grid Files', '.pygo')]
        options['initialdir'] = '{0}/.pago'.format(os.environ['HOME'])
        options['title'] = u'Choose the .pygo file containing sections to be opened'
        filename = tkFileDialog.askopenfilename(**options)

        if filename:

            # we load the sections and put them in the bmapplot.sections attributes
            self.bmapplot.sections = pypago.pyio.load(filename)

            # we switch the working mode to check sections
            self.radiobox_editmode_lab.set(2)

            # we disable the widgets associated with domain edition
            self.entry_imin.configure(state='disabled')
            self.entry_imax.configure(state='disabled')
            self.entry_jmin.configure(state='disabled')
            self.entry_jmax.configure(state='disabled')

            # we extract the scale factors and compute surface/volume
            # on the subdomain
            self.bmapplot.finalize_domain()

            # we generate section staircases
            self.bmapplot.generate_sections()

            # and we activate the 'edit section' and 'check section' widgets
            self.radiobox_editmode1.configure(state='normal')
            self.radiobox_editmode2.configure(state='normal')

            # we activate the 'load section' menu
            self.file_menu.entryconfig(1, state='normal')
            self.file_menu.entryconfig(2, state='normal')

    def on_change_model(self, event):   # pylint: disable=unused-argument

        """
        Function that is called when the model is changed
        """

        # we simply change the modelproc attribute
        self.modelproc = self.combobox_model.get()

    def on_change_segment(self, event):   # pylint: disable=unused-argument

        """
        Function that is called when the segment is changed
        """

        # we change the segment
        self.bmapplot.model_sections[self.bmapplot.indsec].dire[self.bmapplot.segment] = self.combobox_segment.get()

    def on_del_button(self):

        """
        Function which is called when the delete button is activated
        """

        # We reset the description of the i and j coordinates
        self.label_i_lab.set('i coordinates:')
        self.label_j_lab.set('j coordinates:')

        # we disable the delete button
        self.button_del.configure(state='disabled')

        # we reset the section name status
        self.label_secname_lab.set('Section name:')

        # we remove the indsec'th sections and model_sections
        self.bmapplot.sections.pop(self.bmapplot.indsec)
        self.bmapplot.model_sections.pop(self.bmapplot.indsec)
        self.bmapplot.secname = None  # we set the secname to None

        # we remove the scatter plot
        self.bmapplot.scatter.remove()

        # we redraw the sections
        self.bmapplot.plot_sections()

    def on_change_status(self):

        """
        Function which is called when we change the edition mode
        """

        # we first disconnect the dompic and editsec events (if possible)
        # prevents from conflicting events
        self.bmapplot.disconnect_events()

        # If we chose to edit the domain
        if self.radiobox_editmode_lab.get() == 0:

            # We deactivate the save and save as item
            self.file_menu.entryconfig(1, state='disabled')
            self.file_menu.entryconfig(2, state='disabled')

            self.reprocess_grid = True

            # We deselect the section
            self.bmapplot.secname = None

            # Here we draw the domain
            self.bmapplot.plot_map_global()
            self.bmapplot.plot_dom()
            self.bmapplot.dompic_event = self.bmapplot.figure.canvas.mpl_connect('pick_event',
                                                                                 self.bmapplot.on_dompic)

            # we enable the widgets for the domain edition
            self.entry_imin.configure(state='normal')
            self.entry_imax.configure(state='normal')
            self.entry_jmin.configure(state='normal')
            self.entry_jmax.configure(state='normal')

            # we reinitialize the status (i,j coordinates and name) of the section
            self.label_i_lab.set('i coordinates:')
            self.label_j_lab.set('j coordinates:')
            self.label_secname_lab.set('Section name:')

            # We disable all the widgets for the section editions
            self.combobox_segment.configure(state='disabled')
            self.button_del.configure(state='disabled')

        elif self.radiobox_editmode_lab.get() == 1:  # If we are willing to edit the sections

            # We deactivate the save and save as menu items
            self.file_menu.entryconfig(1, state='disabled')
            self.file_menu.entryconfig(2, state='disabled')

            # We deactivate the widgets for the domain edition
            self.entry_imin.configure(state='disabled')
            self.entry_imax.configure(state='disabled')
            self.entry_jmin.configure(state='disabled')
            self.entry_jmax.configure(state='disabled')

            # we plot the sections in 'edit mode'
            self.bmapplot.plot_sections()

        else:  # if we are verifying the section

            # we enable the save and save as menus
            self.file_menu.entryconfig(1, state='normal')
            self.file_menu.entryconfig(2, state='normal')
            self.bmapplot.secname = None  # we reset the selected name

            # we finalize the domain only if the domain has been edited!
            if self.reprocess_grid:
                self.bmapplot.finalize_domain()

            # we regenerate all the sections
            self.bmapplot.generate_sections()

            # We reinitialise the section information widgets
            self.label_i_lab.set('i coordinates:')
            self.label_j_lab.set('j coordinates:')
            self.label_secname_lab.set('Section name:')

            # We deactivate all the edition widgets
            self.entry_imin.configure(state='disabled')
            self.entry_imax.configure(state='disabled')
            self.entry_jmin.configure(state='disabled')
            self.entry_jmax.configure(state='disabled')
            self.combobox_segment.configure(state='disabled')
            self.button_del.configure(state='disabled')

    def on_change_coordinates(self, event):    # pylint: disable=unused-argument

        """
        Function that handles the change coordinates
        """

        # if the new coordinates can be converted to integer
        try:

            imin = np.int(self.entry_imin_lab.get())
            imax = np.int(self.entry_imax_lab.get())
            jmin = np.min([np.int(self.entry_jmin_lab.get()), np.int(self.entry_jmax_lab.get())])
            jmax = np.max([np.int(self.entry_jmin_lab.get()), np.int(self.entry_jmax_lab.get())])

            # if they are within the domain limits, the attributes are changed
            if (imin >= self.bmapplot.lonlim[0]) & \
               (imin <= self.bmapplot.lonlim[-1]) & \
               (imax >= self.bmapplot.lonlim[0]) & \
               (imax <= self.bmapplot.lonlim[-1]) & \
               (jmax >= self.bmapplot.latlim[0]) & \
               (jmax <= self.bmapplot.latlim[-1]) & \
               (jmin >= self.bmapplot.latlim[0]) & \
               (jmin <= self.bmapplot.latlim[-1]):

                self.bmapplot.imin = imin
                self.bmapplot.imax = imax
                self.bmapplot.jmin = jmin
                self.bmapplot.jmax = jmax

            # if they are out of the domain limits, the attributes are resetted to
            # their initial values
            else:
                self.bmapplot.imin = self.bmapplot.lonlim[0]
                self.bmapplot.imax = self.bmapplot.lonlim[-1]
                self.bmapplot.jmin = self.bmapplot.latlim[0]
                self.bmapplot.jmax = self.bmapplot.latlim[-1]

            self.bmapplot.plot_dom()  # we redraw the domain

        except:
            pass

        # We modify the display with the imin, imax, jmin and jmax
        self.entry_imin_lab.set(str(self.bmapplot.imin))
        self.entry_imax_lab.set(str(self.bmapplot.imax))
        self.entry_jmin_lab.set(str(self.bmapplot.jmin))
        self.entry_jmax_lab.set(str(self.bmapplot.jmax))


if __name__ == "__main__":
    PAGOFRAME = GridModel(None)
    PAGOFRAME.title("PyPago Grid Model")
    PAGOFRAME.mainloop()
    sys.exit(0)
