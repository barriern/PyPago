#! /usr/bin/env python

"""
Executable file for domain extraction.

First argument is the name of the Cordinates pygo file,
the second argument is the name of the output file.

.. code-block:: bash
    make_grid.py NEMO mesh_file.nc grid.pygo

.. warning::
    The directory in which this script is executed must
    contain a *param.py* file.

"""

import sys
import pylab as plt
from pypago.disp import PypagoErrors
import pypago.coords
import pypago.grid
import pypago.pyio
import pypago.plot


def extract_grid(coords):
    
    # check that the file has the proper type (Coords object)
    if not isinstance(coords, pypago.coords.Coords):
        message = 'The file content is not a "Coords" object.\n'
        message += 'This program will be stopped.'
        raise PypagoErrors(message)

    # extract the lon/lat
    nlat, nlon = coords.mask.shape
    agree = False    # True when the user agrees with the domain

    # intialises the figure where with the plots
    plt.ion()
    fig = plt.figure(1)
    ax = plt.gca()
    coords.plot_mask()
    plt.show()

    # if agree is false
    while not(agree):

        message = '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n'
        message += 'Enter the index of the southern and northern limits:\n'
        message += ' 1 <= jmin < jmax <= %d\n' %(nlat-2)
        message += ' 0 <= imin, imax <= %d\n'  %(nlon-1)
        print(message)

        # ask for imin, imax, lons and jmax values
        try:
            passed = False
            while not passed:
                try:
                    imin = int(raw_input('imin = '))
                    passed = True
                except:
                    pass
            
            passed = False
            while not passed:
                try:
                    imax = int(raw_input('imax = '))
                    passed = True
                except:
                    pass
            
            passed = False
            while not passed:
                try:
                    jmin = int(raw_input('jmin = '))
                    passed = True
                except:
                    pass
            
            passed = False
            while not passed:
                try:
                    jmax = int(raw_input('jmax = '))
                    passed = True
                except:
                    pass

            if jmax < jmin:
                print("jmin and jmax have been swithed!")
                jmax, jmin = jmin, jmax

            if imax == imin:
                print("imax and imin are equal: they have been set to None")
                imax = None
                imin = None

            # extracts the grid
            grid = pypago.grid.create_grid(coords, jmin, jmax,
                                           imin, imax)

            # we clear the domain lines and redraws it with
            # the new values
            plt.figure(1)
            ax.lines = []
            grid.plot_dom(ax=ax)

            # drawing of the domain mask
            plt.figure(2)
            ax2 = plt.gca()
            plt.cla()
            grid.plot_mask()

            # ask if the user is ok with the domain
            agree = int(raw_input("Is the domain ok? (0 if not) "))

        except:
            pass

    return grid


if __name__ == "__main__":

    argv = sys.argv[1:]

    if len(argv) != 3:
        message = 'The number of input arguments must be 2.'
        message += 'Currently, %d arguments provided.' %len(argv)
        message += 'This program will be stopped.\n'
        message += "Syntax:\n"
        message += "> make_grid.py NEMO mesh_file.nc grid.pygo"
        raise PypagoErrors(message)

    modelname = argv[0]
    filename = argv[1]
    outfile = argv[2]

    # load the coordinate file
    #coords = pypago.pyio.load(filename)
    coords = pypago.coords.create_coord(modelname, filename)

    grid = extract_grid(coords)

    pypago.pyio.save(grid, outfile)
