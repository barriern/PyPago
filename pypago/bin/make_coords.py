#! /usr/bin/env python

"""
Executable file that extracts
the coordinates of a specific model
mesh file.

Fist argument is the model name, second argument
is the file location, third argument is the name
of the output file.

.. code-block:: bash
    make_coords.py NEMO data/ORCA1_mesh_mask.nc coords.pygo

.. warning::
    The directory in which this script is executed must
    contain a *param.py* file.
"""

import sys
from pypago.disp import PypagoErrors
import pypago.coords
import pypago.grid
import pypago.pyio
import pypago.plot

if __name__ == "__main__":

    argv = sys.argv[1:]

    if len(argv) < 3:
        message = 'The number of input arguments must be at least 2.'
        message += 'Currently, %d arguments provided.' %len(argv)
        message += 'This program will be stopped.\n'
        message += "Syntax:\n"
        message += "> make_coords.py NEMO data/ORCA1_mesh_mask.nc coords.pygo"

        raise PypagoErrors(message)

    modelname = argv[0]
    filename = argv[1]
    outfilename = argv[2]

    coords = pypago.coords.create_coord(modelname, filename)

    pypago.pyio.save(coords, outfilename)
