#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Executable function that allows to
extract spatial domains, by giving its section
boundaries.

The use of this function is fully described in
:ref:`definition_of_closed_domain`

First argument is the Pygo grid file, the second argument is the
Pygo gridded sections file, the third argument
is the name of the output file.

.. code-block:: bash
    make_areas.py grid.pygo gridsec.pygo domain.pygo

.. warning::
    The directory in which this script is executed must
    contain a *param.py* file.

"""

# pylint: disable=superfluous-parens

import sys
import pylab as plt
import numpy as np
import scipy.misc as misc

from pypago.disp import PypagoErrors
import pypago.plot
import pypago.pyio
import pypago.misc
import pypago.areas

def _process_sections(grid, gridsec, secnames, wheretogo):

    """
    Function that initialises the mask array that will
    be drawn in the .png file, depending on the sections
    that define the domain and their orientations.

    :param dict data: The dictionary that contains the `grid`
       and `gridsec` variables.
    :param pago_obj areaint: the area object which has been initialised
    :return: a numpy.array that contains the initialised mask
    :rtype: numpy.array
    """

    mask = grid.mask.copy()
    signes = []

    for isec in xrange(0, len(secnames)):
        sec = pypago.misc.findsecnum(gridsec, secnames[isec])
        if wheretogo[isec] == 'in':
            signe = 1
        else:
            signe = -1
        signes.append(signe)

        veci = gridsec[sec].veci
        vecj = gridsec[sec].vecj
        orientation = gridsec[sec].orient*signe
        faces = gridsec[sec].faces

        for l in xrange(0, len(veci)):
            if faces[l] == 'N':
                if orientation[l] == 1:
                    mask[vecj[l]+1, veci[l]] = 2
                else:
                    mask[vecj[l], veci[l]] = 2
            else:
                if orientation[l] == 1:
                    mask[vecj[l], veci[l]] = 2
                else:
                    mask[vecj[l], veci[l]-1] = 2

    return mask, signes


def _definearea(grid, gridsec):

    """
    Function which is called any time a new area is defined.
    It initialises the area as a |pypago| object and it calls the
    then asks the user to define the domain by giving the
    sections and whether the sections are oriented toward or
    outward of the domain. It creates a `.png` file that
    must be edited by an external software (Gimp for instance).
    The grid points within the domain must be set to white.

    :param dict data: The dictionary that contains the `grid`
       and `gridsec` variables.
    :return: a |pypago| object that contains all the variables relative to
       the domain
    :rtype: pypago.pyio.pago_obj
    """

    # initialisation of the output as a pago_obj
    # areaint = pypago.pyio.pago_obj({})
    secgrid_list = [s.name for s in gridsec]

    # Small algorithm that checks that the user
    # input are OK.
    passed = False
    while not passed:
        try:
            areaname = raw_input('name of the area? ')
            if len(areaname) == 0:
                message = "You should provide a section name.\n"
                message += "Please try again"
                print(message)
            else:
                passed = True
        except:
            pass
    
    passed = False
    while not passed:
        try:
            secname = raw_input('give names of the sections ' + \
                                '(separate the names by a space) ')
            secname = secname.split()
            if len(secname) > 0:
                test = np.array([s in secgrid_list for s in secname])
                if np.any(test==0):
                   message = "One ore more user defined sections do not belong "
                   message += "to the list of grid section. Please try again"
                   print(message)
                else:
                    passed = True
            else:
                message = "The length of the section names is 0."
                message += "Please retry."
                print(message)
        except:
            pass

    passed = False
    while not passed:
        try:
            wheretogo = raw_input('give orientation of the sections ' + \
                                  '(in: directed toward the basin/out directed ' + \
                                  'out of the area) ')

            wheretogo = wheretogo.split()
            if len(wheretogo) == len(secname):
                test = np.array([s in ["in", "out"] for s in wheretogo])
                if np.any(test==0):
                   message = "The elements of the wheretogo list must be 'in' or 'out'. "
                   message += "Please try again"
                   print(message)
                else:
                    passed = True
            else:
                message = "The length of the wheretogo list must be the same "
                message += "as the length of the secname list. Try again."
                print(message)
        except:
            pass

    # we process the sections that close the domain
    mask, signes = _process_sections(grid, gridsec, secname, wheretogo)

    misc.imsave('mask_init_' + areaname + '.png', mask)

    domain_ok = False

    while not domain_ok:

        print('Edit the mask_init_' + areaname + '.png file using gimp, ' + \
              'paint or any other software. Fill in white the area enclosed ' + \
              'in your boundaries. Save the new png file as ' + \
              'mask_init_' + areaname + '_bis.png.')

        passed = False
        while not passed:
            raw_input('When done, press any key ')
            try:
                mask = misc.imread('./mask_init_' + areaname + '_bis.png')
                passed = True
            except:
                message = "The ./mask_init_" + areaname + '_bis.png file\n' 
                message += "has not been found.\n"
                print(message)
                pass

        mask[np.nonzero(mask != 255)] = 0
        mask[np.nonzero(mask == 255)] = 1

        # to remove white stripes on land
        mask = mask*grid.mask
    
        plt.figure(2)
        plt.clf()
        pypago.plot.plot_dom_mask(grid, gridsec, mask)

        print('Check that the domain mask is well defined ' + \
              '(i.e. between the section lines). If not, reconsider ' + \
              'the png edition')

        passed = False
        while not passed:
            try:
                domain_ok = int(raw_input('Is it ok? (0 if not) '))
                passed = True
            except:
                pass
        
    i, j = np.nonzero(mask == 1)
    areaint = pypago.areas.Areas(grid, areaname, i, j, secnames=secname, signs=signes)


    return areaint


def extract_areas_fromsec(grid, gridsec):

    """
    Main function for the domain extraction.

    :param pypago.grid.Grid grid: Grid object
    :param list gridsec: List
     of :py:class:`pypago.sections.GridSection` objects

    :return: A list of :py:class:`pypago.areas.Areas` object
    :rtype: list

    """

    plt.ion()

    areas = []

    # initialisation of the plot:
    #
    # - shows the wet and land points
    # - shows the sections, their orientations and names
    #_initplot(grid, gridsec)
    pypago.plot.plot_dom_mask(grid, gridsec, mask=None, ax=None)

    defnewarea = True  # true if we must define sections

    while defnewarea:
        # if we are to define a new area
        # we call the _definearea function
        # and append the output to the final list
        newarea = _definearea(grid, gridsec)
        areas.append(newarea)
        passed = False
        while not passed:
            try:
                defnewarea = int(raw_input('Define another area? (0 or 1) '))
                passed = True
            except:
                pass

    plt.ioff()

    return areas


if __name__ == '__main__':

    argv = sys.argv[1:]
    if len(argv) != 3:
        message = 'The number of inputs arguments must be equal to 3.\n'
        message += ' - The Pygo grid file\n'
        message += ' - The Pygo grids section file\n'
        message += ' - The Pygo output file\n'
        message += '%d arguments were provided.\n' %len(argv)
        message += 'This program will stop'
        raise PypagoErrors(message)

    modgrid = pypago.pyio.load(argv[0])
    modgridsec = pypago.pyio.load(argv[1])
    modareas = extract_areas_fromsec(modgrid, modgridsec)
    pypago.pyio.save(modareas, argv[2])
