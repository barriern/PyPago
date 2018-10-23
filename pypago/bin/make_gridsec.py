#! /usr/bin/env python

"""
Executable function that allows to load section endpoints
and to eventually edit them interactively (bad orientation,
out of land points, etc).

First argument is the grid file, second argument is the
file containing the list of section endpoints, third argument
is the name of the output file.

.. code-block:: bash
    make_sections.py grid.pygo endpoints.pygo gridsec.pygo

.. warning::
    The directory in which this script is executed must
    contain a *param.py* file.

"""

# pylint: disable=superfluous-parens

import sys
import pylab as plt
from pypago.disp import PypagoErrors
import pypago.coords
import pypago.grid
import pypago.pyio
import pypago.plot
import pypago.sections
import numpy as np

def extract_gridsec(grid, section):
    
    # extract the list of sections to analyse
    # gridsec = sections on the domain
    # badsec = indexes of out of dom. sections
    gridsec, badsec = pypago.sections.extract_grid_sections(grid, section)

    # if there is out of domain sections, print a warning
    if len(badsec) > 0:
        message = "@@@@@@@@@@@@@ Warning @@@@@@@@@@@@@\n"
        seclist = ['-%s\n' %section[i].name for i in badsec]
        message += ''.join(seclist)
        message += 'sections are out of the domain. '
        message += 'These sections have been discarded.\n'
        message += "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
        print(message)
        passed = False
        while not passed:
            try:
                keepgoing = int(raw_input('Do you want to continue? (0 if not) '))
                passed = True
            except:
                pass
        if not(keepgoing):
            sys.exit(0)

    # if there is no good section, we stop the program
    if len(gridsec) == 0:
        message = 'No sections could have been defined. '
        message += 'This programm will be stopped. '
        raise PypagoErrors(message)

    # extraction of the good section endpoints
    section = [section[i] for i in range(0, len(section)) if i not in badsec]

    # draw the section staircases
    # for section checking
    plt.ion()
    fig = plt.figure(1)
    ax1 = plt.gca()
    grid.plot_mask()
    for sec in gridsec:
        sec.plotsecfaces()
    plt.show()

    # extracts the list of section names
    section_list = [sec.name for sec in gridsec]

    # agree = True if the sections are well defined.
    passed = False
    while not passed:
        try:
            agree = int(raw_input('Are the sections well defined? (0 if not) '))
            passed = True
        except:
            pass

    # while sections are not well defined.
    while not agree:

        # asking for the section name to correct.
        # ask until the name is within the section_list
        section_name = '-9999999'
        while section_name not in section_list:
            section_name = raw_input('Which section is wrong? ')

        # find the index of the section
        indsec = pypago.misc.findsecnum(gridsec, section_name)

        # secint is the associated gridsec object
        secint = gridsec[indsec]

        # draws the sections on the domain (i, j) values
        fig2 = plt.figure(2)
        plt.cla()
        ax2 = plt.gca()
        grid.plot_mask()
        plt.plot(secint.i, secint.j, marker='o')
        plt.show()

        # print some informations about the section
        print('Section %s:' %secint.name)
        print(' -i: %s' %str(secint.i))
        print(' -j: %s' %str(secint.j))
        print(' -dire: %s' %str(secint.dire))

        # copy the i,j,dire attributes
        newi = secint.i.copy()
        newj = secint.j.copy()
        newdir = secint.dire.copy()

        passed = False   # passed=True if the definition of i worked
        while not passed:
            try:
                # ask whether to change i. if this statement
                # fails, then we loop again
                if int(raw_input('Change i? ')):
                    newi = raw_input(' -i: ').split()
                    newi = [int(n) for n in newi]
                testi = np.array([(temp >= 0) & (temp <= nlon-1) for temp in newi])
                if np.any(testi==False):
                    message = "The values of indi must be within the domain boundary.\n"
                    message += "Please try again."
                    print(message)
                else:
                    passed = True
            except:
                pass
        

        passed = False
        while not passed:
            try:
                if int(raw_input('Change j? ')):
                    newj = raw_input(' -j: ').split()
                    newj = [int(n) for n in newj]
                testj = np.array([(temp >= 0) & (temp <= nlat-1) for temp in newj])
                if np.any(testj==False):
                    message = "The values of indj must be within the domain boundary.\n"
                    message += "Please try again."
                    print(message)
                else:
                    passed = True
            except:
                pass

        passed = False
        while not passed:
            try:
                if int(raw_input('Change dir? ')):
                    newdir = raw_input(' -dire: ').split()
                testdir = np.array([temp in ["NW", "NE", "SE", "SW"] for temp in newdir])
                if np.any(testdir==False):
                    message = "The directions must be NW, NE, SE, SW.\n"
                    message += "Please try again."
                    print(message)
                else:
                    passed = True
            except:
                pass

        # we check that
        if (len(newdir) != len(newi)-1) | (len(newi) != len(newj)):
            message = "The new i and new j lists must have the same length.\n"
            message += "The new dir must be of length len(newi) - 1.\n"
            message += "Currently, len(newi)=%d, len(newj)=%d, len(newdir)=%d" %(len(newi), len(newj), len(newdir))
            message += "Please try again."
            print(message)
            continue

        # modification of the section endpoints and direction
        section[indsec].lon = np.array(grid.lont[newj, newi])
        section[indsec].lat = np.array(grid.latt[newj, newi])
        section[indsec].dire = np.array(newdir)

        # here, we recompute only the modified section (without recalculating everything)
        gridsectmp, badsectmpt = pypago.sections.extract_grid_sections(grid, section[indsec])
        gridsec[indsec] = gridsectmp[0]
        secint = gridsec[indsec]

        # we redraw the two plots
        plt.figure(1)
        plt.cla()
        grid.plot_mask()
        for sec in gridsec:
            sec.plotsecfaces()

        fig2 = plt.figure(2)
        plt.cla()
        ax2 = plt.gca()
        grid.plot_mask()
        plt.plot(secint.i, secint.j, marker='o')
        plt.show()

        agree = int(raw_input('Are the sections now well defined? (0 if not) '))

    return gridsec


if __name__ == "__main__":

    argv = sys.argv[1:]

    if len(argv) != 3:
        message = 'The number of input arguments must be at least 3.'
        message += 'Currently, %d arguments provided.' %len(argv)
        message += 'This program will be stopped.\n'
        message += "Syntax:\n"
        message += "> make_gridsec.py grid.pygo sec.pygo gridsec.pygo"
        raise PypagoErrors(message)

    gridfile = argv[0]
    sectionfile = argv[1]
    savefile = argv[2]

    grid = pypago.pyio.load(gridfile)
    nlat, nlon = grid.mask.shape

    # check that the file has the proper type (Coords object)
    if not isinstance(grid, pypago.grid.Grid):
        message = 'The file content is not a "Grid" object.\n'
        message += 'This program will be stopped.'
        raise PypagoErrors(message)

    section = pypago.pyio.load(sectionfile)

    gridsec = extract_gridsec(grid, section)

    # when we are done, saving of the gridfile
    pypago.pyio.save(gridsec, savefile)
