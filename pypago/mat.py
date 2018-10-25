# -*- coding: utf-8 -*-

"""
Module that handles function to convert :file:`.mat` files into \
:file:`.pygo` files, and conversely.
"""

import scipy.io
import pypago.pyio
import pypago.sections
import pypago.grid

def grid_topygo(filename, modelname, nlon):
    """
    Converts a :file:`.mat`
    containing grid informations into
    a :file:`.pygo` file

    :param str filename: name of the :file:`.mat` file

    """

    gridout = pypago.grid.Grid(None, None, None,
                               None, None)

    fin = scipy.io.loadmat(filename)
    model_grid = fin['MODEL_grid'][0]
    gridout.jmin = model_grid['lats'][0][0][0] - 1
    gridout.jmax = model_grid['latn'][0][0][0] - 1
    gridout.imin = model_grid['lonw'][0][0][0] - 1
    gridout.imax = model_grid['lone'][0][0][0] - 1
    gridout.mask = model_grid['mask'][0]
    gridout.bathy = model_grid['bathy'][0]
    gridout.lont = model_grid['lont'][0]
    gridout.latt = model_grid['latt'][0]
    gridout.filename = filename
    gridout.modelname = modelname
    gridout.nlon = nlon

    output = {}
    output['MODEL_grid'] = fin['MODEL_grid'][0]
    pypago.pyio.save(output, filename.replace('.mat', '.pygo'))

    return gridout

def gridsec_topygo(filename, modelname, nlon):

    """
    Converts a :file:`.mat`
    containing grid section endpoints into
    a :file:`.pygo` file

    :param str filename: name of the :file:`.mat` file

    """

    output = []

    fin = scipy.io.loadmat(filename)

    model_time = fin['MODEL_time'][0]
    model_sec = fin['MODEL_sections']
    model_grid = fin['MODEL_grid'][0]
    jmin = model_grid['lats'][0][0][0] - 1
    jmax = model_grid['latn'][0][0][0] - 1
    imin = model_grid['lonw'][0][0][0] - 1
    imax = model_grid['lone'][0][0][0] - 1

    for mod in model_sec[:]:
        secint = pypago.sections.GridSection(None, None)
        secint.name = mod['name'][0][0]
        secint.faces = mod['faces'][0]
        secint.orient = mod['orient'][0]
        secint.depthvect = mod['depthvect'][0]
        secint.areavect = mod['areavect'][0]
        secint.lvect = mod['lvect'][0]
        secint.lengthvect = mod['lengthvect'][0]
        secint.veci = mod['veci'][0] - 1
        secint.vecj = mod['vecj'][0] - 1
        secint.vecs = mod['vecs'][0]
        secint.vect = mod['vect'][0]
        secint.vecv = mod['vecv'][0]
        secint.jmin = jmin
        secint.jmax = jmax
        secint.imin = imin
        secint.imax = imax
        secint.modelname = modelname
        secint.nlon = nlon

        output.append(secint)

    outputdict = {}
    outputdict['MODEL_sections'] = output
    outputdict['MODEL_time'] = model_time

    pypago.pyio.save(outputdict, filename.replace('.mat', '.pygo'))

    return output

def secmat_topygo(filename):

    """
    Converts a :file:`.mat`
    containing section endpoints into
    a :file:`.pygo` file

    :param str filename: name of the :file:`.mat` file

    """

    fin = scipy.io.loadmat(filename)
    sectionsmat = fin['sections'][0]

    secpygo = []

    for secmat in sectionsmat:

        secint = pypago.sections.Section(secmat['name'][0],
                                         secmat['lon'][0],
                                         secmat['lat'][0],
                                         secmat['dir'])

        secpygo.append(secint)

    pypago.pyio.save(secpygo, filename.replace('.mat', '.pygo'))

    return secpygo


if __name__ == '__main__':

    import pylab as plt

    filename_test = 'natl_ORCA025.L75-GRD88_y2000m07_sections.mat'
    #output = gridsec_topygo(filename, 'NEMO', 100)

    gridout_test = pypago.grid.Grid(None, None, None,
                                    None, None)

    fin_test = scipy.io.loadmat(filename_test)
    model_grid_test = fin_test['MODEL_grid'][0]
    gridout_test.jmin = model_grid_test['lats'][0][0][0] - 1
    gridout_test.jmax = model_grid_test['latn'][0][0][0] - 1
    gridout_test.imin = model_grid_test['lonw'][0][0][0] - 1
    gridout_test.imax = model_grid_test['lone'][0][0][0] - 1
    gridout_test.mask = model_grid_test['mask'][0]
    gridout_test.bathy = model_grid_test['bathy'][0]
    gridout_test.lont = model_grid_test['lont'][0]
    gridout_test.latt = model_grid_test['latt'][0]

    gridout_test.plot_mask()
    plt.show()
