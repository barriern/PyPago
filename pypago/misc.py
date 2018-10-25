# -*- coding: utf-8 -*-

"""
Module that contains miscellaneous functions (functions to determine the
indexes of a section or area within a list) and the class that allows to define
|matlab|-like structures
"""

from __future__ import print_function
import inspect
import numpy as np
from pypago.disp import PypagoErrors


def extract_mask(array):

    """ True where missing data (masked or NaN) """

    output = (np.isnan(array) | np.ma.getmaskarray(array))
    return output


def make_percentile_cmap(array, perc):

    """
    Returns colormaps limits using percentiles. The
    minimum color value.

    :param array numpy.array: Input array

    :param float perc : Percentile (in percentage)

    :returns: A tuple containing the lower and upper colorbar limits
     associated with the array percentiles.

    :rtype: tuple

    """

    # converts from ND to 1D
    array = np.ravel(array)

    # mask data where missing and extracts only unmasked
    # data
    array = np.ma.masked_where(array != array, array)
    array = array[np.logical_not(np.ma.getmaskarray(array))]

    if array.size == 0:
        cmin, cmax = None, None

    else:
        # computes the percentiles
        cmin = np.percentile(array, perc)
        cmax = np.percentile(array, 100 - perc)

        # rounds the percentiles (unnecessary???)
        # cmin = unit_round(cmin, "floor")
        # cmax = unit_round(cmax, "ceil")

    return cmin, cmax


def extract_str(attrnames, dataclass):

    """
    Extract a string which will be used in the __str__ method.
    The string list the attributes and display either their shape
    (for numpy arrays) or their values

    :param list attrnames: List of attributes to displau.
    :param object dataclass: Input object
    :return: A string
    :rtype: str

    """

    output = ''
    for attr in attrnames:
        val = getattr(dataclass, attr)
        if isinstance(val, np.ndarray):
            output += '  -%s: %s\n' % (attr, str(val.shape))
        else:
            output += '  -%s: %s\n' % (attr, str(val))

    return output


def extract_attrlist(dataclass):

    """
    Extract the names of the oject attributes which are not instance
    method and which do not start by a _ character

    :param object dataclass: Input object
    :return: A list of strings
    :rtype: list
    """

    attrnames = ([attr for attr in dir(dataclass) if attr[0] != '_'])
    attrnames = [attr for attr in attrnames if not inspect.ismethod(getattr(dataclass, attr))]
    return attrnames


def finddomnum(model_areas, areaplotname):

    """
    Determines the index of a domain in a list of domains

    :param list model_areas: list of the Pago areas where to find the area
    :param str areaplotname: name of the area
    :return: index of the area within the list
    :rtype: int
    """

    for inddom in xrange(0, len(model_areas)):
        if model_areas[inddom].areaname == areaplotname:
            return inddom

    message = r'Domain {0} does not exist'.format(areaplotname)
    error = PypagoErrors(message)
    raise error


def findsecnum(model_sections, secplotname):

    """
    Determines the index of a section in a list of sections

    :param list model_sections: list of the Pago sections where to find the
       section
    :param str secplotname: name of the section
    :return: index of the section within the list
    :rtype: int

    :raise error: if `secplotname` not found in `MODEL_sections`
    """

    for indsec in xrange(0, len(model_sections)):
        if model_sections[indsec].name == secplotname:
            return indsec

    message = r'Section {0} does not exist'.format(secplotname)
    error = PypagoErrors(message)
    raise error
