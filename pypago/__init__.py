# -*- coding: utf-8 -*-
"""
PAGO License TDB

.. todo:

   description and others : get them from setup.py
   or like in setup.py
   same in conf.py

"""

from __future__ import print_function

import os

import pkg_resources  # part of setuptools
try:
    __version__ = pkg_resources.require("pypago")[0].version
except:
    VERSION_FILE = os.path.join('{0}/../'.format(os.path.dirname(__file__)), 'VERSION')
    with open(VERSION_FILE, 'r') as infile:
        __version__ = infile.read().strip()

__description__ = "Python Physical Analysis of Gridded Ocean"
__author_email__ = "pago-dev@groupes.renater.fr"
