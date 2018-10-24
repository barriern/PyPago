#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Setup for packaging |pago|"
"""

__docformat__ = "restructuredtext en"

import os
from setuptools import setup, find_packages

VERSION_FILE = 'VERSION'
version = open(VERSION_FILE).read().strip()

setup(
    name="pypago",
    version=version,
    author="PAGO team",
    author_email="pago-dev@groupes.renater.fr",
    maintainer='Nicolas Barrier',
    maintainer_email='nicolas.barrier@ird.fr',
    description="Python Physical Analysis of Gridded Ocean",
    license="CeCILL", 
    keywords="ocean; grid model; transport; sections; budgets; heat; freshwater; ROMS; NEMO; IPSL; CNRM; GFDL",
    include_package_data=True,
    url="https://sourcesup.renater.fr/pago",
    packages=find_packages(),
    install_requires=['docutils>=0.12',
                      'sphinx>=1.3.1',
                      'pylint>=1.4.2',
                      'pyenchant>=1.6.6',
                      'pep8>=1.6.2',
                      'pyflakes>=0.9.2',
                      'check-manifest>=0.25',
                      'numpy>=1.9',
                      'netCDF4>=1.1', 
                      'matplotlib>=1.4',
                      'basemap>=1.0',
                     ],
    requires=['numpy(>=1.9.2)',
              'matplotlib(>=1.43)',
              'basemap(>=1.0.7)',
              'netcdf4(>1.1.9)',
             ],
    long_description=open('README.md').read(),

    classifiers = [
        #"Development Status :: 5 - Production/Stable",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: Free To Use But Restricted",  # .. todo::
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 2.7",
    ],

    # ++ test_suite =
    # ++ download_url
    platforms=['linux', 'mac osx'],
    scripts = ['pypago/bin/make_grid.py', 
               'pypago/bin/make_gridsec.py',
               'pypago/bin/make_areas.py',
               'pypago/bin/make_coords.py',
               'pypago/guis/gui_sections_edition.py',
               'pypago/guis/gui_grid_model.py']
)
#
# todo
# ====
#
# validation with :samp:` python setup.py check` and
# :samp:`python setup.py test`
#
# check evolution of 'sphinxcontrib-bibtex after >=0.3.3
#
# identify scripts to be place in bin/, to be documented with man
#
# add requirement for developer (doc builder, commiter) : doc8,
# hunspell (with its hunspell-dict-en_us dictionary)
#
# think about uninstall
#
# requires or install_requires
#
# make a difference between requirement for doc production and execution
#
# make it work for pip install
#
# make it work for pip register
#
# license
#
# complete with pipy classifiers
# cf. http://pypi.python.org/pypi?%3aaction=list_classifiers
#
# check for the next version of setuptools if sphinx latexpdf builder available
#
# link to svn repository
#
# EVOLUTIONS
# ==========
#
# - fplod 20151030T145305Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * version uniquely defined in :file:`VERSION`
#     thanks to https://packaging.python.org/en/latest/single_source_version/
#     and others
#   * project name not anymore pypago. replace by PAGO
#
# - fplod 20151028T171338Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * move usage to packaging.rst
#
# - fplod 20150918T145617Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * add requirement pyenchant
#
# - fplod 20150914T085629Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * give up usage of sphinxcontrib-bibtex
#
# - fplod 20150811T123543Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * usage of find_packages to avoid hard coded path of pypago sub-directories
#   * add |matlab| files see :file:`matlab/__init__.py` and file:`MANIFEST.in`
#
# - fplod 20150810T113415Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * add wheel distribution making (this will replace egg)
#     thanks to https://hynek.me/articles/sharing-your-labor-of-love-pypi-quick-and-dirty/
#     not so dirty ...
#
# - fplod 20150807T155911Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * get inspired by pylint version of setup files
#     cf. https://bitbucket.org/logilab/pylint/src
#
# - fplod 20150806T125248Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * add requirement of pep8, pyflakes, pylint
#
#     .. note:: pep8 bug
#
#        pep8 1.6.2 do not detect syntax error in
#        https://sourcesup.renater.fr/scm/viewvc.php/trunk/pypago/pypago.py?view=markup&revision=53&root=pago&pathrev=53  # pylint: disable=line-too-long
#
# - fplod 20150805T120655Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * more precise requirements
#   * fix typo for sphinxcontrib.bibtex (module vs package name)
#
# - fplod 20150804T160106Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * better pylint mark
#   * read README.txt to fill long_description
#     see output of :samp:`python setup.py --long_description`
#
# - fplod 20150731T102757Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * new topics, add audience thanks to inspiring
#     https://forge.ifremer.fr/scm/viewvc.php/trunk/setup.py?view=markup&root=vacumm
#
# - fplod 20150730T160440Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * thanks to MANIFEST.in, HTML and latexpdf produced by sphinx
#     are in the tar.gz file
#
# - fplod 20150730T123052Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * apply https://docs.python.org/2/distutils/examples.html
#     recommendations for packages
#   * add command for production of sphinx html builder (see also setup.cfg)
#
# - fplod 20150728T154914Z guest-242.locean-ipsl.upmc.fr (Darwin)
#
#   * 2d round - missing doc, test
#
# - fplod 20150618T114221Z callisto.locean-ipsl.upmc.fr (Linux)
#
#   * 1st draft thanks to
#     http://pythonhosted.org/an_example_pypi_project/setuptools.html
