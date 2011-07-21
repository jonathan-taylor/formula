""" This file contains defines parameters for formula that we use to fill
settings in setup.py, the formula top-level docstring, and for building the
docs.  In setup.py in particular, we exec this file, so it cannot import formula
"""

# formula version information.  An empty _version_extra corresponds to a
# full release.  '.dev' as a _version_extra string means this is a development
# version
_version_major = 0
_version_minor = 1
_version_micro = 0
_version_extra = '.dev'
# _version_extra = ''

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
__version__ = "%s.%s.%s%s" % (_version_major,
                              _version_minor,
                              _version_micro,
                              _version_extra)

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: BSD License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

description  = 'Symbolic statistical formulae'

# Note: this long_description is actually a copy/paste from the top-level
# README.txt, so that it shows up nicely on PyPI.  So please remember to edit
# it only in one place and sync it correctly.
long_description = """
=======
Formula
=======

This package contains an implementation of symbolic statistical models in
Python.

The implementation is of a full algebra of factors and numerical variables in
statistical models.

It has some similarities to model formulae in R, but differs in that the rules
for the creation of design matrices are more completely defined by the model
algebra.

An earlier version of this package is in nipy_.

Mailing Lists
=============

Please see the developer's list here::

    http://mail.scipy.org/mailman/listinfo/nipy-devel

Code
====

You can find our sources and single-click downloads:

* `Main repository`_ on Github.
* Documentation_ for all releases and current development tree.
* Download as a tar/zip file the `current trunk`_.
* Downloads of all `available releases`_.

.. _nipy: http://nipy.org
.. _main repository: http://github.com/jonathan-taylor/formula
.. _documentation: http://github.com/jonathan-taylor/formula
.. _current trunk: http://github.com/jonathan-taylor/formula/archives/master
.. _available releases: http://github.com/jonathan-taylor/formula/downloads

License
=======

formula is licensed under the terms of the BSD license.  Please the COPYING file
in the formula distribution.
"""

# versions for dependencies
NUMPY_MIN_VERSION='1.2'
SYMPY_MIN_VERSION='0.7.0'

# Main setup parameters
NAME                = 'formula'
MAINTAINER          = "Jonathan Taylor"
MAINTAINER_EMAIL    = "nipy-devel@neuroimaging.scipy.org"
DESCRIPTION         = description
LONG_DESCRIPTION    = long_description
URL                 = "http://github.com/jonathan-taylor/formula"
DOWNLOAD_URL        = "http://github.com/jonathan-taylor/formula/downloads"
LICENSE             = "BSD license"
CLASSIFIERS         = CLASSIFIERS
AUTHOR              = "Jonathan Taylor"
AUTHOR_EMAIL        = "nipy-devel@neuroimaging.scipy.org"
PLATFORMS           = "OS Independent"
MAJOR               = _version_major
MINOR               = _version_minor
MICRO               = _version_micro
ISRELEASE           = _version_extra == ''
VERSION             = __version__
PROVIDES            = ["formula"]
REQUIRES            = ["numpy (>=%s)" % NUMPY_MIN_VERSION,
                       "sympy (>=%s)" % SYMPY_MIN_VERSION]
