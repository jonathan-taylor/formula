# Init for formula
import os

from .info import __version__, long_description as __doc__

from .parts import Term, Factor, Formula
from .ancova import ANCOVA
from .convenience import terms, make_recarray

from .pkg_info import get_pkg_info as _get_pkg_info
get_info = lambda : _get_pkg_info(os.path.dirname(__file__))
