""" Test R utilities """

from ..rutils import checkfor_r, r_cmd, r_check_library

from nose.tools import assert_true, assert_false, assert_equal, assert_raises
from nose.plugins.skip import SkipTest

def setup():
    if not checkfor_r():
        raise SkipTest('Need R for tests')


def test_check_library():
    # library checks
    assert_true(r_check_library('base'))
    assert_false(r_check_library('i_prefer_python'))

