# Wrapper to import correct version of ipython_directive
try:
    from ipython_directive_0p11 import *
except ImportError:
    from ipython_directive_0p10 import *
