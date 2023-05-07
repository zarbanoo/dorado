# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa

from .target import *
from .ceres import *
from .stack import *
from .core import *
from .timeseries import *
from .dorphot import *
# from .graf import * # TODO decide if this is needed

# from .fournax import *

# ----------------------------------------------------------------------------

# __all__ = []
# from .config import *
# from .example_mod import *   # noqa
# Then you can be explicit to control what ends up in the namespace,
# __all__ += ['do_primes']   # noqa
# or you can keep everything from the subpackage with the following instead
# __all__ += example_mod.__all__


# .. warning::
# :mod:`word`