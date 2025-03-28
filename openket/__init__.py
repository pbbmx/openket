# -*- coding: utf-8 -*-

"""
openket.
"""

# Standard imports
import logging
import sys
import warnings
from typing import List, Optional

# Basic logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

# Logger for the package
logger = logging.getLogger(__name__)

# Package metadata
__version__ = "0.1.0"

# Check for Scipy
try:
    import scipy
except ImportError:
    warnings.warn("scipy not found: Time evolution functions will not work.")
else:
    del scipy

#Import modules
from .core.diracobject import Ket, Bra, Operator, CreationOperator, AnnihilationOperator
from .core.metrics import *
from .core.gates import *
from .core.evolution import *


# Startup message
logger.info(f"openket v{__version__} initialized successfully.")

# Dependency verification
def _check_dependencies() -> Optional[ImportError]:
    """
    Verifies that all required dependencies are installed.
    """
    try:
        import sympy
        import matplotlib
    except ImportError as e:
        logger.error("Missing required dependencies. Please install them.")
        return e
    return None

# Run dependency verification when the package is imported
_check_dependencies()