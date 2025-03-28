.. _install:

**************
Installation
**************

.. _prerequisites:

==============
Prerequisites
==============

* Python 3.6 or higher. Check with

.. code-block:: bash

   python --version

* Git (to clone the repository).

* Updated pip

.. code-block:: bash

   pip install --upgrade pip


.. _quick-start:

===========
Quick Start
===========

You should be able to download openket from `OpenKet GitHub repository <https://github.com/pbbmx/openket.git>`_.
To install, first clone the repository:

.. code-block:: bash

   git clone https://github.com/pbbmx/openket.git
   cd openket

Then, in the root directory of the proyect execute:

.. code-block:: bash

   pip install -r requirements.txt
   pip install .


It is not recommended to install any packages directly into the system Python environment;
consider using ``pip`` or ``conda`` virtual environments to keep your operating system space clean,
and to have more control over Python and other package versions.

Finally, verify installation:

.. code-block:: python

   import openket
   print(openket.__version__)



=====================
General Requirements
=====================

The following packages are required:

.. cssclass:: table

+----------------+--------------+
| Package        | Version      |
+================+==============+
| **Python**     | 3.6+         |
+----------------+--------------+
| **mpmath**     | 1.3+         |
+----------------+--------------+
| **NumPy**      | 2.1+         |
+----------------+--------------+
| **Sympy**      | 1.13+        |
+----------------+--------------+
| **SciPy**      | 1.10+        |
+----------------+--------------+