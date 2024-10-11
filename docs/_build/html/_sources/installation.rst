.. _install:

**************
Installation
**************

.. _quick-start:

===========
Quick Start
===========

In the `OpenKet GitHub repository <https://github.com/pbbmx/openket.git>`_ you should be able to download OpenKet package.
Then, in the root directory of the proyect execute:


.. code-block:: bash

   pip install .

It is not recommended to install any packages directly into the system Python environment;
consider using pip or conda virtual environments to keep your operating system space clean,
and to have more control over Python and other package versions.

Finally, yout must include in your .py file the module, for example:

.. code-block:: python

   from openket.core.diracobject import Ket