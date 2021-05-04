
Ghost detection tool for JWST NIRISS data
=========================================

Author: Takahiro Morishita

.. image:: ./figure/demo.png

Purpose
-------

NIRISS WFSS and direct imaging modes are known to produce optical ghosts when there are bright sources in the observed field of view, as summarized `here <https://jwst-docs.stsci.edu/near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-gr150-grisms#NIRISSGR150Grisms-Ghosts>`__.
Users are advised to apply this code to final NIRISS images, so that they can avoid confusing real sources with optical ghosts in their final products.

By providing _i2d.fits image and source catalog to this script, you can:

- identify possible ghost in the input image
- get a source catalog with a new flag column, ``is_this_ghost``.


Usage
-----

Ghost detection in a calibrated image (i.e. _i2d.fits from the JWST pipeline).

.. code-block:: bash

    python detect_ghost_image3.py [image] [catalog]

Optional arguments:

- --rlim: Search radius from the predicted coordinates of a ghost, in pixel.
- --frac_ghost: Fraction flux of a ghost compared to the source.
- --o: Output directory. Default is set to the working directory.
- --f_mirage: Is the input image created by Mirage? If not (i.e. on-sky data), set this False.
- --keyword_flux: Column name for flux in ``catalog``. Default is source_sum (one that comes with photutils.).


Determine ghost axis point (GAP) coordinates based on a calibrated image (i.e. _i2d.fits from the JWST pipeline).

.. code-block:: bash

    python get_gap.py [image] [catalog]

Optional arguments:

- --nmc: Number of MCMC iterations (3000 default).
- --nwalkers: Number of walkers (20 default).
- --check_flux: Use flux ratio for posterior calculation (True default).


Appendix: Simulation of ghosts in a NIRISS scene
------------------------------------------------

Ghost simulation module has been implemented in Mirage. Take a look at the documentation `here <https://mirage-data-simulator.readthedocs.io/en/latest/ghosts.html>`__.

An example notebook is also available from this repository, which demonstrates a case with a custom input file for ghosts.


Result with a custom fits stamp:

.. image:: ./figure/demo_custom.png
