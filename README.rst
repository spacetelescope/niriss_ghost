
Ghost identification tool for JWST NIRISS data
==============================================

Author: Takahiro Morishita

.. image:: ./figure/demo.png

Purpose
-------

NIRISS WFSS and direct imaging modes are known to produce optical ghosts when there are bright sources in the observed field, as summarized `here <https://jwst-docs.stsci.edu/near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-gr150-grisms#NIRISSGR150Grisms-Ghosts>`__.
Users can identify possible ghosts and mask them out by applying the script provided in this repository.

By providing _cal.fits image (IMAGE2 product) and source catalog to the script, you can:

- identify possible ghosts in the input image
- get a source catalog with a new flag column, ``is_this_ghost``.
- change values of the pixles of the detected ghosts in the DQ array to ``DO_NOT_USE``, so the pipeline IMAGE3 step will ignore these contaminated pixels.


.. figure:: ./figure/DQ_masking.png
    :width: 800
    :align: center

    *How this works. Ghost detection, and then DQ tweaking here, can be applied to IMAGE2 products 
    before these go into the IMAGE3 step.*



Installation
------------

.. code-block:: bash

    python setup.py install

or 

.. code-block:: bash

    pip install niriss_ghost

Usage
-----

Ghost detection in a calibrated image (i.e. _cal.fits from the JWST pipeline). See ./example/Photometry_IMAGE2_product.ipynb for how to get photometry catalog and segmentation map on _cal.fits.

.. code-block:: bash

    python detect_ghost_image2.py [image] [catalog]

Input arguments:

- image: Input image(s). If more than one, then comma separated string e.g., image1_cal.fits,image2_cal.fits
- catalog: Input catalog(s). If more than one, then comma separated string e.g., image1.cat,image2.cat

Optional arguments:

- --rlim: Search radius from the predicted coordinates of a ghost, in pixel.
- --frac_ghost: Fraction flux of a ghost compared to the source.
- --o: Output directory. Default is set to the working directory.
- --f_mirage: Is the input image created by Mirage? If not (i.e. on-sky data), set this False.
- --keyword_id: Column name for source id in ``catalog``. Default is id.
- --keyword_flux: Column name for flux in ``catalog``. Default is source_sum.
- --keyword_xcent: Column name for x-pixel-position in ``catalog``. Default is xcentroid.
- --keyword_ycent: Column name for y-pixel-position in ``catalog``. Default is ycentroid.
- --f_tweak_dq: Change DQ arrays of the positions of the detected ghosts. You need the segmentation map of the provided catalog (_seg.fits).
- --segmap: Segmentation map of the provided catalog, if f_tweak_dq==True. (Default: image.replace('.fits', '_seg.fits'))

Alternatively, you can run the script in your python script;

.. code-block:: bash

    from niriss_ghost import detect_ghost_image2
    list_images = ['image1_cal.fits']
    list_catalogs = ['image1.cat']
    detect_ghost_image2.run(list_images, list_catalogs)


Caveat
------

- This script currently supports only _cal.fits images.
- Due to recent changes in the photutils package, the column keywords used in this script may not match with those in the input catalog. 
If this happens, a warning will appear. Users may specify those keywords by using --keyword_* argumens (see above, `Optional arguments`).


Appendix: Simulation of ghosts in NIRISS scenes
-----------------------------------------------

Ghost simulation module has been implemented in Mirage. Take a look at the documentation `here <https://mirage-data-simulator.readthedocs.io/en/latest/ghosts.html>`__.

An example notebook is also available from this repository, which demonstrates a case with a custom input file for ghosts.


.. figure:: ./figure/demo_custom.png
    :width: 800
    :align: center

    *Result with a custom fits stamp.*