.. _GettingStarted

*******************
Getting Started with Dorado
*******************

Install
-------


Install dorado by running the command: ``pip install dorado``.

    .. note:: The `dorado <https://pypi.org/project/dorado/>`_ project on `PyPI <https://pypi.org/>`_ . 

Pip should install all relevent dependencies for dorado into the python environment 
pip was run from. The dorado team recommends using the latest version of python unless otherwise 
mentioned for the best compatibility.


Dependencies
-------------

Dorado aims to be as lightweight as possible and utilize as little dependencies as it can. 

Currently DORADO relies on:  

1.  `numpy <http://www.numpy.org/>`_

2.  `matplotlib <https://matplotlib.org/>`_

3.  `astropy <https://www.astropy.org/index.html>`_

4.  `CCDPROC <https://ccdproc.readthedocs.io/en/latest/index.html#>`_ >> `CCDprocX <https://pypi.org/project/ccdprocx/>`_

5.  `scipy <https://www.scipy.org/>`_

6.  `photutils <https://photutils.readthedocs.io/en/stable/index.html>`_

7.  `astroquery <https://astroquery.readthedocs.io/en/latest/#>`_

8.  `astroalign <https://astroalign.quatrope.org/en/latest/?badge=latest>`_

9.  `tqdm <https://tqdm.github.io/>`_

10. `lightkurve <http://docs.lightkurve.org/>`_


    .. note:: CCDPROC requires the ``astroscrappy`` package for install. Currently, astroscrappy requires 
                Visual C++ redistributable for build and function. This install step has a variety of potential issues
                users can encounter so DORADO instead uses a fork of CCDPROC with astroscrappy removed known as ccdprocx.

Babies first Dorado script
==========================

Adding Dorado into your flow is fairly simple and can be broken into 3 main steps; 
import, initialize, and call. 

.. code:: python

        ## import
        import dorado

        ## initialize
        # create a target object 'target of interest' toi
        toi = 'target_name'
        Dorado.mktrgt(toi)
        # create a control target object
        control = 'control_target_name'
        Dorado.mktrgt(control)
        # create a series object
        Dorado.mkceres('2021-01-01+02', name = 'ceres01', target = toi)

        ## call
        # calibrate the series red filter data
        Dorado.dorphot.calibrate('ceres01','R')
        # align the series red filter data
        Dorado.dorphot.align('ceres01', 'R')
        # perform differential photometry on the target in the red filter 
        # data using the control as a reference
        Dorado.dorphot.apPhot('ceres01', 'R', toi, control)

        ## finish by saving
        # save the resulting ceres data 
        Dorado.savewrk('ceres01')
        # record the target results 
        Dorado.targets[Dorado.target_keys[toi]].record('ceres01')


Next: :doc:`Dorado Core</dorado/core>`
