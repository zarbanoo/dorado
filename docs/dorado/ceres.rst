Ceres and stack classes: Observational data handling
====================================================

Ceres
------

To interface with an observational series, Dorado offers the Ceres class. The Ceres class contains a list 
of Stack objects and the coresponding metadata including :class:`~dorado.Target` objects associated with
the campaign. As stated before, Astronomical data is collected and stored in many different ways and formats;
Dorado aims to be an open and community friendly pipeline tool, to do this, dorado needs a modular and 
personalizable way to integrate data into the dorado system. 

Ceres objects are by convention stored in the core dorado instance and referenced by key strings for ease of use.


The first step to working with a ceres object is to import dorado
.. code:: python

        ## import
        import dorado

You may want to store a target object associated with the campaign in the ceres instance. Lets quickly
initialize one.

.. code:: python
        ## initialize
        # create a target object 'target of interest' toi
        toi = 'target_name'
        Dorado.mktrgt(toi)

Next up comes initializing the ceres object into the core dorado instance. We will hand it an observation
date of `2021-01-01+02` which points dorado to a date folder with the same name in the target data folder.
In this case that folder has the path `home/.dorado/data/raw/2021-01-01+02` with `home/.dorado/data/raw/` 
being the default data path. We will nickname the ceres instance as `ceres_name` so we can later reference it 
by its key string. Lastly, we will append our target object into the ceres instance.

.. code:: python
        # create a series object
        Dorado.mkceres('2021-01-01+02', name = 'ceres_name', target = toi)

lastly, we wish to work with our new ceres instance. To calibrate the red filter data within our ceres 
instance we can use the following. 

.. code:: python
        ## call
        # calibrate the series red filter data
        Dorado.dorphot.calibrate('ceres_name','R')

.. note:: The above syntax for calibration is assuming that :class:`~dorado.aicoPhot` is set as the 
    active `dorado.dorphot` instance.



Stack
------

The :class:`~dorado.Stack` hosts idividual stacks of data in dorado such as data specific to a target
in a single filter. The reader class conventionally handles Stack instances by reading them into ceres 
instances and labelling them; conventionally, Stack key strings are generated from the filter name 
found in the data. In the above example, the filter name was `R` so we told the calibration method to 
calibrate the `R` stack in the `ceres_name` ceres instance.

In the future, more advanced documentation on the stack class will be available for those wishing to 
produce custom pipeline classes or work with data in a more fundamental way.

Previous: :doc:`Reader</dorado/reader>` || Next: :doc:`Targets</dorado/target>`