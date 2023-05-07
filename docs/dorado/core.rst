Getting to know the core of dorado
======================


Dorado sorta breaks the rules when it comes to python by elevating its existence in the python runtime/namespace 
to the same level as python itself. Once imported by any script in a python runtime, a single instance of dorado will
be accesable anywhere by any script, package or function without further import or passing. One can think of this similar
to how once python starts running, functions like print() and int() do not need to be imported or passed. For those curious,
this simply involved inserting an instance dorado into the current run of pythons 'built_in' namespace at import. 
Some programmers may be upset by this move and feel that dorado is 'polluting the namespace' but we feel that the benefit
of the efficiency and ease of use that this move grants is worth temporarily adding a single additional entry to 
the python namespace.

.. code:: python

        import dorado
    


While Dorado is pretty flexible in terms of functionality and use, Dorado does apply some strict procedures
to keep things organized and universal; this is most apparent within ``.dorado`` directory which will be created 
either upon install of Dorado or upon the first time Dorado is imported. The ``.dorado`` directory will be located 
within your home directory, the location of the home directory varies with operating system. 

This directory when created has the following format:

::

    .dorado
    ├── data          
    │   ├── wrk
    │   └── bias
    │   └── flats
    │   └── raw
    │   └── graphical
    │   └── projects
    │   └── targets        
    ├── logs
    ├── cache         
    │   └── astrometryNet


Lets look at the ``/data`` folder, raw astronomical data is stored in subfolders within the ``/raw`` folder.
Subfolders in the raw data directory should take a naming scheme based on the date of observation
``/YYYY-MM-DD+DD`` where the ``DD+DD`` portion denotes the date before and after local midnight during the 
observing night. The inside of the nights folder can take the form of one of two basic file structures; single
directory, where all files are stored within the folder with no subdirectories, and multi level directory,
where there is a folder for lights, bias, and flatfield frames. If multiple filters were used to observe, multilevel 
directory format allows for each filter to have its own subdirectory in the lights and flatfield directories.

The ``/bias`` and ``/flats`` directory store the base bias and flatfield images produced by Dorado for each relevent 
date. In the event that an observing session lacked these calibration frames, Dorado will consult these two directories
to retrieve the closest in date calibration frames to use for calibration. 

    .. note:: In the near future this functionality will be expanded upon to allow constraints on the max date separation 
        allowed between the retrieved calibration frames and the date of observation; as well as the ability to compare these 
        calibration frames over time to see how your imaging setup changes and thus allow you to refine your max date separation
        constraints.

The ```/wrk`` directory contains the working folders for your processed data, this is the default location for any
calibrated, aligned, or otherwise worked with images. Data saved will follow the ``/YYYY-MM-DD+DD`` subdirectory
naming scheme.

At the current epoch, the ``/graphical`` and ``/projects`` subdirectories are unused but in the near future ``/graphical`` will  
store any finalized stacked or colour merged images, graphs, or timeseries plots, and ``/projects`` will be a place to store data 
pertaining to ongoing projects that span more than one observing session and target. In the future the community will be consulted 
as to how best organize and interact with the project directory.

The ``/targets`` directory is a directory used to store target data used by Dorado for future use in the context of the 
:doc:`Target</dorado/target>` class object. 

The ``/logs`` folder is th future home of log files produced by Dorado during operation, this may either be during observations
or during subsequent data operations. These logs can then be used to debug pipelines and Dorado or used for future reference
as to what actions were performed on a dataset. 

The ``/cache`` folder is a place for Dorado to temporarily store files to the hard disk while performing certain tasks such as uploading
frames to 'Astrometry.net' and downloading the results suring runtime. The cache may also be used to store temporary log files
or any other file objects that are temporary in nature and need a place to exist where it is easily acessible while also easy to remove
via cache cleaning commands.

=============

Reference/API
=============

.. automodapi:: dorado.core

-------

Previous: :doc:`Getting Started</dorado/GettingStarted>` || Next: :doc:`Reader</dorado/reader>`





