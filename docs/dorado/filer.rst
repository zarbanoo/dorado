Getting to know Filer
======================


In Dorado, Filer is a class object that facilitates interactions between the user and dore Dorado operations
such as interfacing with the ``.dorado`` directory (more on that below). To get staarted with Filer, and thus Dorado,
import the Filer class and initialize an instance of Filer :

.. code:: python

        # import Clippy from Dorado
        from dorado.filer import Filer

        # initialize an instance
        clippy = Filer()

Here we named our instance of the Filer class 'clippy' after the infamous Microsoft office digital assistant
of the late 90's. As the user, you can name your instance whatever you please, such as ``filer`` or ``dortool``.
Filer is akin to an automatic file and directory assistant, Filer can create, locate, navigate, read, and manage astronomical data
directories to reduce user error and frusteration and increase consistency .

The Filer class you now have in your script will serve as the backbone of our example Dorado pipeline.
While Dorado is pretty flexible in terms of functionality and use, Dorado does apply some strict procedures
to keep things organized and universal; this is most apparent within Filer and manifest itself as the
``.dorado`` directory which will be created either upon install of Dorado or upon the first time Filer is 
initialized. The ``.dorado`` directory will be located within your home directory, the location of the home directory
varies with operating system. This directory when created has the following format:

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
date. In the event that an observing session lacked these calibration frames, Filer will consult these two directories
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

The ``/targets`` directory is a directory used to store target data used by Dorado for future use in th context of the 
:doc:`Target</dorado/target>` class object. 

The ``/logs`` folder is th future home of log files produced by Dorado during operation, this may either be during observations
or during subsequent data operations. These logs can then be used to debug pipelines and Dorado or used for future reference
as to what actions were performed on a dataset. 

The ``/cache`` folder is a place for Filer to temporarily store files to the hard disk while performing certain tasks such as uploading
frames to 'Astrometry.net' and downloading the results suring runtime. The cache may also be used to store temporary log files
or any other file objects that are temporary in nature and need a place to exist where it is easily acessible while also easy to remove
via cache cleaning commands.


Previous: :doc:`Getting Started</dorado/GettingStarted>` || Next: :doc:`Ceres</dorado/ceres>`



