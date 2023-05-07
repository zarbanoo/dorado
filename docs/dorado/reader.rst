Intro to the reader class system
================================= 

The first iteration of Dorado was known as Draco (Data reduction at the Allan Carswell Observatory) way back in 2018/2019.
It was specific to the way the Allan I. Carswell Observatory (AICO) handled undergraduate variable star research which was
previously handled using IRAF. Dorado grew as AICO grew and sought to accommodate the plethora of changes including data 
operations involving the newly installed Planewave CDK-1000 1m telescope. 

Dorado at current epoch uses a 'modular' data reader to create, locate, navigate, read, and manage astronomical data within 
the dorado directory. Users can create their own data reader specific to their data pipeline needs or use one of the prebuilt
data readers provided in the dorado package such as :class:`~dorado.aico_reader`. 

dirscan
--------

dirscan is used to search and sort the data in a desired directory and return the results to possibly be used in functions
like mkceres. In many cases dirscan will be called internally without the end user ever needing
to directly inteface with the function itself.


mkceres
-------

Astronomical data is collected and stored in many different ways and formats and for dorado to truly be an open and community 
friendly pipeline tool, dorado needs a modular and personalizable way to integrate data into the dorado system. One of the 
universal data objects in dorado is the :class:`~dorado.Ceres` class which handles data associated with a nights observation 
campaign, including individual stacks of data held in :class:`~dorado.Stack` data objects. The reader class facilitates the 
reading of the data into classes like :class:`~dorado.Ceres` and :class:`~dorado.Stack` where it can be operated on in a consistant
manor. 

mkflat and mkbias
------------------

mkflat and mkbias are other methods which may be present in a reader class depending on the desired pipeline. mkflat and mkbias
utilize the results of dirscan to produce a base flat frame and base bias frame which can further be stored and used in a
:class:`~dorado.Ceres` or :class:`~dorado.Stack` instance.


savewrk
-------

Just as data is read into dorado, the reader class is responsible for exporting data in an appropriate format. Different pipeline
models possible with dorado may utilize the savewrk method to export data in wildly different ways.

**Further information on readers is coming as the reader system matures.**

-------

The inbuilt aico_reader class 
==============================

The built in reader class included with the dorado package is :class:`~dorado.aico_reader`. This class handles 
the directory scanning logic specific to the file saving patterns used at AICO, and then searches for and sorts the files,
before reading into dorado. 



Previous: :doc:`Dorado Core</dorado/core>` || Next: :doc:`Ceres</dorado/ceres>`

