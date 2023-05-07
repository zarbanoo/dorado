Target handling and special target classes
====================================================

While looking at the night sky in general is a completely normal and natural
endevour, typically in astronomical data operations, one has a target or targets
in mind. 

The Target class
----------------

The :class:`~dorado.Target` class handles storing and interfacing with astronomical 
targets in dorado. This base class takes a target name and queries Simbad for a match;
Simbad then returns the targets coordinates. 

.. code:: python
    toi = dorado.Target('target name') # toi == target of interest
    print(toi.coords) # astropy.coordinates.SkyCoord()

    # or create a target within the dorado core 
    Dorado.mktrgt('target name')

There are many types of targets one may wish to work with, and some of these target types
may benifit from having additional metadata, the ability to store child targets (such as
a group of variable star targets within a glubular cluster target class), or 
methods specific to the target object type (such as surface brightness for extended targets 
like galaxies or nebulous regions). 

The :class:`~dorado.Target` class is thus a base class to be extended to suit a target types
needs in the future.



The Fournax class
-----------------

The :class:`~dorado.Fournax` class is the first built in extension to the :class:`~dorado.Target`
class. Fournax is an abbreviation of Fourier numerical astronomy extension,
 its name is a backronym styled to match the constellation 'fornax'. 

.. code:: python
    toi = dorado.Fournax('target name') # toi == target of interest
    print(toi.coords) # astropy.coordinates.SkyCoord()

This class aims to be a base for stellar target types whose brightness varies in time such as stars which are transited
by exoplanets or delta scuti variable stars. The current implementation favors intrinsic
variable stars and not stars whose brightness changes due to occultation like exoplanet transits.

Future updates of dorado will expand the diversity and functionality of Fournax and other target type
specific classes.




Previous: :doc:`Ceres</dorado/ceres>` || Next: :doc:`Dorphot</dorado/dorphot>`