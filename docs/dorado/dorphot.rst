Basic Photometry and Dorphot
====================================================

Base photometry classes in dorado act as an interface between the core of dorado and the observational
data stored in the stacks of a :class:`~dorado.Ceres` instance. There are many things one may wish to do
to such data depending on its level of calibration and the target of study; thus, dorado again takes a
modular approach.


The aicophot class
----------------

The :class:`~dorado.aicoPhot` class is a dorphot class that like :class:`~dorado.aico_reader` was created
to be both an example class for others to model their extension classes after, as well as serving the data
processing needs of the Allan I. Carswell Observatory (AICO) as discussed in :doc:`reader</dorado/reader>`.

If we already have a :class:`~dorado.Target` (toi) and wish to operate on a :class:`~dorado.Ceres` instance
containing raw data in a Johnson red filter labelled 'R' in the '.fits' header or HDU.header we should start by
setting our dorphot instance:

.. code:: python
    Dorado.dorphot = dorado.aicoPhot()

The next logical step is to calibrate the data for flatfield and bias effects:

.. code:: python
    Dorado.dorphot.calibrate('ceres_name','R')

The :class:`~dorado.aico_reader` should have located the bias frames for the observation and constructed
a base bias frame using 'mkbias', it likewise should have constructed base flat frames for each 
:class:`~dorado.Stack`. 'What about Dark frames?' 







.. code:: python
    # align the series red filter data
    Dorado.dorphot.align('ceres_name', 'R')
    # perform differential photometry on the target in the red filter 
    # data using the control as a reference
    Dorado.dorphot.apPhot('ceres_name', 'R', toi, control)


Previous: :doc:`Ceres</dorado/target>` || Next: :doc:`Dorphot</dorado/timeseries>`