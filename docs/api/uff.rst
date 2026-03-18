uff Data Classes
================

The ``uff`` data classes represent the Ultrasound File Format (UFF) and are
used to store ultrasound data, probes, scans, and related structures. All
data classes can be written to and read from UFF files.

Data Classes
------------

.. autoclass:: uff.channel_data
   :members:

.. autoclass:: uff.beamformed_data
   :members:

Probe Classes
-------------

.. autoclass:: uff.probe
   :members:

.. autoclass:: uff.linear_array
   :members:

.. autoclass:: uff.curvilinear_array
   :members:

.. autoclass:: uff.matrix_array
   :members:

.. autoclass:: uff.curvilinear_matrix_array
   :members:

Scan Classes
------------

.. autoclass:: uff.scan
   :members:

.. autoclass:: uff.linear_scan
   :members:

.. autoclass:: uff.linear_scan_3D
   :members:

.. autoclass:: uff.sector_scan
   :members:

.. autoclass:: uff.sector_scan_3D
   :members:

Wave and Pulse
--------------

.. autoclass:: uff.wave
   :members:

.. autoclass:: uff.pulse
   :members:

.. autoclass:: uff.wavefront
   :members:

Support Classes
---------------

.. autoclass:: uff.apodization
   :members:

.. autoclass:: uff.phantom
   :members:

.. autoclass:: uff.point
   :members:

.. autoclass:: uff.transform
   :members:

.. autoclass:: uff.window
   :members:

Functions
---------

.. autofunction:: uff.read_object

.. autofunction:: uff.write_object

.. autofunction:: uff.version

.. autofunction:: uff.index
