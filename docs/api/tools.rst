tools Functions
===============

The ``tools`` functions provide utilities for signal processing,
visualization, downloading datasets, and other helper tasks.

Signal Processing
-----------------

.. autofunction:: tools.band_pass

.. autofunction:: tools.high_pass

.. autofunction:: tools.low_pass

.. autofunction:: tools.rf2iq

.. autofunction:: tools.estimate_frequency

Beamforming
-----------

.. autofunction:: tools.matlab_beamformer

.. autofunction:: tools.matlab_gpu_beamformer

.. autofunction:: tools.calculate_unified_delay_model

Scan Conversion
---------------

.. autofunction:: tools.scan_convert

.. autofunction:: tools.scan_convert_na

.. autofunction:: tools.scan_convert_all_frames

Measurement
-----------

.. autofunction:: tools.measure_contrast_ratio

.. autofunction:: tools.measure_contrast_circles

.. autofunction:: tools.dynamic_range_test

.. autofunction:: tools.power_spectrum

.. autofunction:: tools.power_spectrum2

Statistics
----------

.. autofunction:: tools.weigthed_mean

.. autofunction:: tools.weigthed_var

.. autofunction:: tools.weigthed_std

.. autofunction:: tools.uniform_fov_weighting

Data Management
---------------

.. autofunction:: tools.download

.. autofunction:: tools.hash

Visualization
-------------

.. autofunction:: tools.plot_circle

.. autofunction:: tools.histogram_match

.. autofunction:: tools.viridis

.. autofunction:: tools.magma

.. autofunction:: tools.plasma

.. autofunction:: tools.inferno

Utilities
---------

.. autofunction:: tools.rotate_points

.. autofunction:: tools.check_memory

.. autofunction:: tools.getAvailableMemory

.. autofunction:: tools.text_progress_bar

.. autofunction:: tools.workbar

.. autofunction:: tools.dialog_timeout

.. autofunction:: tools.sector_MLA

.. autofunction:: tools.makehtmldoc
