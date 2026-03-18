postprocess
===========

A postprocess modifies beamformed data. This includes adaptive beamforming,
compounding, and image enhancement techniques.

| **Input:** ``beamformed_data`` |rarr| **Output:** ``beamformed_data``

Compounding
-----------

.. autoclass:: postprocess.coherent_compounding
   :members:

.. autoclass:: postprocess.incoherent_compounding
   :members:

Adaptive Beamforming
--------------------

.. autoclass:: postprocess.coherence_factor
   :members:

.. autoclass:: postprocess.generalized_coherence_factor
   :members:

.. autoclass:: postprocess.generalized_coherence_factor_OMHR
   :members:

.. autoclass:: postprocess.phase_coherence_factor
   :members:

.. autoclass:: postprocess.capon_minimum_variance
   :members:

.. autoclass:: postprocess.eigenspace_based_minimum_variance
   :members:

.. autoclass:: postprocess.delay_multiply_and_sum
   :members:

.. autoclass:: postprocess.simplified_delay_multiply_and_sum
   :members:

.. autoclass:: postprocess.short_lag_spatial_coherence
   :members:

Displacement Estimation
-----------------------

.. autoclass:: postprocess.autocorrelation_displacement_estimation
   :members:

.. autoclass:: postprocess.modified_autocorrelation_displacement_estimation
   :members:

Image Enhancement
-----------------

.. autoclass:: postprocess.non_local_means_filtering
   :members:

.. autoclass:: postprocess.wiener
   :members:

.. autoclass:: postprocess.median
   :members:

.. autoclass:: postprocess.scan_converter
   :members:

Gray Level Transforms
---------------------

.. autoclass:: postprocess.gray_level_transform
   :members:

.. autoclass:: postprocess.polynomial_gray_level_transform
   :members:

.. autoclass:: postprocess.scurve_gray_level_transform
   :members:

Utilities
---------

.. autoclass:: postprocess.max
   :members:

.. autoclass:: postprocess.stack
   :members:
