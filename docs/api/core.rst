Core Classes
============

Root-level classes and enumerations that form the foundation of the USTB
framework. These classes are defined at the repository root (not inside a
``+package`` folder).

Base Classes
------------

**uff** (``uff.m``)
   Base class for all UFF (Ultrasound File Format) objects. Provides common
   functionality for naming, serialization, and HDF5 read/write support.
   Inherits from ``handle``.

**process** (``process.m``)
   Base class for all USTB processing steps. Provides the ``go()`` method
   interface and common properties like ``channel_data``, ``scan``, and
   apodization settings.

**midprocess** (``midprocess.m``)
   Base class for mid-processing (beamforming) algorithms. Inherits from
   ``process``. Subclasses include ``midprocess.das``.

**postprocess** (``postprocess.m``)
   Base class for post-processing algorithms. Inherits from ``process``.
   All post-processors in the ``+postprocess`` package inherit from this class.

**preprocess** (``preprocess.m``)
   Base class for pre-processing algorithms. Inherits from ``process``.
   All pre-processors in the ``+preprocess`` package inherit from this class.

**pipeline** (``pipeline.m``)
   Processing pipeline that chains multiple mid-process and post-process steps.
   Provides the ``go()`` method which accepts a cell array of processing objects
   and executes them in sequence.

Enumerations
------------

**code** (``code.m``)
   Enumeration for selecting the beamformer implementation.

   - ``code.matlab`` -- Pure MATLAB implementation
   - ``code.mex`` -- MEX C implementation
   - ``code.matlab_gpu`` -- MATLAB GPU implementation
   - ``code.mex_gpu`` -- MEX CUDA GPU implementation

**dimension** (``dimension.m``)
   Enumeration for selecting the beamforming dimension.

   - ``dimension.transmit`` -- Beamform in transmit only
   - ``dimension.receive`` -- Beamform in receive only
   - ``dimension.both`` -- Beamform in both transmit and receive
   - ``dimension.none`` -- No beamforming

**spherical_transmit_delay_model** (``spherical_transmit_delay_model.m``)
   Enumeration for selecting the spherical transmit delay model.

   - ``spherical_transmit_delay_model.spherical`` -- Full spherical model
   - ``spherical_transmit_delay_model.hybrid`` -- Hybrid model (default)
   - ``spherical_transmit_delay_model.plane`` -- Plane wave approximation

Utility Functions
-----------------

**data_path** (``data_path.m``)
   Returns the path to the USTB data directory used for storing downloaded
   datasets.

**ustb_path** (``ustb_path.m``)
   Returns the root path of the USTB installation.
