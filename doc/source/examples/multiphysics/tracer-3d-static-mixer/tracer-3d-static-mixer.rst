================================================
Tracer Transport Through Helical Static Mixer
================================================

This example is the follow up of the `3D Static Mixer Using RBF Sharp-Immersed Boundary` example in ``/examples/sharp-immersed-boundary/3d-rbf-static-mixer/``.

As their name implies, static mixers' main purpose is mixing fluids. In this example, we build from the previous example by activating the tracer physics and reusing the velocity solution to advect the tracer at a reduced computing cost.


+-----------------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/static_mixer_stl_casing_arrows.png                                                                      |
|     :align: center                                                                                                          |
|     :width: 800                                                                                                             |
|     :name: Surface grid representation of a helix static mixer with its casing.                                             |
|                                                                                                                             |
|     Surface grid representation of a helix static mixer.                                                                    |
|                                                                                                                             |
+-----------------------------------------------------------------------------------------------------------------------------+

----------------------------------
Features
----------------------------------

- Solvers: ``lethe-fluid-sharp``
- Transient transport (and mixing) of tracer
- Two step simulation: first with active fluid dynamics, second with only tracer active (reusing fluid dynamics solution)
- Complex surface-grid-defined static solid modeled with the sharp immersed boundary method


----------------------------
Files Used in This Example
----------------------------
All files mentioned below are located in the example's folder (``/examples/sharp-immersed-boundary/3d-rbf-static-mixer/``).

* RBF geometry file: ``/lethe_sharp_simulation/RBF_helix.output``. The extension is ``.output``, because it was named from the `bitpit <https://github.com/optimad/bitpit>`_ perspective;
* Composite geometry file: ``/lethe_sharp_simulation/mixer_long.composite``
* Parameter file: ``/lethe_sharp_simulation/flow_in_long_mixer.prm``


-----------------------
Description of the Case
-----------------------

In this example, we build upon the static mixer example by adding the usage of the tracer physics. This allows us to get a sense of the mixing capabilities of the helical mixer.

The solid to be used and its associated files are the same as in the previous example.


----------------------------
Dealing with Two Simulations
----------------------------

Two simulations are run in this example, but they can both be launched from a single command by running ``script_tracer.sh``. The first simulation has both fluid dynamics and tracer physics active (although nothing happens with the tracer), and the second simulation deactivates fluid dynamics and solves only the tracer physics, with activating an inlet Dirichlet boundary condition.

The script first clears the terminal, then removes the output directories from previous simulations. It then runs the first simulation, using MPI. Here we use 12 processes but the user can change the script to fewer or more, as their equipment allows. The whole folder is then copied to another name, and it will serve as basis for the second simulation. This second simulation is then launched in a similar fashion.

.. code-block:: text

    #!/bin/bash
    # RELEASE
    clear

    rm -rf output_cfd
    rm -rf output_cfd_tracer
    mpirun -np 12 lethe-fluid-sharp flow_in_long_mixer_no_tracer.prm

    cp -r output_cfd output_cfd_tracer
    mpirun -np 12 lethe-fluid-sharp flow_in_long_mixer_split_tracer.prm



---------------
Parameter File
---------------

Simulation Control
~~~~~~~~~~~~~~~~~~



**First simulation** The parameter stay the same as in the previous example, except the ``output path`` which is set to ``output_cfd`` for clarity.

.. code-block:: text

    subsection simulation control
      set method      = bdf1
      set time end    = 40e-4
      set time step   = 1e-4
      set output path = ./output_cfd/
      set output name = output
    end

**Second simulation** The ``time end`` is increase to ``5000`` to give enough time for the tracer to pass through the whole domain. The ``time step`` is increased to ``1`` to reduce computing time. The ``output path`` is changed to ``output_cfd_tracer`` to signify that the tracer physics is active.

.. code-block:: text

    subsection simulation control
      set method      = bdf1
      set time end    = 5000
      set time step   = 1
      set output path = ./output_cfd_tracer/
      set output name = output
    end

Physical Properties
~~~~~~~~~~~~~~~~~~~

From the previous example, we added the tracer diffusivity. It has no physical meaning in this context.

.. code-block:: text

    subsection physical properties
      subsection fluid 0
        set kinematic viscosity = 0.01
        set tracer diffusivity  = 0.01
      end
    end


Multiphysics
~~~~~~~~~~~~

**First simulation** Both ``fluid dynamics`` and ``tracer`` physics are enabled. ``use time average velocity field`` is added and explicitly disabled to ensure that the ``average velocity field`` mentioned in the post-processing section is not used yet.

.. code-block:: text

    subsection multiphysics
      set fluid dynamics                  = true
      set tracer                          = true
      set use time average velocity field = false
    end

**Second simulation** From the assumption that steady state is reached for the pressure and velocity fields, ``fluid dynamics`` is deactivated to reduce computing costs. ``use time average velocity field`` is now enabled to signify that the velocity solution to use is now the pre-computed average.

.. code-block:: text

    subsection multiphysics
      set fluid dynamics                  = false
      set tracer                          = true
      set use time average velocity field = true
    end


Post-Processing
~~~~~~~~~~~~~~~

``calculate average velocities`` is enabled, but only for the last three time steps (after ``initial time = 38e-4``). This ensures that the average only applies when the steady-state is reached.

.. code-block:: text

    subsection post-processing
      set verbosity               = verbose
      set calculate pressure drop = true
      set calculate flow rate     = true
      set inlet boundary id       = 0
      set outlet boundary id      = 1

      set calculate average velocities  = true
      set initial time                  = 38e-4
    end


Restart
~~~~~~~

This section is at the core of the two-simulation mechanism. It saves checkpoints in the first simulation, and uses the last checkpoint to restart with new parameters that take advantage of the tracer physics.


**First simulation** The ``checkpoint`` is enabled, at every time step.

.. code-block:: text

    subsection restart
      set checkpoint = true
      set frequency  = 1
      set filename   = restart
      set restart    = false
    end

**Second simulation** The ``checkpoint`` is now disabled, but ``restart`` is enabled to allow using the first simulation as a starting point.

.. code-block:: text

    subsection restart
      set checkpoint = false
      set frequency  = 1
      set filename   = restart
      set restart    = true
    end



Boundary Conditions Tracer
~~~~~~~~~~~~~~~~~~~~~~~~~~

The boundary conditions for the tracer physics need to be defined.

#. Only 1 of them is defined, as the outlet boundary condition is assumed implicitly to apply no constraint.
#. The same idea goes for the wall boundary conditions: we want the concentration gradient to be 0, which is the natural condition of finite elements.
#. The inlet boundary condition is of Dirichlet type and needs to be changed between the first (``set Function expression = 0``) and second (``set Function expression = y>0?1:0``)simulations. The second boundary condition will allow testing of the mixing effects in the `Y` direction.

.. code-block:: text

  subsection boundary conditions tracer
    set number = 1
    subsection bc 0
      set id   = 0
      set type  = dirichlet
      subsection dirichlet
        # To use in the first simulation
        set Function expression = 0
        # To use in the second simulation
        set Function expression = y>0?1:0
      end
    end
  end



--------
Results
--------

TODO

+-----------------------------------------------------------------------------------------------------------------------------+
|  .. figure:: images/long_static_mixer_medium_thick_p_v.png                                                                  |
|     :align: center                                                                                                          |
|     :width: 800                                                                                                             |
|     :name: Streamlines in the static mixer colored by velocity magnitude and pressure                                       |
|                                                                                                                             |
|     Streamlines in the static mixer colored by velocity magnitude and pressure                                              |
|                                                                                                                             |
+-----------------------------------------------------------------------------------------------------------------------------+


