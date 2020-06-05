^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package control_box_rst
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

0.0.7 (2020-06-05)
------------------
* Fixed issue: Wrong sparse Jacobian pattern in solver LM
* Stage functions: added dedicated method for state and time dependency
* Lsq form added to hybrid cost functions
* Added missing header for eigenvalue computation
* Contributors: Christoph Rösmann

0.0.6 (2020-05-24)
------------------
* Hybrid cost functions added: Minimum time and control/state quadratic form
* Discretization grids: The time difference is now initialized to dt_ref for proper reference caching
* Contributors: Christoph Rösmann

0.0.5 (2020-05-13)
------------------
* Changed minimum CMake version to 3.1
* Add default parameter to SolverIpopt::initialize() in case Ipopt is not found..
* Removed obsolete gmock include
* Fix setLastControlRef and graph consistency.
* Internal libgtest target renamed. Fixed compilation issue for yaml-cpp test
* Contributors: Christoph Rösmann, Maximilian Krämer

0.0.4 (2020-02-21)
------------------
* Removed exec_depend on catkin to allow the use of the same release for ros1 and ros2
* Contributors: Christoph Rösmann

0.0.3 (2020-02-21)
------------------
* Changed build_type to cmake in package.xml according to REP-136
* Set some optional findPackage cmake scripts to quiet
* Contributors: Christoph Rösmann

0.0.2 (2020-02-19)
------------------
* Changed ros dependency from eigen3 to eigen
* Contributors: Christoph Rösmann

0.0.1 (2020-02-19)
------------------
* First release
* Contributors: Christoph Rösmann
