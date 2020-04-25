# Control-Box RST

## Status

ROS package:
- ROS Melodic (*melodic-devel*): [![Melodic Status](http://build.ros.org/buildStatus/icon?job=Mdev__control_box_rst__ubuntu_bionic_amd64)](http://build.ros.org/job/Mdev__control_box_rst__ubuntu_bionic_amd64/)


## Documentation 

Build and installation instructions as well as further documentation are provided in the [project wiki](https://github.com/rst-tu-dortmund/control_box_rst/wiki).

## Authors

* Christoph Rösmann <christoph.roesmann@tu-dortmund.de>
* Maximilian Krämer <maximilian.kraemer@tu-dortmund.de>

## Citing the Software

*Since a lot of time and effort has gone into the development, please cite at least one of the following publications if you are using the software for published work.*

**Standard MPC and Hypergraph**

- C. Rösmann, M. Krämer, A. Makarow, F. Hoffmann und T. Bertram: Exploiting Sparse Structures in Nonlinear Model Predictive Control with Hypergraphs, IEEE/ASME International Conference on Advanced Intelligent Mechatronics (AIM), New Zealand, July 2018. 

**Time-Optimal MPC and Hypergraph**

- C. Rösmann: [Time-optimal nonlinear model predictive control, Direct transcription methods with variable discretization and structural sparsity exploitation](http://dx.doi.org/10.17877/DE290R-20283). Dissertation, Technische Universität Dortmund, Oct. 2019.

**Uniform Grid Time-Optimal MPC**

- C. Rösmann, F. Hoffmann und T. Bertram: Timed-Elastic-Bands for Time-Optimal Point-to-Point Nonlinear Model Predictive Control, European Control Conference (ECC), Austria, July 2015.
- C. Rösmann, F. Hoffman und T. Bertram: Convergence Analysis of Time-Optimal Model Predictive Control under Limited Computational Resources, European Control Conference (ECC), Denmark, June 2016.

**Non-Uniform Grid Time-Optimal MPC**

- C. Rösmann, A. Makarow, F. Hoffmann und T. Bertram: Sparse Shooting at Adaptive Temporal Resolution for Time-Optimal Model Predictive Control, IEEE Conference on Decision and Control (CDC), Australia, December 2017.
- C. Rösmann, A. Makarow, F. Hoffmann und T. Bertram: Time-Optimal Nonlinear Model Predictive Control with Minimal Control Interventions, IEEE Conference on Control Technology and Applications (CCTA), Hawai'i, August 2017.

<a href="https://www.buymeacoffee.com/croesmann" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/lato-orange.png" alt="Buy Me A Coffee" height="31px" width="132px"></a>

## License

Copyright (c) 2020,
TU Dortmund - Institute of Control Theory and Systems Engineering.
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Some third-party dependencies are included that are licensed under different terms:
 - *Eigen*, MPL2 license, http://eigen.tuxfamily.org

Optional included third-party dependencies (selected during configuration)
 - *gRPC*, Apache License 2.0, https://grpc.io
 - *protobuf*, in parts BSD, https://developers.google.com/protocol-buffers 
 - *qcustomplot*, GPLv3, https://www.qcustomplot.com
 - *googletest*, BSD-3-Clause, https://github.com/google/googletest
 - *yaml-cpp*, MIT License, https://github.com/jbeder/yaml-cpp

Optional third-party dependencies (optional linking)
 - *Ipopt*, EPL 1.0, https://github.com/coin-or/Ipopt
 - *OSQP*, Apache 2.0, https://osqp.org/
 - *Qt*, different licensing options available, https://www.qt.io/





