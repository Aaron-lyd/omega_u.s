# Neutral density variables
## Software for creating neutral density variables

### gamma^SCV

The pressure-invariant neutral density variable (Lang et al., 2020) which is created by using SCV (submesoscle coherent vortices) methods  to connect the observation bottle with th reference data.

### gamma^n

Neutral density based on Jackett and McDougall (1997).

# Neutral Surfaces
## Software for approximately neutral surfaces and geostrophic streamfunctions

### Approximately Neutral Surfaces

In the absence of irreversible mixing, fluid parcels in the ocean move such that they are always at their level of neutral buoyancy with their environment.  This constrains the direction of flow lie in a plane, called the neutral tangent plane. This plane is well-defined at every point in the ocean.  A neutral surface is an extensive 2D surface that is everywhere aligned with the neutral tangent plane (McDougall 1987).  However, neutral surfaces are not well-defined 2D surfaces -- a consequence of the non-linear equation of state for the density of seawater, and the resulting non-zero neutral helicity in the ocean (McDougall and Jackett 1988).  Hence, physical oceanographers craft approximately neutral surfaces, which are well-defined extensive 2D surfaces that are everywhere approximately aligned with the neutral tangent plane.  The following approximately neutral surfaces may be calculated with this software package. 

#### Omega surfaces

Omega surfaces (Klocker et al., 2009;  Stanley et al. 2021) are highly accurate approximately neutral surfaces, formed from an iterative procedure that solves a global least squares problem to minimize the neutrality error. 

#### Omega_s surfaces

The ANS that formed from an iterative procedure that solves a global least squares problem to minimize the slope error. 

#### Omega_s^2 surfaces

The ANS that formed from an iterative procedure that solves a global least squares problem to minimize the fictitious diapyncal diffusivity. 

#### Omega_u.s surfaces

The ANS that formed from an iterative procedure that solves a global least squares problem to minimize the square of the spurious diasurface velocity.

#### Omega_u.sTz surfaces

The ANS that formed from an iterative procedure that solves a global least squares problem to minimize the square of the spurious diasurface velocity times the veritial gradient of the temperature.

#### Omega_u.sSz surfaces

The ANS that formed from an iterative procedure that solves a global least squares problem to minimize the square of the spurious diasurface velocity times the veritial gradient of the salinity.

#### Omega_u.s+s^2 surfaces

The combination of Omega_u.s and Omega_s^2.



## Contents:
- `./GSW_3_06/                        `- GSW toolbox V3.6 based on TEOS-10
- `./gamma/                           `- create gamma^SCV neutral density variables
- `./lib_add/                         `- libraries
- `./run/                             `- example scripts
- `./omega-surface/                   `- create omega series of surfaces
- `./README.md                        `- this file

## Requirements:
MATLAB 2016b or higher (tested on 2017b, 2018b, and 2020a) with the Optimization Toolbox


## Installation:
1. To run the code, neutral-surfaces (by geoffstanley) repository needs to be installed firstly
2. Run "addpath(genpath(<<The path to this repository>>))"



## References:

Klocker, A., McDougall, T.J., Jackett, D.R., 2009. A new method for forming approximately neutral surfaces. Ocean Science 5, 155–172. https://doi.org/10.5194/os-5-155-2009

Lang, Y., Stanley, G. J., McDougall, T. J., & Barker, P. M. (2020). A pressure-invariant neutral density variable for the world's oceans. Journal of Physical Oceanography, 50(12), 3585-3604.

Lang, Y., Stanley, G. J., & McDougall, T. J., (2023). Spurious dianeutral advection and methods for its minimization. Journal of Physical Oceanography.

McDougall, T.J., 1987. Neutral Surfaces. Journal of Physical Oceanography. https://doi.org/10.1175/1520-0485(1987)017<1950:NS>2.0.CO;2

McDougall, T.J., Jackett, D.R., 1988. On the helical nature of neutral trajectories in the ocean. Progress in Oceanography 20, 153–183. https://doi.org/10.1016/0079-6611(88)90001-8

Montgomery, R., 1937. A suggested method for representing gradient flow in isentropic surfaces. Bull. Amer. Meteor. Soc 18, 210–212.

Stanley, G.J., 2019a. Neutral surface topology. Ocean Modelling 138, 88–106. https://doi.org/10.1016/j.ocemod.2019.01.008

Stanley, G. J., McDougall, T. J., & Barker, P. M. (2021). Algorithmic improvements to finding approximately neutral surfaces. Journal of Advances in Modeling Earth Systems, 13, e2020MS002436. https://doi.org/10.1029/2020MS002436 

Wüst, G., 1935. The stratosphere of the Atlantic ocean. Scientific Results of the German Atlantic Expedition of the Research Vessel “Meteor” 1925–27 6.

## Copyright:
MIT License

Copyright (c) 2023 Yandong Lang, Geoff Stanley, and Trevor McDougall

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:  

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  

Author(s) : Yandong Lang and Geoff Stanley

Email     : yandong.lang@unsw.edu.au; gstanley@uvic.ca

Email     : aaronlangyandong@gmail.com; geoffstanley@gmail.com
