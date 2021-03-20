# Copyright (C) 2020-2021 by Andrew Hoffman <hoffmaao@uw.edu>
#
# This file is part of the two-stage glacier model.
#
# the two-stage glacier model is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The full text of the license can be found in the file LICENSE in the
# icepack source directory or at <http://www.gnu.org/licenses/>.

r"""Physical constants
This module constains physical constants used throughout the library such
as the acceleration due to gravity, the universal gas constant, the density
of ice and water, etc.
"""

#: number of seconds in a year, for unit conversions
year = 365.25 * 24 * 60 * 60

#: exponent in the nonlinear consitutive law for ice
glen_flow_law = 3.0

#: acceleration due to gravity (m / yr^2)
gravity = 9.81

#: density of ice
ice_density = 917

#: density of seawater
water_density = 1024