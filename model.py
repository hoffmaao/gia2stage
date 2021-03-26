# Copyright (C) 2020-2021 by Andrew Hoffman <hoffmaao@uw.edu>
#
# This file is part of the two-stage code.
#
# the two-stage glacier model is free software: you can redistribute it 
# and/or modify it under the terms of the GNU General Public License as 
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# The full text of the license can be found in the file LICENSE in the
# gravity-inversion source directory or at <http://www.gnu.org/licenses/>.

import numpy as np
from constants import (gravity as g,
	glen_flow_law as n,
	ice_density as ρ_I,
	water_density as ρ_W,
	bedrock_rigidity as D
	)

def bed(b0,m,L):
	r"""
	return bed elevation at grounding line
	"""
	return b0+m*L

def grounding_thickness(b):
	r"""
	return grounding line thickness
	"""

	return -ρ_W/ρ_I*b

def interior_flux(L,H,ν,α,n):
	r"""
	
	"""

	return ν*H**α/(L**n)

def grounding_flux(Ω,b0,m,L,β):
	r"""
	return grounding zone flux
	"""

	return Ω*grounding_thickness(bed(b0,m,L))**β

def dhdt(S,H,L,b,Q,Qg):
	r"""
	return change in thickness
	"""
	dhdt = S - Qg/L - (H/(grounding_thickness(b)*L))*(Q-Qg)
	return dhdt

def dLdt(Q,Qg,hg):
    r"""
    description of the change in glacier length
    """
    dLdt = (Q-Qg)/hg
    return dLdt

def dmdt(H,L,H0,L0):
	r"""

	"""
	
	Vf=(L0-L)*g*ρ_I
	dmdt= Vf/D*L0
	return dmdt