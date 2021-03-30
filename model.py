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
import sympy
from constants import (gravity as g,
	glen_flow_law as n,
	ice_density as ρ_I,
	water_density as ρ_W,
	bed_density as ρ_b,
	bedrock_rigidity as D
	)

def bed(b0, m, L):
	r"""
	return bed elevation at grounding line
	"""
	return b0+m*L

def grounding_thickness(b):
	r"""
	return grounding line thickness
	"""

	return -ρ_W/ρ_I*b

def interior_flux(L, H, ν, α, n):
	r"""
	the interior flux of the ice sheet.
	"""

	return ν*H**α/(L**n)

def grounding_flux(Ω, b0, m, L, β):
	r"""
	return grounding zone flux
	"""

	return Ω*grounding_thickness(bed(b0,m,L))**β

def dhdt(S, H, L, b, Q, Qg):
	r"""
	return change in thickness
	"""
	dhdt = S - Qg/L - (H/(grounding_thickness(b)*L))*(Q-Qg)
	return dhdt

def dLdt(Q, Qg, hg):
    r"""
    description of the change in glacier length
    """
    dLdt = (Q-Qg)/hg
    return dLdt

def dbdt(H, H0, b, b0, S):
	r"""
	"""
	σ=ρ_I * g * (H - H0)
	return -σ/D


def stress_profile(H, L, m, b, x):
	r"""
	returns the stress profile of 
	"""
	x = sympy.symbols('x', real=True, positive=True)
	Hg = grounding_thickness(bed(b,m,L))
	σ_I = ρ_I * g * ((H-Hg) * sympy.sqrt((L-x)/L)+ Hg)
	σ_W = ρ_W * g * -(m * x + b)
	return sympy.Piecewise((σ_I,x<=L),(σ_W,x>L))



def dmdt(H, L, m, H0, L0, m0, b0, D):
	r"""

	"""
	x = sympy.symbols('x', real=True, positive=True)
	Lmax=np.max((L0,L))
	Vf = sympy.integrate(stress_profile(H,L,m,b0,x) - stress_profile(H0,L0,m0,b0,x),(x,0,Lmax))
	return -Vf/Lmax/D