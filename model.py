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
	bedrock_density as ρ_s,
	bedrock_rigidity as D
	)

def bed(b, m, L):
	r"""
	calculate bed elevation at grounding zone

	b : float
	    interior bed elevation
	m : float
	    bedrock slope
	L : float
	    glacier length
	"""

	return b+m*L

def grounding_thickness(b):
	r"""
	calculate ice thickness at grounding line

	b : float
	    bedrock elevtion at grounding zone
	"""

	return -ρ_W/ρ_I*b

def interior_flux(L, H, ν, α, n):
	r"""
	calculate interior flux of the ice sheet.

	L : float
	    glacier length
	H : float
	    interior thickness
	ν : float
	    interior coefficient
	α : float
	    interior exponent
	n : float
	    glen exponent
	"""

	return ν*H**α/(L**n)

def grounding_flux(Ω, b0, m, L, β):
	r"""
	calculate grounding zone flux

	Ω : float
	    ocean flux coefficient
	b0 : float
	     interior bed elevation
	m : float
	    bedrock slope
	L : float
	    glacier length
	β : float
	    grounding zone exponent
	"""

	return Ω*grounding_thickness(bed(b0,m,L))**β

def dhdt(S, H, L, b, Q, Qg):
	r"""
	calculate change in thickness

	S : float
	    surface mass balance.(m/a)
	H : float
	    interior ice thickness.(m)
	L : float
	    glacier length.(m)
	b : float
	    interior bedrock elevation.(m)
	Q : float
	    interior flux
	Qg : float
	     grounding zone flux
	"""

	return S - Qg/L - (H/(grounding_thickness(b)*L))*(Q-Qg)

def dLdt(Q, Qg, hg):
    r"""
    calculate change in glacier length

    Q : float
        interior flux
    Qg : float
         grounding zone flux
    hg : float
         grounding zone ice thickness
    """

    return (Q-Qg)/hg

def stress_profile(H, L, m, b, x):
	r"""
	calculate the stress profile assuming Weertman equilibrium profile.

	H : float 
	    glaicer thickness
	L : float
	    glacier length
	m : float
	    bedrock slope
	b : float
	    interior bed elevation
	x : sympy dependent variable
	    distance
	"""

	Hg = grounding_thickness(bed(b,m,L))
	σ_I = ρ_I * g * ((H-Hg) * sympy.sqrt((L-x)/L)+ Hg)
	σ_W = ρ_W * g * -(m * x + b)

	return sympy.Piecewise((σ_I,x<=L),(σ_W,x>L))

def stress(H, L, m, H0, L0, m0, b0):
	r"""
	calculate the average stress difference along the profile

	H : float
	    interior ice thickness
	L : float
	    glacier length
	m : float
	    bedrock slope
	H0 : float
	     equalibrium interior thickness
	L0 : float
	     equalibrium glacier length
	m0 : float
	     equalibrium bedrock slope
	b0 : float
	     interior bedrock elevation
	"""

	x = sympy.symbols('x', real=True, positive=True)
	Lmax=np.max((L0,L))
	Vf = sympy.integrate(stress_profile(H,L,m,b0,x) - stress_profile(H0,L0,m0,b0,x),(x,0,Lmax))
	return Vf / Lmax


def deflection_angle(H, L, m, H0, L0, m0, b0):
	r"""
	calculate the equilibrium angle

	H : float
	    interior ice thickness
	L : float
	    glacier length
	m : float
	    bedrock slope
	H0 : float
	     equalibrium interior thickness
	L0 : float
	     equalibrium glacier length
	m0 : float
	     equalibrium bedrock slope
	b0 : float
	     interior bedrock elevation
	"""

	x = sympy.symbols('x', real=True, positive=True)
	Lmax=np.max((L0,L))
	Vf = sympy.integrate(-ρ_s * g * ((m * x + b0) - (m0 * x + b0)) +  stress_profile(H,L,m,b0,x) - stress_profile(H0,L0,m0,b0,x),(x,0,Lmax))
	return Vf / (ρ_s * g) / Lmax**2

def dmdt(w ,m , m0, τb):
	r"""
	calculate dmdt 

	w : float
	    current equalibrium slope
	m : float
	    bed slope
	m0 : float
	     initial equalibrium.
	τb : float
         timescale of deflection
	"""

	return -1/τb*(w+m-m0)


def elastic(σ,E,γ):
	r"""
	calculate the elastic deformation due to an applied stress.

	σ : float
	    stress disequalibrium
	E : float
	    bedrock elasticity
	"""

	return γ*σ/D

def viscoelastic(σ,σp,E,η,γ,dt):
	r"""
	calculate the viscoelastic deformation rate due to an applied stress.

	σ : float
	    current length average stress
	σp : float
	    previous length average stress
    E : float
        bedrock elasticity
    η : float
        bedrock viscosity
    dt : float
        time step
	"""

	return γ/η*(η/E*(σ-σp)/dt+σ)
