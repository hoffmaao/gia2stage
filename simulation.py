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
# two-stage glacier model source directory or at <http://www.gnu.org/licenses/>.

import numpy as np
from sympy import integrate
import model
import forcing
from constants import (gravity as g,
	year as year,
	glen_flow_law as n,
	ice_denisty as ρ_I,
	water_density as ρ_W
	)


def m(n):
	return 1/n

def α(n):
	return 2*n+1

def β(n):
	return (m(n)+n+3)/(m(n)+1)

def ν(n,C):
	return (ρ_I*g/C)**n 



def initialize_forcing(N, Sbar, Sσ, Obar, Oσ, δS=0.0, δO=0.0, startS=0.0, endS=0.0, startO=0.0, endO=0.0):
	r"""
	return ocean melt and surface mass balance forcing
	"""

	return forcing.Ω(N,Obar,Oσ,δO,startO,endO), forcing.smb(N,Sbar,Sσ,δS,startS,endS)

def equilibrium(L,H,b0,m0,C,S,Ωbar,dt=1.0,N=100000):
	r"""
	return equilibrium thickness and length
	"""

    for t in range(N):
        b = model.bed(b0,m,L)
        Hg = model.grounding_thickness(b)
        Q = model.interior_flux(L,H,ν(n,C),α(n),n)
        Qg = model.grounding_flux(Ω,b0,m(n),L,β(n))
        dh = dhdt(S,H,L,Q,Qg)*dt*year
        dL = dLdt(Q,Qg,Hg)*dt*year
        h = h + dh
        L = L + dL

    return H, L

def simulation(N, H0, L0, m0, b0, Ω, S):
	r"""
	return timeseries of thickness, glacier, length, slope and fluxes.
	"""

	L = np.zeros(N)
	H = np.zeros(N)
	b = np.zeros(N)
	m = np.zeros(N)
	Q = np.zeros(N)
	Qg = np.zeros(N)

	Htmp=H0.copy()
	Ltmp=L0.copy()
	mtmp=m0.copy()



	for t in range(N):
		m[t] = mtmp
    	H[t] = Htmp
    	L[t] = Ltmp

    	b= model.bed(b0,mtmp,Ltmp)
    	Hg= model.grounding_thickness(b)
    	Q[t]=model.interior_flux(L,H,ν(n,C),α(n),n)
    	Qg[t]=model.grounding_flux(Ω,b0,m(n),L,β(n))

    	dh=model.dhdt(S[t],Htmp,Ltmp,Q[t],Qg[t])*dt
    	dL=model.dLdt(Q[t],Qg[t],Hg)*dt
    	dm=model.dmdt(htmp,Ltmp,L0,mtmp)*dt
    	mtmp=mtmp + dm
    	Htmp = Htmp + dH
    	Ltmp = Ltmp + dL

    return L,H,m,Q,Qg

