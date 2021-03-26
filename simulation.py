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
import model
import forcing
from constants import (gravity as g,
	glen_flow_law as n,
	ice_density as ρ_I,
	water_density as ρ_W,
	glen_flow_law as n,
	glen_coefficient as A,
	friction_coefficient as C,
	theta as θ,
	year as year
	)


def m(n):
	return 1/n

def γ(n):
	n

def α(n):
	return 2*n+1

def β(n):
	return (m(n)+n+3)/(m(n)+1)

def ν(n, C):
	return (ρ_I*g/C)**n 



def initialize_forcing(Sbar, Sσ, Obar, Oσ, N, δS=0.0, δO=0.0, startS=0, endS=0, startO=0, endO=0):
	r"""
	return ocean melt and surface mass balance forcing
	"""

	return forcing.Ω(N,Obar,Oσ,δO,startO,endO), forcing.smb(N,Sbar,Sσ,δS,startS,endS)

def equilibrium(H, L, b0, m0, Sbar, Ωbar, N=100000, dt=1.0):
	r"""
	return equilibrium thickness and length
	"""

	Sbar=Sbar/year
	for t in range(N):
		b = model.bed(b0,m0,L)
		Hg = model.grounding_thickness(b)
		Q = model.interior_flux(L,H,ν(n,C),α(n),n)
		Qg = model.grounding_flux(Ωbar,b0,m0,L,β(n))
		dh = model.dhdt(Sbar,H,L,b,Q,Qg)*dt*year
		dL = model.dLdt(Q,Qg,Hg)*dt*year
		H = H + dh
		L = L + dL
	return H, L

def simulation(H0, L0, m0, b0, Ω, S, N, dt=1.0):
    r"""
    return timeseries of thickness, glacier, length, slope and fluxes.
    """

    L = np.zeros(N)
    H = np.zeros(N)
    b = np.zeros(N)
    bx = np.zeros(N)
    Q = np.zeros(N)
    Qg = np.zeros(N)
    Htmp=H0
    Ltmp=L0
    mtmp=m0
    for t in range(N):
    	bx[t] = mtmp
    	H[t] = Htmp
    	L[t] = Ltmp
    	b= model.bed(b0,mtmp,Ltmp)
    	Hg= model.grounding_thickness(b)
    	Q[t]=model.interior_flux(Ltmp,Htmp,ν(n,C),α(n),n)
    	Qg[t]=model.grounding_flux(Ω[t],b0,mtmp,Ltmp,β(n))
    	dH=model.dhdt(S[t]/year,Htmp,Ltmp,b,Q[t],Qg[t])*dt*year
    	dL=model.dLdt(Q[t],Qg[t],Hg)*dt*year
    	dm=model.dmdt(Htmp,Ltmp,H0,L0)*dt*year
    	mtmp=mtmp + dm
    	Htmp = Htmp + dH
    	Ltmp = Ltmp + dL

    return H,L,bx,Q,Qg

