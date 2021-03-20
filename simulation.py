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
from sympy import integrate
import model
import forcing
from constants import (gravity as g,
	year as year,
	glen_flow_law as n,
	ice_denisty as ρ_I,
	water_density as ρ_W
	)


def equilibrium(L,H,b0,m0,S0=Sbar,Ωbar=Ωbar,rho_w=ρ_W,rho_i=rho_I,α=α,β=β,nts=n):
    S0 = S0/year
    Larray=np.zeros(nts)
    harray=np.zeros(nts)
    for t in range(nts):
        # update fluxes
        b = model.bed(b0,m,L)
        Hg = model
        Q = ν*h0**α/(L0**n)
        Qg = Ωbar*Hg**β
        # update geometry
        dh = dhdt(S0,h0,L0,Q,Qg)*dt*year
        dxg = dLdt(Q,Qg,hg)*dt*year
        h0 = h0 + dh
        L0 = L0 + dxg
        harray[t]=h0
        Larray[t]=L0
    return harray,Larray


    h=h0[-1]
xg=L0[-1]
bt=b0
mt=bx

def simulation(n,H0,L0,m0,b0,Ω,S,α=,β=,ν=):
	L = np.zeros(n)
	H = np.zeros(n)
	b = np.zeros(n)
	m = np.zeros(n)
	Q = np.zeros(n)
	Qg = np.zeros(n)

	Htmp=H0.copy()
	Ltmp=L0.copy()
	mtmp=m0.copy()





	for t in range(n):

    	dm=model.dmdt(htmp,Ltmp,L0,mtmp)

    	mtmp=mtmp+dm*dt
    	b= model.bed(b0,mtmp,Ltmp)
    	Hg= model.grounding_thickness(b)
    	Q[t]=model.interior_flux(L,H,ν,α,n)
    	Qg[t]=model.grounding_flux(Ω,b0,m,L,β)

    	dh=model.dhdt(S[t],Htmp,Ltmp,Q[t],Qg[t])*dt*year
    	dL=model.dLdt(Q[t],Qg[t],Hg)*dt*year
    	Htmp = Htmp + dh
    	Ltmp = Ltmp + dL
    	m[t] = mtmp
    	H[t] = Htmp
    	L[t] = Ltmp

    return L,H,m,Q,Qg

