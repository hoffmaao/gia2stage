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

import tqdm
import copy
import numpy as np
import model
import linear
import forcing
from constants import (gravity as g,
	ice_density as ρ_I,
	water_density as ρ_W,
	bedrock_density as ρ_s,
	glen_flow_law as n,
	glen_coefficient as A,
	friction_coefficient as C,
	theta as θ,
	year as year
	)


def m(n):
	return 1/n

def γ(n):
	return n

def α(n):
	return 2*n+1

def β(bedn,n):
	return (m(bedn)+n+3)/(m(bedn)+1)

def ν(n, C):
	return (ρ_I*g/C)**n
def λ(ρ_W,ρ_I):
    return ρ_I/ρ_W



def initialize_forcing(Sbar, Sσ, Obar, Oσ, N, δS=0.0, δO=0.0, startS=0, endS=0, startO=0, endO=0):
	r"""
	return ocean melt and surface mass balance forcing
	"""

	return forcing.Ω(N,Obar,Oσ,δO,startO,endO), forcing.smb(N,Sbar,Sσ,δS,startS,endS)

def equilibrium(H, L, m0, b0, Ωbar, Sbar,C=C,bedn=n,alphan=n,n=n, N=1000000, dt=1.0,):
	r"""
	return equilibrium thickness and length
	H : float
         initial interior ice thickness
    L : float
         initial glacier length
    m0 : float
         initial bedrock slope
    b0 : float
         interior bedrock elevation
    Ωbar : float
        average discharge coefficient
    Sbar : float
        average surface mass balace
    N : float
        number of timesteps
	"""

	for t in tqdm.trange(N):
		b = model.bed(b0,m0,L)
		Hg = model.grounding_thickness(b)
		Q = model.interior_flux(L,H,ν(bedn,C),α(bedn),bedn)
		Qg = model.grounding_flux(Ωbar,b0,m0,L,β(bedn,n))
		dh = model.dhdt(Sbar,H,L,b,Q,Qg)*dt*year
		dL = model.dLdt(Q,Qg,Hg)*dt*year
		H = H + dh
		L = L + dL
	return H, L

def simulation(H0, L0, m0, b0, Ω, S, N, bedn=n,alphan=n,n=n, C=C, dt=1.0):
    r"""
    return timeseries of ice thickness, glacier length, 
    bedrock slope and fluxes.

    H0 : float
         initial interior ice thickness
    L0 : float
         initial glacier length
    m0 : float
         initial bedrock slope
    b0 : float
         interior bedrock elevation
    Ω : float
        discharge coefficient
    S : float
        surface mass balace
    N : float
        number of timesteps

    """

    L = np.zeros(N)
    H = np.zeros(N)
    M = np.zeros(N)
    Q = np.zeros(N)
    Qg = np.zeros(N)
    Htmp=H0
    Ltmp=L0
    σ=0.0
    σp=0.0

    #btmp=b0
    for t in tqdm.trange(N):
    	H[t] = Htmp
    	L[t] = Ltmp
    	b = model.bed(b0,m0,Ltmp)
    	Hg= model.grounding_thickness(b)
    	Q[t]=model.interior_flux(Ltmp,Htmp,ν(bedn,C),α(bedn,),n)
    	Qg[t]=model.grounding_flux(Ω[t],b0,m0,Ltmp,β(bedn,n))
    	dH=model.dhdt(S[t],Htmp,Ltmp,b,Q[t],Qg[t])*dt*year
    	dL=model.dLdt(Q[t],Qg[t],Hg)*dt*year
    	Htmp = Htmp + dH
    	Ltmp = Ltmp + dL

    return H,L,Q,Qg




def elastic_simulation(H0, L0, m0, b0, Ω, S, N, E, γ=1.0, dt=1.0):
    r"""
    return timeseries of ice thickness, glacier length, 
    bedrock slope and fluxes.

    H0 : float
         initial interior ice thickness
    L0 : float
         initial glacier length
    m0 : float
         initial bedrock slope
    b0 : float
         interior bedrock elevation
    Ω : float
        discharge coefficient
    S : float
        surface mass balace
    N : float
        number of timesteps
    E : float
        elasticity
    γ : float
        aspect ratio
    """

    L = np.zeros(N)
    H = np.zeros(N)
    M = np.zeros(N)
    Q = np.zeros(N)
    Qg = np.zeros(N)
    Htmp=H0
    Ltmp=L0
    mtmp=m0
    #btmp=b0
    for t in tqdm.trange(N):
    	M[t] = mtmp
    	H[t] = Htmp
    	L[t] = Ltmp
    	b = model.bed(b0,mtmp,Ltmp)
    	Hg= model.grounding_thickness(b)
    	Q[t]=model.interior_flux(Ltmp,Htmp,ν(n,C),α(n),n)
    	Qg[t]=model.grounding_flux(Ω[t],b0,mtmp,Ltmp,β(n))
    	dH=model.dhdt(S[t],Htmp,Ltmp,b,Q[t],Qg[t])*dt*year
    	dL=model.dLdt(Q[t],Qg[t],Hg)*dt*year
    	σ=model.stress(Htmp,Ltmp,mtmp, H0, L0, m0, b0)
    	dm=model.elastic(σ,E,γ)
    	mtmp =mtmp + dm
    	Htmp = Htmp + dH
    	Ltmp = Ltmp + dL

    return H,L,M,Q,Qg


def viscoelastic_simulation(H0, L0, m0, b0, Ω, S, N, E, η, γ=1.0, dt=1.0):
    r"""
    return timeseries of ice thickness, glacier length, 
    bedrock slope and fluxes.

    H0 : float
         initial interior ice thickness
    L0 : float
         initial glacier length
    m0 : float
         initial bedrock slope
    b0 : float
         interior bedrock elevation
    Ω : float
        discharge coefficient
    S : float
        surface mass balace
    N : float
        number of timesteps
    E : float
        elasticity
    γ : float
        aspect ratio
    """

    L = np.zeros(N)
    H = np.zeros(N)
    M = np.zeros(N)
    Q = np.zeros(N)
    Qg = np.zeros(N)
    Htmp=H0
    Ltmp=L0
    mtmp=m0
    σ=0.0
    σp=0.0

    #btmp=b0
    for i in tqdm.trange(N):
    	M[t] = mtmp
    	H[t] = Htmp
    	L[t] = Ltmp
    	b = model.bed(b0,mtmp,Ltmp)
    	Hg= model.grounding_thickness(b)
    	Q[t]=model.interior_flux(Ltmp,Htmp,ν(n,C),α(n),n)
    	Qg[t]=model.grounding_flux(Ω[t],b0,mtmp,Ltmp,β(n))
    	dH=model.dhdt(S[t],Htmp,Ltmp,b,Q[t],Qg[t])*dt*year
    	dL=model.dLdt(Q[t],Qg[t],Hg)*dt*year
    	σp=copy.deepcopy(σ)
    	σ=model.stress(Htmp,Ltmp,mtmp,H0,L0,m0,b0)
    	dm=model.viscoelastic(σ,σp,E,η,γ,dt*year)*dt*year
    	mtmp=mtmp + dm
    	Htmp = Htmp + dH
    	Ltmp = Ltmp + dL

    return H,L,M,Q,Qg

def Oerlemans_simulation(H0, L0, m0, b0, Ω, S, N, τ, dt=1.0):
    r"""
    return timeseries of ice thickness, glacier length, 
    bedrock slope and fluxes.

    H0 : float
         initial interior ice thickness
    L0 : float
         initial glacier length
    m0 : float
         initial bedrock slope
    b0 : float
         interior bedrock elevation
    Ω : float
        discharge coefficient
    S : float
        surface mass balace
    τ : float
        timescale of rebound
    E : float
        elasticity
    γ : float
        aspect ratio
    """

    L = np.zeros(N)
    H = np.zeros(N)
    M = np.zeros(N)
    Q = np.zeros(N)
    Qg = np.zeros(N)
    Htmp=H0
    Ltmp=L0
    mtmp=m0
    σ=0.0
    σp=0.0

    #btmp=b0
    for t in tqdm.trange(N):
    	M[t] = mtmp
    	H[t] = Htmp
    	L[t] = Ltmp
    	b = model.bed(b0,mtmp,Ltmp)
    	Hg = model.grounding_thickness(b)
    	Q[t] = model.interior_flux(Ltmp,Htmp,ν(n,C),α(n),n)
    	Qg[t] = model.grounding_flux(Ω[t],b0,mtmp,Ltmp,β(n))
    	dH = model.dhdt(S[t],Htmp,Ltmp,b,Q[t],Qg[t])*dt*year
    	dL = model.dLdt(Q[t],Qg[t],Hg)*dt*year

    	w = model.deflection_angle(Htmp,Ltmp,mtmp,H0,L0,m0,b0)
    	dm = model.dmdt(w,mtmp,m0,τ)*dt
    	mtmp = mtmp + dm
    	Htmp = Htmp + dH
    	Ltmp = Ltmp + dL

    return H,L,M,Q,Qg

def linear_simulation(Hbar, Lbar, mbar, hgbar, Qgbar, b0, Ω, S, τ, κ, N, γ=1.0, dt=1.0):
    r"""
    return timeseries of ice thickness, glacier length, 
    bedrock slope and fluxes.

    Hbar : float
        initial interior ice thickness
    Lbar : float
        initial glacier length
    mbar : float
        initial bedrock slope
    hgbar : float
        average grounding zone thickness
    Qgbar : float
        average grounding zone flux
    b0 : float
        interior bedrock elevation
    Ω : float
        discharge coefficient
    S : float
        surface mass balace
    τ : float
        timescale of rebound
    γ : float
        aspect ratio
    """

    L = np.zeros(N)
    H = np.zeros(N)
    M = np.zeros(N)
    dt=dt*year


    for t in tqdm.trange(1,N):
        L[t]=linear.Lprime(H[t-1],L[t-1],M[t-1],Hbar,Lbar,mbar,hgbar,Qgbar,α(n),γ,ν(n,C),β(n),λ(ρ_W,ρ_I),b0,Ω[t],S[t],κ,τ,dt)
        H[t]=linear.Hprime(H[t-1],L[t-1],M[t-1],Hbar,Lbar,mbar,hgbar,Qgbar,α(n),γ,ν(n,C),β(n),λ(ρ_W,ρ_I),b0,Ω[t],S[t],κ,τ,dt)
        M[t]=linear.bxprime(H[t-1],L[t-1],M[t-1],Hbar,Lbar,mbar,hgbar,Qgbar,α(n),γ,ν(n,C),β(n),λ(ρ_W,ρ_I),b0,Ω[t],S[t],κ,τ,dt)
    return L,H,M




