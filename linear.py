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
# gia2stage source directory or at <http://www.gnu.org/licenses/>.


import numpy as np
from constants import (gravity as g,
	glen_flow_law as n,
	ice_density as ρ_I,
	water_density as ρ_W,
	bedrock_density as ρ_s,
	bedrock_rigidity as D
	)



def AH(Hbar,Lbar,mbar,α,γ,ν,b0):
	r"""
	"""
	
	return -(Hbar**(-1+α)*Lbar**(-γ*α*ν)/((b0+Lbar*mbar)*γ))


def AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω):
	r"""
	"""
	term1 = -(-Hbar**α*Lbar**(-1-γ)*γ*ν + mbar*β*γ*Ω*(-(b0 + Lbar*mbar)*γ)**(-1+β))/((b0 + Lbar*mbar)*λ)
	term2 = (mbar*(Hbar**α*Lbar**(-γ)*ν + Ω*(b0+Lbar*mbar)*λ)**β)/((b0+Lbar*mbar)**2*λ)

 	return term1+term2

def Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω):
	r"""
	"""
    term1 = -((Lbar*β*(-b0+Lbar*mbar)*λ)**(-1+β)*Ω)/(b0+Lbar*mbar)
    term2 = (Lbar*(Hbar**α)*Lbar**(-γ*ν)-(-(b0+Lbar*mbar)*λ)**β*Ω)/(b0+Lbar*mbar)**2*λ

	return term1+term2

def BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω):
	r"""
	"""
	term1 = Hbar**α*Lbar**(-1-γ)*α*ν/((b0+Lbar*mbar)*λ)
	term2 = Hbar**α*Lbar**(-γ*ν) -(-(b0+Lbar*mbar)*λ)**β*Ω/(Lbar*(b0+Lbar*mbar)*λ)

	return term1+term2

def BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω):
	r"""
	"""
	term1 = -((Lbar*β*(-(b0+Lbar*mbar)*λ)**(-1+β)*Ω))/Lbar
	term2 = ((-b0+Lbar*mbar)*λ)**β*Ω/(Lbar**2)
	term3 = Hbar*(-Hbar**α*Lbar**(-1-γ)*γ*ν+mbar*β*λ*(-(b0+Lbar*mbar)*λ)**(-1+β)*Ω)/(Lbar*(b0+Lbar*mbar)*λ)
	term4 = -(Hbar*mbar*(Hbar**α*Lbar**(-γ*ν)-(-(b0+Lbar*mbar)*λ)**β*Ω))/(Lbar*(b0+Lbar*mbar)**2*λ)
	term5 = -(Hbar*(Hbar**α*Lbar**(-γ*ν)-(-(b0+Lbar*mbar)*λ)**β*Ω))/(Lbar**2*(b0+Lbar*mbar)*λ)

	return term1 + term2 + term3 + term4 + term5

def Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω):
	r"""
	"""
	term1 = Hbar*β*(-(b0+Lbar*mbar)*λ)**(-1+β)*Ω/(b0+Lbar*mbar)
	term2 = β*λ*(-(b0+Lbar*mbar)*λ)**(-1+β)*Ω
	term3 = Hbar*(Hbar**α*Lbar**(-γ*ν)-(-(b0+Lbar*mbar)*λ)**β)*Ω/((b0+Lbar*mbar)**2*λ)

	return term1 + term2 + term3

def CH(Lbar,κ,τ):
	r"""
	"""

	return 2*Lbar*κ/(3*τ)

def CL(Lbar,Hbar,mbar,κ,τ,b0):
	r"""
	"""

	return (1/3*Lbar*mbar*κ - 1/3*(-b0+2*Hbar+Lbar*mbar)*κ)/τ

def Cbx(Lbar,κ,τ):
	r"""
	"""

	return (1 - (Lbar**2*κ/3))/τ


def Hprime(L,H,bx,Hp,Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω,dt):
	r"""
	"""

	F = AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*L+Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*bx+ 1/Lbar*(Hbar/hgbar-1)*Qg+S
	return 1/(1-dt*AH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω))*(dt*F+Hp)

def Lprime(L,H,bx,Lp,Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω,hgbar,Qg,dt):
	r"""
	"""
	
	G = BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*H + Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*bx -1/hg*Qg
	return 1/(1-dt*BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω))*(dt*G+Lp)

def bxprime(L,H,bx,bxp,Hbar,Lbar,mbar,κ,τ,hg,Qg,dt):
	r"""
	"""

	J = CH(Lbar,κ,τ)*H + CL(Lbar,Hbar,mbar,κ,τ,b0)*L
	return 1/(1-dt*Cbx(Lbar,κ,τ))*(dt*J+bxp)


