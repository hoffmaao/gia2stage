# Copyright (C) 2020-2021 by Andrew Hoffman <hoffmaao@uw.edu>
#
# This file is part of the GIA two-stage code.
#
# the two-stage glacier model is free software: you can redistribute it 
# and/or modify it under the terms of the GNU GenerAL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) PuBL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)ic License as 
# puBL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)ished by the Free Software Foundation, either version 3 of the License,
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
	
	return -(Hbar**(-1+α)*Lbar**(-γ)*α*ν)/((b0+Lbar*mbar)*γ)


def AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω):
    r"""
    """
    term1 = -(-Hbar**α*Lbar**(-1-γ)*γ*ν + mbar*β*γ*Ω*(-(b0 + Lbar*mbar)*γ)**(-1+β))/((b0 + Lbar*mbar)*λ)
    term2 = (mbar*(Hbar**α*Lbar**(-γ)*ν + Ω*(b0+Lbar*mbar)*λ)**β)/((b0+Lbar*mbar)**2*λ)
    return term1+term2

def Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω):
    r"""
    """
    term1 = -(Lbar*β*(-(b0+Lbar*mbar)*λ)**(-1+β)*Ω)/(b0+Lbar*mbar)
    term2 = (Lbar*(Hbar**α)*Lbar**(-γ)*ν-(-(b0+Lbar*mbar)*λ)**β*Ω)/((b0+Lbar*mbar)**2*λ)

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


def Hprime(H,L,bx,Hbar,Lbar,mbar,hgbar,Qgbar,α,γ,ν,β,λ,b0,Ω,S,κ,τ,dt):
	r"""
	"""
	
	Qgi=Ω*Qgbar
	ah=AH(Hbar,Lbar,mbar,α,γ,ν,b0)
	al=AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	abx=Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	bh=BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	bl=BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	bbx=Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	ch=CH(Lbar,κ,τ)
	cl=CL(Lbar,Hbar,mbar,κ,τ,b0)
	cbx=Cbx(Lbar,κ,τ)



	F = (dt*((-1 + dt*(bl + cbx + bbx*cl*dt - bl*cbx(Lbar,κ,τ)*dt))*Hbar +
		dt*(AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0)*dt - AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ)*dt)*Lbar)*Qgi +
	hgbar*((-1 + dt*(BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + Cbx(Lbar,κ,τ) + Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0)*dt - BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ)*dt))*H*Lbar + 
		dt*(-Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Lbar*(CL(Lbar,Hbar,mbar,κ,τ,b0)*dt*L + bx - BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt*bx) + 
			AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Lbar*((-1 + Cbx(Lbar,κ,τ)*dt)*L - Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt*bx) + (1 - 
				dt*(Cbx(Lbar,κ,τ) + Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0)*dt) + BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt*(-1 + Cbx(Lbar,κ,τ)*dt))*(Qgi - 
				Lbar*S))))/((-1 + dt*(AH(Hbar,Lbar,mbar,α,γ,ν,b0) + BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + 
	Cbx(Lbar,κ,τ) + (AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ) + Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0) - BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ) - 
		AH(Hbar,Lbar,mbar,α,γ,ν,b0)*(BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + Cbx(Lbar,κ,τ)))*dt + (-Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ) + AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ) + Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0) - 
		AH(Hbar,Lbar,mbar,α,γ,ν,b0)*Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0) - AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ) + AH(Hbar,Lbar,mbar,α,γ,ν,b0)*BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ))*dt**2))*hgbar*Lbar)

	return F

def Lprime(H,L,bx,Hbar,Lbar,mbar,hgbar,Qgbar,α,γ,ν,β,λ,b0,Ω,S,κ,τ,dt):
	r"""
	"""

	Qgi=Ω*Qgbar
	AH=AH(Hbar,Lbar,mbar,α,γ,ν,b0)
	AL=AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	Abx=Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	BH=BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	BL=BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	Bbx=Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	CH=CH(Lbar,κ,τ)
	CL=CL(Lbar,Hbar,mbar,κ,τ,b0)
	Cbx=Cbx(Lbar,κ,τ)


	G = ((1 - dt*(Cbx(Lbar,κ,τ) + Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ)*dt) + AH(Hbar,Lbar,mbar,α,γ,ν,b0)*dt*(-1 + Cbx(Lbar,κ,τ)*dt))*Lbar*(-hgbar*L +
		dt*Qgi) - Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt*((1 - AH(Hbar,Lbar,mbar,α,γ,ν,b0)*dt)*hgbar*Lbar*bx + 
		CH(Lbar,κ,τ)*dt*(hgbar*H*Lbar + dt*Hbar*Qgi - dt*hgbar*Qgi +
			dt*hgbar*Lbar*S)) + BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt*((-1 + Cbx(Lbar,κ,τ)*dt)*hgbar*H*Lbar + dt*(-1 + Cbx(Lbar,κ,τ)*dt)*Hbar*Qgi +
		dt*hgbar*(-Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Lbar*bx - (-1 + Cbx(Lbar,κ,τ)*dt)*(Qgi - Lbar*S))))/((-1 + 
		dt*(AH(Hbar,Lbar,mbar,α,γ,ν,b0) + BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) +
			Cbx(Lbar,κ,τ) + (AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ) + Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0) - BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ) - 
				AH(Hbar,Lbar,mbar,α,γ,ν,b0)*(BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + Cbx(Lbar,κ,τ)))*dt + (-Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ) + AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ) + Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0) - 
				AH(Hbar,Lbar,mbar,α,γ,ν,b0)*Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0) - AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ) + AH(Hbar,Lbar,mbar,α,γ,ν,b0)*BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ))*dt**2))*hgbar*Lbar)

	return G

def bxprime(H,L,bx,Hbar,Lbar,mbar,hgbar,Qgbar,α,γ,ν,β,λ,b0,Ω,S,κ,τ,dt):
	r"""
	"""

	Qgi=Ω*Qgbar
	AH=AH(Hbar,Lbar,mbar,α,γ,ν,b0)
	AL=AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	Abx=Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	BH=BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	BL=BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	Bbx=Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)
	CH=CH(Lbar,κ,τ)
	CL=CL(Lbar,Hbar,mbar,κ,τ,b0)
	Cbx=Cbx(Lbar,κ,τ)


	J = ((-1 + AH(Hbar,Lbar,mbar,α,γ,ν,b0)*dt)*Lbar*((1 - BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt)*hgbar*bx + 
      CL(Lbar,Hbar,mbar,κ,τ,b0)*dt*(hgbar*L - dt*Qgi)) - 
   BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt**2*(-AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*hgbar*Lbar*bx + CL(Lbar,Hbar,mbar,κ,τ,b0)*dt*Hbar*Qgi + 
      CL(Lbar,Hbar,mbar,κ,τ,b0)*hgbar*(H*Lbar - dt*Qgi + dt*Lbar*S)) + 
   CH(Lbar,κ,τ)*dt*(dt*((-1 + BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt)*Hbar + AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt*Lbar)*Qgi + 
      hgbar*((-1 + BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt)*H*Lbar - AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt*Lbar*L - 
         dt*(-1 + BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*dt)*(Qgi - Lbar*S))))/((-1 + 
     dt*(AH(Hbar,Lbar,mbar,α,γ,ν,b0) + BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + 
        Cbx(Lbar,κ,τ) + (AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ) + Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0) - BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ) - 
           AH(Hbar,Lbar,mbar,α,γ,ν,b0)*(BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω) + Cbx(Lbar,κ,τ)))*dt + (-Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ) + AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CH(Lbar,κ,τ) + Abx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0) - 
           AH(Hbar,Lbar,mbar,α,γ,ν,b0)*Bbx(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*CL(Lbar,Hbar,mbar,κ,τ,b0) - AL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*BH(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ) + AH(Hbar,Lbar,mbar,α,γ,ν,b0)*BL(Hbar,Lbar,mbar,α,γ,ν,β,λ,b0,Ω)*Cbx(Lbar,κ,τ))*dt**2))*hgbar*Lbar)

	return J


