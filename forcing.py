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
	ice_denisty as ρ_I,
	water_density as ρ_W
	)


def noiseS(nts):
	noiseS=np.random.randn(int(nts))
	noiseS=noiseS/np.std(noiseS)
	return noiseS

def noiseO(nts):
	noiseO = np.random.randn(int(nts))
	noiseO = noiseO/np.std(noiseO)
	return noiseO

def smb(nts,year=1.0,Sbar = 0.5,sigS=.5,deltaS=0.0,start=0.0,end=0.0):
	anom = np.linspace(0.0,1.0,end-start)           
	S = deltaS*Sbar*np.concatenate([np.zeros(startyrS), anom, np.ones(nts-endyrS)]) + Sbar + sigS*Sbar*noiseS(nts)
	S = S/year
	return S

def omega_bar(nts,theta=0.7,C=7.624e6,A_glen=4.22e-25):
	m=1/n
	alpha=2*n+1
	gamma=n
	beta=(m+n+3)/(m+1)            ### buttressing parameter ###
	lam=ρ_W/ρ_I
	omega_bar=(A_glen*(ρ_I*g)**(n+1)*(theta*(1-lam**-1))**n*(4**n*C)**-1)**(1/(m+1))
	return omega_bar

def omega(nts,year=1.0,sigO=0.0,deltaO=0.0,start=0.0,end=0):
	anom=np.linespace(0,1,end-start)
	return omega=deltaO*omega_bar(nts)*np.concatenate([np.zeros(startyrO), anom, np.ones(nts-endyrO)]) + omega_bar(nts) + sigO*omega_bar*noiseO(nts)

