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
	)


def noise(N):
	r"""
	return array of white noise
	"""
	
	noise=np.random.randn(int(N))
	noise=noise/np.std(noise)
	return noise

def anomaly(start, end):
	r"""
	return start and end anomaly
	"""

	return np.linspace(0.0,1.0,end-start) 

def smb(N, Sbar, Sσ, δS, start, end):
	r"""
	return surface mass balance array
	"""

	δ = anomaly(start,end)           
	return δS*Sbar*np.concatenate([np.zeros(start), δ, np.ones(N-end)]) + Sbar + Sσ*Sbar*noise(N)

def Ωbar(C, A, θ, m):
	r"""
	return the melt rate.
	"""

	λ=ρ_W/ρ_I
	return (A*(ρ_I*g)**(n+1)*(θ*(1-λ**-1))**n*(4**n*C)**-1)**(1/(m+1))

def Ω(N, Ωbar, Ωσ, δO, start, end):
	r"""
	
	"""

	δ=anomaly(start,end)
	return Ωbar*(δO*np.concatenate([np.zeros(start), δ, np.ones(N-end)]) + 1.0 + Ωσ*noise(N))

