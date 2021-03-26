import numpy as np
import forcing
import simulation
import plot
from constants import (gravity as g,
	glen_flow_law as n,
	ice_density as ρ_I,
	water_density as ρ_W,
	glen_flow_law as n,
	glen_coefficient as A,
	friction_coefficient as C,
	theta as θ
	)
import matplotlib.pyplot as plt



# Simulation parameters
N = 100000
Sbar = .5
Sσ = Sbar/1.
Obar = forcing.Ωbar(C,A,θ,simulation.m(n))
Oσ = 0.0
dt = 1.0

# Initialize and forcing
Ω,smb=simulation.initialize_forcing(Sbar, Sσ, Obar, Oσ, N)

fig,ax,plots=plot.makefig01(1)
plot.figure01(fig,ax,smb,N,dt,label=r'smb.($ma^{-1}$)')
plt.savefig('figures/smb.png')
plt.close()

fig,ax,plots=plot.makefig01(1)
plot.figure01(fig,ax,Ω,N,dt,label=r'Ω')
plt.savefig('figures/omega.png')

b0=-100
m0=-2e-3

H0,L0=simulation.equilibrium(1413,184000,b0,m0,Sbar,Obar)
print(H0)
print(L0)

H,L,m,Q,Qg=simulation.simulation(H0,L0,m0,b0,Ω,smb,N)

fig,ax1,ax2,plots=plot.makefig02(1)
plot.figure02(fig,ax1,ax2,H,L/1000,N,dt)
plt.savefig('figures/simulation.png')
plt.close()