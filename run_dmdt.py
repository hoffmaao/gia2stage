import numpy as np
import forcing
import simulation
import plot

from constants import (gravity as g,
	glen_flow_law as n,
	ice_density as ρ_I,
	water_density as ρ_W,
	bedrock_density as ρ_s,
	glen_flow_law as n,
	glen_coefficient as A,
	friction_coefficient as C,
	theta as θ,
	year as year
	)
import matplotlib.pyplot as plt



# Simulation parameters
N = 1000000
Sbar = .5
Sσ = .5
Obar = forcing.Ωbar(C,A,θ,simulation.m(n))
Oσ = 0.0
dt = 1.0



# Initialize and forcing
Ω,smb=simulation.initialize_forcing(Sbar,Sσ,Obar,0.0,N)
Ω2,smb2=simulation.initialize_forcing(Sbar,0.0,Obar,1.0,N)

fig,ax,plots=plot.makefig01(1)
plot.figure01(fig,ax,smb,N,dt,label=r'smb.($ma^{-1}$)')
plt.savefig('figures/smb.png')
plt.close()

fig,ax,plots=plot.makefig01(1)
plot.figure01(fig,ax,Ω2,N,dt,label=r'Ω')
plt.savefig('figures/omega.png')

b0=-100
m0=-2e-3

H0,L0=simulation.equilibrium(1413,184000,m0,b0,Obar,Sbar/year)

print('starting simulation 1')
τ = np.array([3000])

from scipy import signal


fig1,ax1,ax2,plots=plot.makefig02(1)
fig2,ax,plots=plot.makefig01(1)
fig3,ax3,plots=plot.makefig01(1)



H,L,Q,Qg=simulation.simulation(H0,L0,m0,b0,Ω,smb,N)

x=np.linspace(0,N*dt,N)
color = 'black'

ax2.set_xlabel('year.(a)')
ax1.set_ylabel('thickness.(m)')
ax1.plot(x,H, color=color)
ax1.grid(True, which="both")

ax2.set_ylabel('length.(km)')  # we already handled the x-label with ax
ax2.plot(x,L/1000, color=color)
ax2.grid(True, which="both")
fig1.tight_layout()  # otherwise the right y-label is slightly clipped

freqs, psd = signal.welch(L/1000,1/dt,'hanning', N/8, N/16)
line=ax3.loglog(freqs,psd,color=color)
ax3.set_ylabel("power.(km$^2$yr$^{-1}$)")
ax3.set_xlabel("frequency.(yr$^{-1}$)")
ax3.grid(True, which="both")
fig3.tight_layout()

for i in range(len(τ)):
    H,L,M,Q,Qg=simulation.Oerlemans_simulation(H0,L0,m0,b0,Ω,smb,N,τ[i])
    plot.figure02(fig1,ax1,ax2,H,L/1000,N,dt)
    plot.figure01(fig2,ax,M,N,dt,label=r'slope')
    plot.figure03(fig3,ax3,L/1000,N,dt)



fig1.savefig('figures/simulationSMB01.png')
fig2.savefig('figures/simulationSMB02.png')
fig3.savefig('figures/simulationSMB03.png')