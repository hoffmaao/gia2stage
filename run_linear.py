import numpy as np
import forcing
import simulation
import plot
import model

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
from scipy import signal


# Simulation parameters
N = 500000
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

fig1,ax1,ax2,plots=plot.makefig02(1)
fig2,ax,plots=plot.makefig01(1)
fig3,ax3,plots=plot.makefig01(1)

b0 = -100
mbar = -2e-3
Lbar = 182000
Hbar = 1400
κ = -2*mbar/(Lbar*(-2*b0*ρ_W-2*ρ_I*Hbar+ρ_W*Lbar*mbar))
τ = 3000*year

H,L,Q,Qg=simulation.simulation(Hbar,Lbar,mbar,b0,Ω,smb,N)

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

Hgbar = model.grounding_thickness(model.bed(b0,mbar,Lbar))
Qgbar = model.grounding_flux(1.0,b0,mbar,Lbar,(1/n+n+3)/(1/n+1))
H,L,M=simulation.linear_simulation(Hbar, Lbar, mbar, Hgbar, Qgbar, b0, Ω, smb-Sbar, τ, κ, N)


plot.figure02(fig1,ax1,ax2,Hbar+H, (Lbar+L)/1000,N,dt)
plot.figure01(fig2,ax,M,N,dt,label=r'slope')
plot.figure03(fig3,ax3,L/1000,N,dt)


fig1.savefig('figures/linearSMB01.png')
fig2.savefig('figures/linearSMB02.png')
fig3.savefig('figures/linearSMB03.png')

