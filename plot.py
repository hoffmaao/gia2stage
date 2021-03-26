import matplotlib.pyplot as plt
import numpy as np


def makefig01(plots,w=5,h=3,d=270):
	fig,ax= plt.subplots(figsize=(w,h), dpi=d)
	return fig,ax,plots

def makefig02(plots,w=5,h=3,d=270):
	fig,(ax1,ax2)= plt.subplots(2,1,figsize=(w,h), dpi=d)
	return fig,ax1,ax2,plots


def figure01(fig,ax,data,N,dt,label='thickness.(m)',plot=1):
	x=np.linspace(0,N*dt,N)
	color = next(ax._get_lines.prop_cycler)['color']
	ax.set_xlabel('year.(a)')
	ax.set_ylabel(label, color=color)
	ax.plot(x,data, color=color)
	fig.tight_layout()  # otherwise the right y-label is slightly clipped


def figure02(fig,ax1,ax2,data1,data2,N,dt,label1='thickness.(m)',label2='length.(km)'):

	x=np.linspace(0,N*dt,N)
	color = next(ax1._get_lines.prop_cycler)['color']

	ax2.set_xlabel('year.(a)')
	ax1.set_ylabel('thickness.(m)')
	ax1.plot(x,data1, color=color)

	ax2.set_ylabel('length.(km)')  # we already handled the x-label with ax
	ax2.plot(x, data2, color=color)
	fig.tight_layout()  # otherwise the right y-label is slightly clipped
