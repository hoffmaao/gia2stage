import matplotlib.pyplot as plt



def makefig01(plots,w=5,h=3,d=270):
	fig,ax= plt.subplots(figsize=(w,h), dpi=d)
	return fig,ax,plots

def makefig02(plots,w=5,h=3,d=270):
	fig,ax1,ax2= plt.subplots(2,1,figsize=(w,h), dpi=d)
	return fig,ax1,ax2,plots


def figure01(fig,ax,N,dt,data,label='thickness.(m)',plot=1):
	t=np.linspace(0,N*dt,n)
	plt.xlabel('iterate')
	x=np.linspace(0,n*dt,N)
	color = next(ax._get_lines.prop_cycler)['color']
	ax.set_xlabel('iterate')
	ax.set_ylabel(label, color=color)
	ax.plot(x,data, color=color)
	fig.tight_layout()  # otherwise the right y-label is slightly clipped


def figure02(fig,ax1,ax2,N,dt,data1,data2,label1='thickness.(m)',label2='length.(km)'):
	t=np.linspace(0,N*dt,N)
	plt.xlabel('iterate')
	x=np.linspace(0,N*dt,N)
	color = next(ax._get_lines.prop_cycler)['color']

	ax1.set_xlabel('iterate')
	ax1.set_ylabel('thickness.(m)', color=color)
	ax1.plot(x,data1, color=color)
	ax1.tick_params(axis='y', labelcolor=color)

	ax2.set_ylabel('length.(km)', color=color)  # we already handled the x-label with ax
	ax2.plot(x, data2, color=color)
	ax2.tick_params(axis='y', labelcolor=color)
	fig.tight_layout()  # otherwise the right y-label is slightly clipped
