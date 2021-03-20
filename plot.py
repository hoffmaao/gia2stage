import matplotlib.pyplot as plt

w = 5
h = 3
d = 270


def figure01(n,dt,h,L,w=5,h=3,d=270,filename='figures/fig01.png'):
	t=np.linspace(0,nts*dt,nts)
	fig,ax1= plt.subplots( figsize=(w,h), dpi=d)
	plt.xlabel('iterate')
	x=np.linspace(0,n*dt,n)
	color = 'tab:red'
	ax1.set_xlabel('iterate')
	ax1.set_ylabel('thickness.(m)', color=color)
	ax1.plot(x,h, color=color)
	ax1.tick_params(axis='y', labelcolor=color)

	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
	color = 'tab:blue'
	ax2.set_ylabel('length.(km)', color=color)  # we already handled the x-label with ax
	ax2.plot(x, L/1000, color=color)
	ax2.tick_params(axis='y', labelcolor=color)
	fig.tight_layout()  # otherwise the right y-label is slightly clipped
	plt.savefig(fielname)