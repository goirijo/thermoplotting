from pylab import *

rc('axes', linewidth=5)
rc('axes', grid=True)
rc('font', weight='heavy')
rc('font', size=16)
rc('xtick.major', size=6) 
rc('xtick.major', width=3) 
rc('ytick.major', size=6) 
rc('ytick.major', width=3) 

plt.xlim([0.5,1.0])
plt.ylim([400,1400])

fontsize=16
ax=gca()

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('heavy')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('heavy')

plt.show()
