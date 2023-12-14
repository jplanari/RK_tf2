import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
    'font.size': 50
})

dts = [0.8,1.2,1.8,2.3,2.6]
linestyle = ['->k','-^k','-sk','-dk','--k']
tmax = float(sys.argv[3])

ms=20
lw=4

fig,ax = plt.subplots(figsize=(20,20))
for i in range(len(dts)):
    print(tmax/dts[i]-0.009)
    n = int(np.round(tmax/dts[i]-0.0999999999))
    A = np.loadtxt(sys.argv[1]+"/results/"+sys.argv[2]+"/testPositivity_"+sys.argv[2]+"_h_"+str(dts[i])+"_N_"+str(n)+".dat",unpack=True)
    ax.plot(A[0],A[1],linestyle[i],linewidth=lw,markersize=ms,label='$h=$'+str(dts[i]))

ax.set_xlabel("$t[s]$",fontsize=55)
ax.set_ylabel("$y(t)$",fontsize=55)
ax.legend(fontsize=60)
ax.set_xlim(0,15)
ax.set_ylim(-1,1)
ax.grid()
plt.savefig(sys.argv[1]+"_"+sys.argv[2]+".pdf",bbox_inches='tight')

