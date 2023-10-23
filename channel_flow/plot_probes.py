import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
    'font.size': 35
})

def filename(variable,point,simname):
    return variable+point+"-"+simname+".pdf"

A = np.loadtxt(sys.argv[1],unpack=True)
simname = sys.argv[2]
point = sys.argv[3]


t = A[0]
p = A[1]
u = A[2]
v = A[3]
w = A[4]

fig1,ax1 = plt.subplots(figsize=(20,20))
ax1.plot(t,p,'-dk',linewidth=2.0,markersize=5)
ax1.grid()
ax1.set_xlabel("$t$ [s]",fontsize=40)
ax1.set_ylabel("$p$",fontsize=40)
plt.savefig(filename("p/",point,simname),bbox_inches='tight')

fig2,ax2 = plt.subplots(figsize=(20,20))
ax2.plot(t,u,'-dk',linewidth=2.0,markersize=5)
ax2.grid()
ax2.set_xlabel("$t$ [s]",fontsize=40)
ax2.set_ylabel("$u_x$",fontsize=40)
plt.savefig(filename("ux/",point,simname),bbox_inches='tight')

fig3,ax3 = plt.subplots(figsize=(20,20))
ax3.plot(t,v,'-dk',linewidth=2.0,markersize=5)
ax3.grid()
ax3.set_xlabel("$t$ [s]",fontsize=40)
ax3.set_ylabel("$u_y$",fontsize=40)
plt.savefig(filename("uy/",point,simname),bbox_inches='tight')

fig4,ax4 = plt.subplots(figsize=(20,20))
ax4.plot(t,w,'-dk',linewidth=2.0,markersize=5)
ax4.grid()
ax4.set_xlabel("$t$ [s]",fontsize=40)
ax4.set_ylabel("$u_z$",fontsize=40)
plt.savefig(filename("uz/",point,simname),bbox_inches='tight')
