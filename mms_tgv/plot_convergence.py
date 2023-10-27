import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
    'font.size': 35
})

EULER = np.loadtxt("euler.dat",unpack=True)
RK2 = np.loadtxt("rk2.dat",unpack=True)
RK3 = np.loadtxt("rk3.dat",unpack=True)
RK4 = np.loadtxt("rk4.dat",unpack=True)

fig, ax = plt.subplots(figsize=(20,20))
ax.loglog(EULER[0],EULER[1],'-ok',label='Euler',linewidth=2.0,markersize=10)
ax.loglog(RK2[0],RK2[1],'-dk',label='Heun RK2',linewidth=2.0,markersize=10)
ax.loglog(RK3[0],RK3[1],'-<k',label='Heun RK3',linewidth=2.0,markersize=10)
ax.loglog(RK4[0],RK4[1],'-sk',label='Standard RK4',linewidth=2.0,markersize=10)
ax.loglog([1e-1,1e-2],[5e-4,5e-5],'--k')
ax.loglog([1e-1,1e-2],[1e-7,1e-9],'--k')
ax.loglog([1,1e-1],[2e-8,2e-11],'--k')
ax.loglog([1,1e-1],[4e-11,4e-15],'--k')
ax.legend(fontsize=35)
ax.set_xlabel('$\Delta t$',fontsize=40)
ax.set_ylabel('$L_\infty(|u-u_a|)$',fontsize=40)
ax.text(3e-2,3e-4,"$\Delta t$",fontsize=35)
ax.text(3e-2,3e-8,"$\Delta t^2$",fontsize=35)
ax.text(2e-1,1e-9,"$\Delta t^3$",fontsize=35)
ax.text(3e-1,3e-12,"$\Delta t^4$",fontsize=35)
ax.grid()

plt.savefig("convergence_tf2-rk.pdf",bbox_inches='tight')
