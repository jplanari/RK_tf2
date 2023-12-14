import numpy as np
import matplotlib.pyplot as plt
import sys

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Computer Modern Roman",
    'font.size': 50
})


def z(ev,phi,dt,s):
    real = 0.0
    imag = 0.0

    for k in range(s+1):
        real = real + (-1.0)**k*dt**k*np.cos(k*phi)/np.math.factorial(k)
        imag = imag + (-1.0)**(k+1)*dt**k*np.sin(k*phi)/np.math.factorial(k)

    zs = real**2+imag**2
    return [real,imag,zs]

phi = np.linspace(0,np.pi/2,1000)
dt = np.linspace(0.001,10,10000)

real_stability = np.zeros_like(phi,dtype=float)
imag_stability = np.zeros_like(phi,dtype=float)

real_pos = np.zeros_like(phi,dtype=float)
imag_pos = np.zeros_like(phi,dtype=float)

real_phase = np.zeros_like(phi,dtype=float)
imag_phase = np.zeros_like(phi,dtype=float)



s=int(sys.argv[1])

for i in range(len(phi)):
    vals = z(-1.0,phi[i],dt,s)
    dts = dt[np.where(vals[2] >= 1.0)[0][0]]
    
    real_stability[i] = -dts*np.cos(phi[i])
    imag_stability[i] = dts*np.sin(phi[i])

    try:
        dtr = dt[np.where(vals[0] <= 0.0)[0][0]]
        if dtr == 0.001:
            dtr = 10
    except IndexError:
        dtr = 10
    
    real_pos[i] = -dtr*np.cos(phi[i])
    imag_pos[i] = dtr*np.sin(phi[i])

    try:
        dti = dt[np.where(vals[1] <= 0.0)[0][0]]
        if dti == 0.001:
            dti = 10
    except IndexError:
        dti = 10

    real_phase[i] = -dti*np.cos(phi[i])
    imag_phase[i] = dti*np.sin(phi[i])


r_stab = real_stability**2+imag_stability**2
r_pos = real_pos**2+imag_pos**2
r_phase = real_phase**2+imag_phase**2


fig, ax = plt.subplots(figsize=(20,20))
ax.plot(real_stability,imag_stability,'-k',linewidth=5.0,label='Stability region')
ax.plot(real_pos,imag_pos,'--k',linewidth=5.0,label='Positivity region')
ax.plot(real_phase,imag_phase,'-.k',linewidth=5.0,label='Phase region')
ax.set_xlim(-3.5,0)
ax.set_ylim(0,3.5)
ax.set_xlabel("Re",fontsize=55)
ax.set_ylabel("Im",fontsize=55)
ax.grid()

aval=0.8
col='skyblue'

C1 = r_stab>r_phase
C2 = r_phase>r_pos
C3 = r_stab>r_pos

condition5=real_pos<-1.0

conditions_pos = (~C1&C3)|(C1&C2)
conditions_phase = (C3&~C2)|(~C3&C1)
conditions_stab = (C2&~C3)|(~C3&~C1)

ax.fill_between(real_pos, 0, imag_pos,
        where=conditions_pos, color=col, alpha=aval, label='Go-to zone')
ax.fill_between(real_phase, 0, imag_phase, where=conditions_phase,
        color=col, alpha=aval)
ax.fill_between(real_stability, 0, imag_stability,
        where=conditions_stab, color=col, alpha=aval)

ax.fill_between(real_pos, 0, imag_pos, where=(C3&condition5), color='white')
# Blank the region between Stability region and Phase region
#condition_blank = (imag_stability < imag_phase)
#ax.fill_between(real_stability, 0, imag_phase, where=condition_blank, color='white')

ax.legend(fontsize=60)
plt.savefig("plots/regions_"+str(s)+".pdf",bbox_inches='tight')



