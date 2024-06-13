import numpy as np
import math
from calc_pack import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


k = data.k
a = math.ceil((360-50)/k)

phi = np.array([0]*a, dtype="float64")
p = np.array([0]*a, dtype="float64")
V = np.array([0]*a, dtype="float64")
T = np.array([0]*a, dtype="float64")
dL = np.array([0]*a, dtype="float64")
dT = np.array([0]*a, dtype="float64")
dQw = np.array([0]*a, dtype="float64")
dQx = np.array([0]*a, dtype="float64")
cv = np.array([0]*a, dtype="float64")
phi1 = np.array([0]*a, dtype="float64")
p1 = np.array([0]*a, dtype="float64")
V1 = np.array([0]*a, dtype="float64")
T1 = np.array([0]*a, dtype="float64")
dL1 = np.array([0]*a, dtype="float64")
dT1 = np.array([0]*a, dtype="float64")
dQw1 = np.array([0]*a, dtype="float64")
cv1 = np.array([0]*a, dtype="float64")

phi[0] = data.phi0
p[0] = data.p1
T[0] = data.T1
V[0] = calc_Volume(phi[0])
R = data.R
m = calc_m(V[0])
G_1 = calc_G_1(m)
dt = data.deltat

L = 0
shk = np.arange(start=0, stop=a, step=k)
for i in range(0, a-1, k):
    dH = 0
    dU = 0
    cv[i] = calc_cv(T[i], phi[i], m)
    dL[i] = calc_dL(p[i], phi[i])
    dQw[i] = calc_dQw(phi[i], T[i], V[i])
    dQx[i] = calc_dQx(phi[i])

    dT[i] = (dQw[i]+dQx[i]+dH-dU-dL[i])/m/cv[i]

    phi[i+1] = phi[i]+k
    V[i+1] = calc_Volume(phi[i+1])
    T[i+1] = T[i]+dT[i]*dt
    p[i+1] = m*R*T[i+1]/V[i+1]
    dL[i] = (V[i+1]-V[i])*(p[i]+p[i+1])/2
    L = L+dL[i]

# Графики

fig = plt.figure(figsize=(20, 20))
gs = GridSpec(ncols=3, nrows=3, figure=fig)
phi_graph = plt.subplot(gs[0, 0])
phi_graph.plot(shk, phi)

phi_graph.set_ylabel("phi")
phi_graph.set_title('Phi')
plt.grid()
dL_graph = plt.subplot(gs[0, 1])
dL_graph.plot(shk, dL)

dL_graph.set_ylabel("dL")
dL_graph.set_title("dL")
plt.grid()
dV_graph = plt.subplot(gs[0, 2])
dV_graph.plot(shk, V)

dV_graph.set_ylabel("dV")
dV_graph.set_title("dV")
plt.grid()

dQx_graph = plt.subplot(gs[1, 1])
dQx_graph.plot(shk, dQx)

dQx_graph.set_ylabel("dQx")
dQx_graph.set_title("dQx")
plt.grid()
cv_graph = plt.subplot(gs[1, 2])
cv_graph.plot(shk[:-1], cv[:-1])

cv_graph.set_ylabel("cv")
cv_graph.set_title("cv")
plt.grid()
T_graph = plt.subplot(gs[2, 0])
T_graph.plot(shk[:-1], T[:-1])

T_graph.set_ylabel("T")
T_graph.set_title("T")
plt.grid()
dQw_graph = plt.subplot(gs[2, 1])
dQw_graph.plot(shk[:-1], dQw[:-1])

dQw_graph.set_ylabel("dQw")
dQw_graph.set_title("dQw")
plt.grid()
pV_graph = plt.subplot(gs[2, 2])
pV_graph.plot(V, p)
pV_graph.set_xlabel("V")
pV_graph.set_ylabel("p")
pV_graph.set_title("pV")
plt.grid()

Vh = data.Vh
Hu = data.Hu
i = data.i
n = data.n
omega = data.omega

tau = 4                                 # Четыре такта     
p_ind = L/Vh                            # Индикаторное давление
N_ind = p_ind * i * Vh * n/60 * 2/tau   # Индикаторная мощность    
f = 3600/(720*np.pi/180/omega)
G = (G_1)*f                             # Подача топлива
g_ind = G*i/N_ind                       # Индикаторная подача топлива    
ef_ind = 3600/Hu/g_ind                  # Индикаторный КПД    

# Проведём расчёт параметров без сгорания топлива 

phi1[0] = data.phi0
p1[0] = data.p1
T1[0] = data.T1
V1[0] = calc_Volume(phi1[0])

R = data.R
m = calc_m(V1[0])
dt = data.deltat
print(f"L = {L}", end="\n")
L = 0
for i in range(0, a-1):
    dH = 0
    dU = 0
    cv1[i] = calc_cv(T1[i], phi1[i], m)
    dL1[i] = calc_dL(p1[i], phi1[i])
    dQw1[i] = calc_dQw(phi1[i], T1[i], V1[i])
    dQx1 = 0
    dT1[i] = (dQw1[i]+dQx1+dH-dU-dL1[i])/m/cv1[i]
    phi1[i+1] = phi1[i]+k
    V1[i+1] = calc_Volume(phi1[i+1])

    T1[i+1] = T1[i]+dT1[i]*dt
    p1[i+1] = m*R*T1[i+1]/V1[i+1]

# Графики

p_graph = plt.subplot(gs[1, 0])
p_graph.plot(shk, p)
p_graph.plot(shk, p1)
p_graph.set_xlabel("k")
p_graph.set_ylabel("p")
p_graph.set_title("p")
p_graph.legend(['С учетом сгорания', 'Без учета сгорания'])
plt.grid()
fig2 = plt.figure(figsize=(20, 10))
gs2 = GridSpec(ncols=3, nrows=1, figure=fig2)
P_graph = plt.subplot(gs[0, 0])
P_graph.plot(phi, p)
P_graph.set_xlabel("phi")
P_graph.set_ylabel("P")
P_graph.set_title('Давление от угла поворота коленчатого вала')
plt.grid()
T_graph = plt.subplot(gs[0, 1])
T_graph.plot(phi, T)
T_graph.set_xlabel("phi")
T_graph.set_ylabel("T(phi)")
T_graph.set_title('Температура от угла поворота коленчатого вала')
plt.grid()
PP_graph = plt.subplot(gs[0, 2])
PP_graph.plot(V, p)
PP_graph.set_xlabel("V")
PP_graph.set_ylabel("P(V)")
PP_graph.set_title('Индикаторная диаграмма от угла поворота коленчатого вала')
plt.grid()

# Вывод параметров

print(f"p_ind = {p_ind}",  f"N_ind = {N_ind}",
      f"g_ind = {g_ind}", f"ef_ind = {ef_ind}", sep="\n")
plt.show()
