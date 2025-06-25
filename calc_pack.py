import numpy as np
import math
from dataclasses import dataclass
from typing import Dict, Tuple

@dataclass
class Data:
    """Класс параметров двигателя"""
    D: float = 123e-3
    S: float = 145e-3
    E: float = 20
    lamda: float = 0.25
    i: int = 4
    xz: float = 0.98
    phiz: float = 70
    mz: float = 3.2
    n: float = 1500
    alphasg: float = 550
    alphasj: float = 1700
    alpha: float = 1.73
    phi0: float = -120
    phisg: float = -50
    T1: float = 355
    p1: float = 1e5
    Hu: float = 44.0e6
    gC: float = 0.855
    gH: float = 0.145
    gO: float = 0.0
    R: float = 287.1
    k: int = 1

    def __post_init__(self):
        self.A = np.pi*self.D**2/4
        self.Vh = self.A*self.S
        self.Vc = self.Vh/(self.E-1)
        self.Va = self.Vc + self.Vh
        self.omega = np.pi*self.n/30
        self.deltat = self.k*np.pi/(180*self.omega)

@dataclass
class SimulationResults:
    """Результаты моделирования"""
    # Основные массивы
    phi: np.ndarray
    p: np.ndarray
    V: np.ndarray
    T: np.ndarray
    dL: np.ndarray
    dQw: np.ndarray
    dQx: np.ndarray
    cv: np.ndarray
    
    # Для сравнения (без сгорания)
    phi1: np.ndarray
    p1: np.ndarray
    V1: np.ndarray
    T1: np.ndarray
    
    # Индикаторные параметры
    p_ind: float
    N_ind: float
    g_ind: float
    ef_ind: float
    L: float
    
    # Дополнительные параметры
    shk: np.ndarray
    m: float
    G_1: float

def calc_Volume(phi: float, data: Data) -> float:
    R = data.S/2
    L = data.lamda*R
    phi_rad = np.deg2rad(phi)
    beta = np.arcsin(data.lamda*np.sin(phi_rad))
    return data.Vc + data.A*(R*(1-np.cos(phi_rad)) + L*(1-np.cos(beta)))

def calc_m(V: float, data: Data) -> float:
    return data.p1*V/(data.R*data.T1)

def calc_G_1(m: float, data: Data) -> float:
    return m/(14.5*data.alpha)

def calc_dL(p: float, phi: float, data: Data) -> float:
    R = data.S/2
    phi_rad = np.deg2rad(phi)
    return p*data.A*data.omega*R*(np.sin(phi_rad)+data.lamda/2*np.sin(2*phi_rad))

def calc_dQw(phi: float, T: float, V: float, data: Data) -> float:
    Twx = 135+273.15
    Twp = 325+273.15
    Twg = 325+273.15
    Awp = data.A
    Awg = data.A + data.Vc/data.A*2*np.pi*data.D/2
    Awx = 2*np.pi*data.D*(V-data.Vc)/data.A
    alpha = data.alphasg if phi < 0 else data.alphasj
    return alpha*(Awp*(Twp-T) + Awg*(Twg-T) + Awx*(Twx-T))

def calc_dQx(phi: float, data: Data) -> float:
    if phi < data.phisg:
        return 0.0
    betta = phi - data.phisg
    c = np.log(1-data.xz)
    dxdt = (-np.exp(c*betta**(data.mz+1)/data.phiz**(data.mz+1)) * 
            c*(data.mz+1)*betta**data.mz/data.phiz**(data.mz+1) * 
            180/np.pi/data.omega)
    return data.Hu*data.G_1*dxdt

def calc_cv(T: float, phi: float, m: float, data: Data) -> float:
    M_AIR = 28.97e-3
    if phi <= data.phisg:
        return (20.6 + 0.002638*(T-273.15))/M_AIR
    
    ksi = phi - data.phisg
    c = np.log(1-data.xz)
    x = 1-np.exp(c*(ksi/(data.phiz-26))**(data.mz+1))
    mt = data.G_1*x
    
    m_CO2 = data.gC*mt*3.67
    m_H2O = data.gH*mt*9
    m_O2 = data.gO*mt*2
    mall = mt + m
    
    g_CO2 = m_CO2/mall
    g_H2O = m_H2O/mall
    g_VOZD = 1 - g_CO2 - g_H2O
    
    cv_air = (20.6 + 0.002638*(T-273.15))/M_AIR
    cv_CO2 = (27.941 + 0.019*(T-273.15) - 5.487e-6*(T-273.15)**2)/44.01e-3
    cv_H2O = (24.953 + 0.05359*(T-273.15))/18.01e-3
    
    return g_VOZD*cv_air + g_CO2*cv_CO2 + g_H2O*cv_H2O

def run_simulation(data: Data) -> SimulationResults:
    a = math.ceil((360-50)/data.k)
    shk = np.arange(start=0, stop=a, step=data.k)
    
    # Инициализация массивов
    arrays = {
        'phi': np.zeros(a), 'p': np.zeros(a), 'V': np.zeros(a),
        'T': np.zeros(a), 'dL': np.zeros(a), 'dQw': np.zeros(a),
        'dQx': np.zeros(a), 'cv': np.zeros(a),
        'phi1': np.zeros(a), 'p1': np.zeros(a), 'V1': np.zeros(a),
        'T1': np.zeros(a), 'dL1': np.zeros(a), 'dQw1': np.zeros(a),
        'cv1': np.zeros(a)
    }
    
    # Начальные условия
    arrays['phi'][0] = data.phi0
    arrays['p'][0] = data.p1
    arrays['T'][0] = data.T1
    arrays['V'][0] = calc_Volume(arrays['phi'][0], data)
    m = calc_m(arrays['V'][0], data)
    G_1 = calc_G_1(m, data)
    data.G_1 = G_1
    
    # Основной расчёт
    L = 0
    for i in range(a-1):
        arrays['cv'][i] = calc_cv(arrays['T'][i], arrays['phi'][i], m, data)
        arrays['dL'][i] = calc_dL(arrays['p'][i], arrays['phi'][i], data)
        arrays['dQw'][i] = calc_dQw(arrays['phi'][i], arrays['T'][i], arrays['V'][i], data)
        arrays['dQx'][i] = calc_dQx(arrays['phi'][i], data)
        
        dT = (arrays['dQw'][i] + arrays['dQx'][i] - arrays['dL'][i])/(m*arrays['cv'][i])
        
        arrays['phi'][i+1] = arrays['phi'][i] + data.k
        arrays['V'][i+1] = calc_Volume(arrays['phi'][i+1], data)
        arrays['T'][i+1] = arrays['T'][i] + dT*data.deltat
        arrays['p'][i+1] = m*data.R*arrays['T'][i+1]/arrays['V'][i+1]
        arrays['dL'][i] = (arrays['V'][i+1]-arrays['V'][i])*(arrays['p'][i]+arrays['p'][i+1])/2
        L += arrays['dL'][i]
    
    # Расчёт без сгорания
    arrays['phi1'][0] = data.phi0
    arrays['p1'][0] = data.p1
    arrays['T1'][0] = data.T1
    arrays['V1'][0] = calc_Volume(arrays['phi1'][0], data)
    
    for i in range(a-1):
        arrays['cv1'][i] = calc_cv(arrays['T1'][i], arrays['phi1'][i], m, data)
        arrays['dL1'][i] = calc_dL(arrays['p1'][i], arrays['phi1'][i], data)
        arrays['dQw1'][i] = calc_dQw(arrays['phi1'][i], arrays['T1'][i], arrays['V1'][i], data)
        
        dT = (arrays['dQw1'][i] - arrays['dL1'][i])/(m*arrays['cv1'][i])
        
        arrays['phi1'][i+1] = arrays['phi1'][i] + data.k
        arrays['V1'][i+1] = calc_Volume(arrays['phi1'][i+1], data)
        arrays['T1'][i+1] = arrays['T1'][i] + dT*data.deltat
        arrays['p1'][i+1] = m*data.R*arrays['T1'][i+1]/arrays['V1'][i+1]
    
    # Индикаторные параметры
    p_ind = L/data.Vh
    N_ind = p_ind * data.i * data.Vh * data.n/60 * 0.5
    f = 3600/(720*np.pi/180/data.omega)
    G = G_1*f
    g_ind = G*data.i/N_ind
    ef_ind = 3600/data.Hu/g_ind
    
    return SimulationResults(
        phi=arrays['phi'], p=arrays['p'], V=arrays['V'], T=arrays['T'],
        dL=arrays['dL'], dQw=arrays['dQw'], dQx=arrays['dQx'], cv=arrays['cv'],
        phi1=arrays['phi1'], p1=arrays['p1'], V1=arrays['V1'], T1=arrays['T1'],
        p_ind=p_ind, N_ind=N_ind, g_ind=g_ind, ef_ind=ef_ind, L=L,
        shk=shk, m=m, G_1=G_1
    )