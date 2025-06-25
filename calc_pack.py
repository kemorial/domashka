import numpy as np
import math 


class Data:
    def __init__(self):  
        self.D = 123e-3    # Диаметр цилиндра (м.)
        self.S = 145e-3     # Ход поршня (м.)
        self.E = 20         # Степень сжатия
        self.lamda = 0.25
        self.i = 4        # Количество цилиндров
        self.xz = 0.98     
        self.phiz = 70     # Продолжительность сгорания (град.)
        self.mz = 3.2
        self.n = 1500      # Частота вращения (мин^-1)
        self.alphasg = 550  # Альфа сжатия
        self.alphasj = 1700  # Альфа сгорания
        self.alpha = 1.73    # Коэффициент избытка воздуха
        self.phi0 = -120       # Угол закрытия впускного клапана
        self.phisg = -50     # Угол опережения зажигания
        self.T1 = 355      # Начальная температура (К)
        self.p1 = 1e5       # Начальное давление (Па)
        self.Hu = 44.0e6     # Низшая теплотворность (Дж/кг)
        self.gC = 0.855         # % Углерода в ДТ
        self.gH = 0.145      # % Водорода в ДТ
        self.gO = 0.0        # % Кислорода в ДТ
        self.R = 287.1        
        A = np.pi*self.D*self.D/4      # Площадь цилиндра (м^2)
        Vh = A*self.S                  # Рабочий объём (м^3)
        Vc = Vh/(self.E-1)             # Остаточный объём
        Va = Vc+Vh                     # Полный объём
        k = 1                           # Шаг по углу (град)
        omega = np.pi*self.n / 30       # Угловая скорость (рад/с)
        dt = k*np.pi/(180*omega)        # Шаг по времени (с)
        self.k = k
        self.Vh = Vh
        self.A = A
        self.Vc = Vc
        self.Va = Va
        self.omega = omega
        self.deltat = dt


data = Data()

# Функция расчёта объёма
def calc_Volume(phi_deg):

    R = data.S/2
    L = data.lamda*R
    phi_deg = np.deg2rad(phi_deg)
    beta = np.arcsin(data.lamda*np.sin(phi_deg))
    Sh = R*(1-np.cos(phi_deg))+L*(1-np.cos(beta))
    V = data.Vc + np.pi*(data.D**2)/4*Sh
    return V

# Функция расчёта массы
def calc_m(V):
    p = data.p1
    T = data.T1
    R = data.R
    m = p*V/(R*T)
    return m

# Функция расчёта цикловой подачи топлива
def calc_G_1(m):
    L0 = 14.5
    alpha = data.alpha
    G_1 = m/(L0*alpha)
    data.G_1 = G_1
    return G_1


def calc_dV(phi):
    A = data.A
    omega = data.omega
    lamda = data.lamda
    R = data.S / 2
    L = R/lamda
    k = data.k
    alpha = phi + k
    alpha = np.deg2rad(alpha)
    phi = np.deg2rad(phi)
    dV = A*omega*R*(np.sin(phi)+lamda/2*np.sin(2*phi))
    return dV


def calc_dL(p, phi):
    dV = calc_dV(phi)
    dL = p*dV
    return dL


def calc_alpha(phi):
    alphasg = data.alphasg
    alphasj = data.alphasj
    if (phi < 0):
        alpha = alphasg
    else:
        alpha = alphasj
    return alpha


def calc_Awx(V):
    A = data.A
    Vc = data.Vc
    D = data.D

    Vx = V-Vc
    dx = Vx/A
    Awx = 2*np.pi*D*dx
    return Awx


def calc_dQw(phi, T, V):
    Twx = 135+273.15  # Температура цилиндра
    Twp = 325+273.15  # Температура поршня
    Twg = 325+273.15  # Температура головки

    Awp = data.A
    Awg = data.A+data.Vc/data.A*2*np.pi*data.D/2

    Awx = calc_Awx(V)
    alpha = calc_alpha(phi)
    dQw = alpha*(Awp*(Twp-T)+Awg*(Twg-T)+Awx*(Twx-T))
    return dQw

# Формула Вибе
def calc_Wibe(phi):
    phiz = data.phiz
    phisg = data.phisg
    xz = data.xz
    m = data.mz
    omega = data.omega

    betta = (phi-phisg)

    c = np.log(1-xz)
    x = 1-np.exp(c*((betta/phiz)**(m+1)))

    dx = ((-((np.exp((c*betta**(m+1))/(phiz**(m+1)))*c*m*betta**m+np.exp((c *
          betta**(m+1))/(phiz**(m+1)))*c*betta**m))/(phiz**(m+1)))/(1/omega))*180/np.pi

    return dx


def calc_dQx(phi):
    Hu = data.Hu
    phisg = data.phisg
    phiz = data.phiz
    mdt = data.G_1

    if (phi >= phisg):
        return Hu*mdt*calc_Wibe(phi)
    else:
        return 0


def calc_cv(T: float, phi: float, m: float, data:Data) -> float:
    """
    Расчет удельной теплоемкости (cv) рабочего тела в цилиндре ДВС.
    
    Параметры:
        T - температура [K]
        phi - текущий угол поворота коленвала [град]
        m - масса воздуха в цилиндре [кг]
        data - параметры двигателя
        
    Возвращает:
        cv - удельная теплоемкость [Дж/(кг·K)]
    """
    # Константы
    M_AIR = 28.97e-3    # Молярная масса воздуха [кг/моль]
    M_CO2 = 44.01e-3    # Молярная масса CO2 [кг/моль]
    M_H2O = 18.01e-3    # Молярная масса H2O [кг/моль]
    
    if phi <= data.phisg:
        # Только воздух до начала сгорания
        cv_air = (20.6 + 0.002638*(T-273.15)) / M_AIR
        return cv_air
    
    # Расчет доли сгоревшего топлива (закон Вибе)
    ksi = phi - data.phisg
    c = np.log(1 - data.xz)
    x = 1 - np.exp(c * (ksi/(data.phiz-26))**(data.mz+1))
    mt = data.G_1 * x
    
    # Массы компонентов продуктов сгорания
    m_CO2 = data.gC * mt * 3.67
    m_H2O = data.gH * mt * 9
    m_O2 = data.gO * mt * 2
    mall = mt + m
    
    # Массовые доли компонентов
    g_CO2 = m_CO2 / mall
    g_H2O = m_H2O / mall
    g_O2 = m_O2 / mall
    g_air = 1 - g_CO2 - g_H2O - g_O2
    
    # Расчет теплоемкостей компонентов
    def cv_poly(T, coeffs):
        """Полиномиальная аппроксимация удельной теплоемкости"""
        return sum(c*(T-273.15)**i for i, c in enumerate(coeffs))
    
    # Коэффициенты для полиномов cv(T) [Дж/(кг·K)]
    # Примерные значения - следует уточнить для конкретного случая
    cv_coeffs = {
        'air': [717, 0.07, -1e-5],
        'CO2': [652, 0.85, -2.5e-4],
        'H2O': [1437, 0.32, -5e-5]
    }
    
    cv_air = cv_poly(T, cv_coeffs['air'])
    cv_CO2 = cv_poly(T, cv_coeffs['CO2'])
    cv_H2O = cv_poly(T, cv_coeffs['H2O'])
    
    # Усредненная теплоемкость смеси
    cv_mix = (g_air * cv_air + 
              g_CO2 * cv_CO2 + 
              g_H2O * cv_H2O)
    
    return cv_mix
