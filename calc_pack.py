import numpy as np
import math 


class Data:
    def __init__(self):  
        self.D = 168e-3    # Диаметр цилиндра (м.)
        self.S = 180e-3     # Ход поршня (м.)
        self.E = 23         # Степень сжатия
        self.lamda = 0.31
        self.i = 12        # Количество цилиндров
        self.xz = 0.99     
        self.phiz = 70     # Продолжительность сгорания (град.)
        self.mz = 3.5
        self.n = 1700      # Частота вращения (мин^-1)
        self.alphasg = 800  # Альфа сжатия
        self.alphasj = 1800  # Альфа сгорания
        self.alpha = 2.15    # Коэффициент избытка воздуха
        self.phi0 = -120       # Угол закрытия впускного клапана
        self.phisg = -60     # Угол опережения зажигания
        self.T1 = 355      # Начальная температура (К)
        self.p1 = 1e5       # Начальное давление (Па)
        self.Hu = 42.7e6     # Низшая теплотворность (Дж/кг)
        self.gC = 0.87         # % Углерода в ДТ
        self.gH = 0.126      # % Водорода в ДТ
        self.gO = 0.4        # % Кислорода в ДТ
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
    l0 = 14.9
    alpha = data.alpha
    G_1 = m/(l0*alpha)
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


def calc_cv(T, phi, m):
    phisg = data.phisg
    phiz = data.phiz
    alpha = data.alpha
    xz = data.xz
    mz = data.mz
    gC = data.gC        # углерода в ДТ
    gH = data.gH       # водорода в ДТ
    gO = data.gO         # кислорода в ДТ
    G_1 = data.G_1

    ksi = (phi-phisg)  # новый угол начинется с 0.
    betta = phisg+phiz  # угол конца сгорания
    M = m/(28.97*1e-3)  # кол-во воздуха в моль

    if (phi <= phisg):

        Mcv = 20.6+0.002638*(T-273.15)
        Cv = Mcv*M
        cv = Cv/m

    else:

        c = np.log(1-xz)
        # закон впрыска топлива позакону Вибе
        x = 1-np.exp(c*((ksi/(phiz-26))**(mz+1)))
        mt = G_1*x
        mall = mt+m  # общая масса горючей смеси, кг
        # масса СО2 равна массе С * 3.67, так как молярная масса в 3,67 раз больше и на одну молекулу углерода приходится 2 молекулы кислорода.
        m_CO2 = gC*mt*3.67
        m_H2O = gH*mt*9
        m_O2 = gO*mt*2
        g_O2 = m_O2/mall
        g_H2O = m_H2O/mall
        g_CO2 = m_CO2/mall
        # допущенние этот воздух обедненый, то есть концентрация азота больше чем 70
        g_VOZD = 1 - g_CO2 - g_H2O - g_O2

        Mcv = 20.6+0.002638*(T-273.15)
        M = m/(28.97*1e-3)  # кол-во воздуха в моль
        Cv = Mcv*M
        cv_VOZD = Cv/m

        # для СО2
        Mcv_CO2 = 27.941+0.019*(T-273.15)-0.000005487 * \
            (T-273.15)*(T-273.15)
        M_CO2 = m_CO2/(44.01*1e-3)  # кол-во CO2 в моль
        Cv_CO2 = Mcv_CO2*M_CO2
        cv_CO2 = Cv_CO2/m_CO2

        # для H2O
        Mcv_H2O = 24.953+0.05359*(T-273.15)
        M_H2O = m_H2O/(18*1e-3)  # кол-во CO2 в моль
        Cv_H2O = Mcv_H2O*M_H2O
        cv_H2O = Cv_H2O/m_H2O

        # для O2
        Mcv_O2 = 20.930+0.00464*(T-273.15)-0.00000084 * \
            (T-273.15)*(T-273.15)
        M_O2 = m_O2/(32*1e-3)  # кол-во CO2 в моль
        Cv_O2 = Mcv_O2*M_O2
        cv_O2 = Cv_O2/m_O2
        cv = g_CO2*cv_CO2 + g_H2O*cv_H2O + g_O2*cv_O2 + g_VOZD*cv_VOZD
    return cv
