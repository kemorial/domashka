<?php
define('pi', 3.1415926535898);

class Data
{
    public function __construct()
    {
        (float) $this->D = 168e-3;   # Диаметр цилиндра (м.)
        (float) $this->S = 180e-3;   # Ход поршня (м.)
        (float) $this->E = 23;  # Степень сжатия
        (float) $this->lamda = 0.31;
        (float) $this->i = 12; # Количество цилиндров
        (float) $this->xz = 0.99;
        (float) $this->phiz = 70;# Продолжительность сгорания (град.)
        (float) $this->mz = 3.5;
        (float) $this->n = 1700; # Частота вращения (мин^-1)
        (float) $this->alphasg = 800; # Альфа сжатия
        (float) $this->alphasj = 1800; # Альфа $сгорания
        (float) $this->alpha = 2.15;# Коэффициент избытка воздуха
        (float) $this->phi0 = -120; # Угол закрытия впускного клапана$
        (float) $this->phisg = -60;# Угол опережения зажигания
        (float) $this->T1 = 355;# Начальная температура (К)
        (float) $this->p1 = 1e5; # Начальное давление (Па)
        (float) $this->Hu = 42.7e6;  # Низшая теплотворность (Дж/кг)
        (float) $this->gC = 0.87;   # % Углерода в ДТ
        (float) $this->gH = 0.126;  # % Водорода в ДТ
        (float) $this->gO = 0.4;  # % Кислорода в ДТ
        (float) $this->R = 287.1;
        (float) $this->A = pi * 168e-3 * 168e-3 / 4;    # Площадь цилиндра (м^2)
        (float) $this->Vh = pi * 168e-3 * 168e-3 / 4 * 180e-3;  # Рабочий объём (м^3)
        (float) $this->Vc = pi * 168e-3 * 168e-3 / 4 * 180e-3 / (23 - 1);  # Остаточный объём
        (float) $this->Va = pi * 168e-3 * 168e-3 / 4 * 180e-3 / (23 - 1) + pi * 168e-3 * 168e-3 / 4 * 180e-3;# Полный объём
        (float) $this->k = 1; # Шаг по углу (град)
        (float) $this->omega = pi * 1700 / 30;   # Угловая скорость (рад/с)
        (float) $this->deltat = 1 * pi / (180 * pi * 1700 / 30); # Шаг по времени (с)
    }
}
$data = new Data();

# Функция расчёта объёма
function calc_Volume(&$data, $phi_deg)
{

    $R = $data->S / 2;
    $L = $data->lamda * $R;
    $phi_deg = deg2rad($phi_deg);
    $beta = asin($data->lamda * sin($phi_deg));
    $Sh = $R * (1 - cos($phi_deg)) + $L * (1 - cos($beta));
    $V = $data->Vc + pi * ($data->D ** 2) / 4 * $Sh;
    return $V;
}
# Функция расчёта массы
function calc_m(&$data, $V)
{
    ;
    $p = $data->p1;
    $T = $data->T1;
    $R = $data->R;
    $m = $p * $V / ($R * $T);
    return $m;
}
# Функция расчёта цикловой подачи топлива
function calc_G_1(&$data, $m)
{
    $l0 = 14.9;
    $alpha = $data->alpha;
    $G_1 = $m / ($l0 * $alpha);
    $data->G_1 = $G_1;
    return $G_1;
}

function calc_dV(&$data, $phi)
{
    $A = $data->A;
    $omega = $data->omega;
    $lamda = $data->lamda;
    $R = $data->S / 2;
    $L = $R / $lamda;
    $k = $data->k;
    $alpha = $phi + $k;
    $alpha = deg2rad($alpha);
    $phi = deg2rad($phi);
    $dV = $A * $omega * $R * (sin($phi) + $lamda / 2 * sin(2 * $phi));
    return $dV;
}

function calc_dL(&$data, $p, $phi)
{
    $dV = calc_dV($data, $phi);
    $dL = $p * $dV;
    return $dL;
}

function calc_alpha(&$data, $phi)
{
    $alphasg = $data->alphasg;
    $alphasj = $data->alphasj;
    $phi < 0 ? $alpha = $alphasg : $alpha = $alphasj;
    return $alpha;
}

function calc_Awx(&$data, $V)
{
    $A = $data->A;
    $Vc = $data->Vc;
    $D = $data->D;

    $Vx = $V - $Vc;
    $dx = $Vx / $A;
    $Awx = 2 * pi * $D * $dx;
    return $Awx;
}

function calc_dQw(&$data, $phi, $T, $V)
{
    $Twx = 135 + 273.15; # Температура цилиндра
    $Twp = 325 + 273.15; # Температура поршня
    $Twg = 325 + 273.15; # Температура головки

    $Awp = $data->A;
    $Awg = $data->A + $data->Vc / $data->A * 2 * pi * $data->D / 2;

    $Awx = calc_Awx($data, $V);
    $alpha = calc_alpha($data, $phi);
    $dQw = $alpha * ($Awp * ($Twp - $T) + $Awg * ($Twg - $T) + $Awx * ($Twx - $T));
    return $dQw;
}
# Формула Вибе
function calc_Wibe(&$data, $phi)
{
    $phiz = $data->phiz;
    $phisg = $data->phisg;
    $xz = $data->xz;
    $m = $data->mz;
    $omega = $data->omega;

    $betta = ($phi - $phisg);

    $c = log(1 - $xz);
    $x = 1 - exp($c * (($betta / $phiz) ** ($m + 1)));

    $dx = ((-((exp(($c * $betta ** ($m + 1)) / ($phiz ** ($m + 1))) * $c * $m * $betta ** $m + exp(($c *
        $betta ** ($m + 1)) / ($phiz ** ($m + 1))) * $c * $betta ** $m)) / ($phiz ** ($m + 1))) / (1 / $omega)) * 180 / pi;

    return $dx;
}

function calc_dQx(&$data, $phi)
{
    $Hu = $data->Hu;
    $phisg = $data->phisg;
    $phiz = $data->phiz;
    $mdt = $data->G_1;

    return $phi >= $phisg ? $Hu * $mdt * calc_Wibe($data, $phi) : 0;
}

function calc_cv(&$data, $T, $phi, $m)
{
    $phisg = $data->phisg;
    $phiz = $data->phiz;
    $alpha = $data->alpha;
    $xz = $data->xz;
    $mz = $data->mz;
    $gC = $data->gC;       # углерода в ДТ
    $gH = $data->gH;     # водорода в ДТ
    $gO = $data->gO;      # кислорода в ДТ
    $G_1 = $data->G_1;

    $ksi = ($phi - $phisg); # новый угол начинется с 0.
    $betta = $phisg + $phiz; # угол конца сгорания
    $M = $m / (28.97 * 1e-3); # кол-во воздуха в моль

    if ($phi <= $phisg) {

        $Mcv = 20.6 + 0.002638 * ($T - 273.15);
        $Cv = $Mcv * $M;
        $cv = $Cv / $m;
    } else {

        $c = log(1 - $xz);
        # закон впрыска топлива позакону Вибе
        $x = 1 - exp($c * (($ksi / ($phiz - 26)) ** ($mz + 1)));
        $mt = $G_1 * $x;
        $mall = $mt + $m; # общая масса горючей смеси, кг
        # масса СО2 равна массе С * 3.67, так как молярная масса в 3,67 раз больше и на одну молекулу углерода приходится 2 молекулы кислорода.
        $m_CO2 = $gC * $mt * 3.67;
        $m_H2O = $gH * $mt * 9;
        $m_O2 = $gO * $mt * 2;
        $g_O2 = $m_O2 / $mall;
        $g_H2O = $m_H2O / $mall;
        $g_CO2 = $m_CO2 / $mall;
        # допущенние этот воздух обедненый, то есть концентрация азота больше чем 70
        $g_VOZD = 1 - $g_CO2 - $g_H2O - $g_O2;

        $Mcv = 20.6 + 0.002638 * ($T - 273.15);
        $M = $m / (28.97 * 1e-3); # кол-во воздуха в моль
        $Cv = $Mcv * $M;
        $cv_VOZD = $Cv / $m;

        # для СО2
        $Mcv_CO2 = 27.941 + 0.019 * ($T - 273.15) - 0.000005487 * ($T - 273.15) * ($T - 273.15);
        $M_CO2 = $m_CO2 / (44.01 * 1e-3);  # кол-во CO2 в моль
        $Cv_CO2 = $Mcv_CO2 * $M_CO2;
        $cv_CO2 = $Cv_CO2 / $m_CO2;

        # для H2O
        $Mcv_H2O = 24.953 + 0.05359 * ($T - 273.15);
        $M_H2O = $m_H2O / (18 * 1e-3);  # кол-во CO2 в моль
        $Cv_H2O = $Mcv_H2O * $M_H2O;
        $cv_H2O = $Cv_H2O / $m_H2O;

        # для O2
        $Mcv_O2 = 20.930 + 0.00464 * ($T - 273.15) - 0.00000084 * ($T - 273.15) * ($T - 273.15);
        $M_O2 = $m_O2 / (32 * 1e-3); # кол-во CO2 в моль
        $Cv_O2 = $Mcv_O2 * $M_O2;
        $cv_O2 = $Cv_O2 / $m_O2;
        $cv = $g_CO2 * $cv_CO2 + $g_H2O * $cv_H2O + $g_O2 * $cv_O2 + $g_VOZD * $cv_VOZD;
    }
    return $cv;
}
