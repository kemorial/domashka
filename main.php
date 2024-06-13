<?php

require_once ('./calc_pack.php');

$k = $data->k;
$a = ceil((360 - 50) / $k);

$phi = array_fill(0, $a, (float) 0);
$p = array_fill(0, $a, (float) 0);
$V = array_fill(0, $a, (float) 0);
$T = array_fill(0, $a, (float) 0);
$dL = array_fill(0, $a, (float) 0);
$dT = array_fill(0, $a, (float) 0);
$dQw = array_fill(0, $a, (float) 0);
$dQx = array_fill(0, $a, (float) 0);
$cv = array_fill(0, $a, (float) 0);
$phi1 = array_fill(0, $a, (float) 0);
$p1 = array_fill(0, $a, (float) 0);
$V1 = array_fill(0, $a, (float) 0);
$T1 = array_fill(0, $a, (float) 0);
$dL1 = array_fill(0, $a, (float) 0);
$dT1 = array_fill(0, $a, (float) 0);
$dQw1 = array_fill(0, $a, (float) 0);
$cv1 = array_fill(0, $a, (float) 0);


$phi[0] = $data->phi0;
$p[0] = $data->p1;
$T[0] = $data->T1;
$V[0] = calc_Volume($data, $phi[0]);
$R = $data->R;
$m = calc_m($data, $V[0]);
$G_1 = calc_G_1($data, $m);
$dt = $data->deltat;

$L = 0;
$shk = [];
for ($i = 0; $i < $a; $i += $k) {
    $shk[] = $i;
}
for ($i = 0; $i < $a; $i += $k) {
    $dH = 0;
    $dU = 0;
    $cv[$i] = calc_cv($data, $T[$i], $phi[$i], $m);
    $dL[$i] = calc_dL($data, $p[$i], $phi[$i]);
    $dQw[$i] = calc_dQw($data, $phi[$i], $T[$i], $V[$i]);
    $dQx[$i] = calc_dQx($data, $phi[$i]);

    $dT[$i] = ($dQw[$i] + $dQx[$i] + $dH - $dU - $dL[$i]) / $m / $cv[$i];

    $phi[$i + 1] = $phi[$i] + $k;
    $V[$i + 1] = calc_Volume($data, $phi[$i + 1]);
    $T[$i + 1] = $T[$i] + $dT[$i] * $dt;
    $p[$i + 1] = $m * $R * $T[$i + 1] / $V[$i + 1];
    $dL[$i] = ($V[$i + 1] - $V[$i]) * ($p[$i] + $p[$i + 1]) / 2;
    $L += $dL[$i];
}


$Vh = $data->Vh;
$Hu = $data->Hu;
$i = $data->i;
$n = $data->n;
$omega = $data->omega;

$tau = 4;# Четыре такта     
$p_ind = $L / $Vh;# Индикаторное давление
$N_ind = $p_ind * $i * $Vh * $n / 60 * 2 / $tau;   # Индикаторная мощность    
$f = 3600 / (720 * pi / 180 / $omega);
$G = ($G_1) * $f;           # Подача топлива
$g_ind = $G * $i / $N_ind;            # Индикаторная подача топлива    
$ef_ind = 3600 / $Hu / $g_ind;           # Индикаторный КПД    

# Проведём расчёт параметров без сгорания топлива 

$phi1[0] = $data->phi0;
$p1[0] = $data->p1;
$T1[0] = $data->T1;
$V1[0] = calc_Volume($data, $phi1[0]);

$R = $data->R;
$m = calc_m($data, $V1[0]);
$dt = $data->deltat;
print ("L = $L" . PHP_EOL);
$L = 0;
for ($i = 0; $i < $a; $i++) {
    $dH = 0;
    $dU = 0;
    $cv1[$i] = calc_cv($data, $T1[$i], $phi1[$i], $m);
    $dL1[$i] = calc_dL($data, $p1[$i], $phi1[$i]);
    $dQw1[$i] = calc_dQw($data, $phi1[$i], $T1[$i], $V1[$i]);
    $dQx1 = 0;
    $dT1[$i] = ($dQw1[$i] + $dQx1 + $dH - $dU - $dL1[$i]) / $m / $cv1[$i];
    $phi1[$i + 1] = $phi1[$i] + $k;
    $V1[$i + 1] = calc_Volume($data, $phi1[$i + 1]);

    $T1[$i + 1] = $T1[$i] + $dT1[$i] * $dt;
    $p1[$i + 1] = $m * $R * $T1[$i + 1] / $V1[$i + 1];
}


# Вывод параметров

print ("p_ind = $p_ind" . PHP_EOL . "N_ind = $N_ind" . PHP_EOL .
    "g_ind = $g_ind" . PHP_EOL . "ef_ind = $ef_ind");
