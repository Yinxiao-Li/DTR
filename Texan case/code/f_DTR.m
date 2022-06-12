function [PVcurratio,PVPG,selfcon,fed] = f_DTR(PVratio, weather)

load CaseTexan;
load Pload;

Ta_year = weather(:,2);
W_year = weather(:,1);
Tas_year = weather(:,4); 

k0 = 1.346 * 10^(-7);
a = 2.883 * 10^(-7);
w0 = 0.419;
b = 0.302;
Cv_year = zeros(8760,1);
ds_year = zeros(8760,1);
ps_year = zeros(8760,1);
w = weather(:,3);

for i = 1:8760
    Cv_year(i,1) = 837 * 1520 + 4.19 * 10^6 * w(i,1);
    ds_year(i,1) = k0 + a * exp(-0.5 * (log(w(i,1) / w0) / b)^2);
    ps_year(i,1) = 1 / Cv_year(i,1) / ds_year(i,1);
end

mpc = mpc1;

Pz = zeros(53,8760);
for i = 30:48
    Pz(i,:) = mpc.bus(i,3) * Pload';
end
Qz = tan(acos(0.98)) * Pz;

PVC = PVratio * mpc.bus(30:48,3);

mpc.gen(2:20,:) = zeros(19,21);
mpc.gen(2:20,1) = (30:48)';
mpc.gen(:,7) = 100;
mpc.gen(:,8) = 1;

mpc.gencost(2:20,:) = zeros(19,7);
mpc.gencost(2:20,4) = 3;
mpc.gencost(2:20,6) = -1;
mpc.gencost(2:20,1) = 2;

% transformer parameters
delta_t = 60;
R = 8;
delta_Tto_R = 55;
delta_Th_R = 25;
Tao_to = 155; 
Tao_h = 5; 
Cto = (1 - exp(-delta_t / Tao_to));
Ch = (1 - exp(-delta_t / Tao_h));

Ir = zeros(14,1); 
for i = 1:14
    Ir(i) = mpc.branch(i+38,6) / 100;
end
Sr = zeros(14,1); 
for i = 1:14
    Sr(i) = mpc.branch(i+38,6);
end

Ths = zeros(14,8761);
delta_Tto = zeros(14,8761);
delta_Ths = zeros(14,8761);
mpc.bus(:,3) = Pz(:,8760);
mpc.bus(:,4) = Qz(:,8760);
result = runpf(mpc);
for i = 1:14
    Ki = (result.branch(i+38,16)^2 + result.branch(i+38,17)^2)^0.5 / Sr(i);
    delta_Tto(i,1) = delta_Tto_R * (Ki^2 * R + 1) / (R + 1);
    delta_Ths(i,1) = delta_Th_R * Ki^2;
    Ths(i,1) = Ta_year(8760,1) + delta_Tto(i,1) + delta_Ths(i,1);
end

selfcon = zeros(length(PVC),8760);
fed = zeros(length(PVC),8760);
PVPG = zeros(19,8760);
PV1 = zeros(8760,1);

MVline1 = [2;3;4;5;6;7;17;28;29;32;33;36;37];
MVline2 = [1;38];
MVI1 = zeros(8760,1);
MVI2 = zeros(8760,1);

LV370 = find(mpc.branch(1:38,6) < 0.31);
LV550 = zeros(1,1);
k = 1;
for i = 1:38
    if mpc.branch(i,6) > 0.31 && mpc.branch(i,6) < 0.5
        LV550(k,1) = i;
        k = k + 1;
    end
end

hd = zeros(23,8761);
hb = zeros(23,8761);
hs = zeros(23,8761);
hc = zeros(23,8761);
HC = zeros(23,8761); 
ZB = 0.48^2 / 100;

%% cable 370
Sc = 185; 
dc = 15.3;  
R0 = 0.164 * 10^(-3);  
Ti = 1.6;
l1 = 0.010832237335199;
U0 = 480/3^0.5;  
f = 60; 
ps = ps_year(8760);
ds = ds_year(8760);

di = dc + 2 * Ti; 
dp = di + 0.2 * 2;  
Tm = 0.035 * dp + 1; 
De = dp + 2 * Tm; 

Rm = R0 * (1 + 0.00403 * (90 - 20)); 
xs2 = 8 * pi * f * 10^(-7) / Rm;
ys = xs2^2 / (192 + 0.8 * xs2^2);
xp2 = 8 * pi * f * 10^(-7) / Rm;  
s = De;  
yp = xp2^2 / (192 + 0.8 * xp2^2) * (dc / s)^2 * (0.312 * (dc / s)^2 + 1.18 / (xp2^2 / (192 + 0.8 * xp2^2) + 0.27));

T1 = 3.5 / 2 / pi * log(di / dc);
T3 = 6 / 2 / pi * log(De / dp);
re = De / 2;
T4 = 1.5 * ps / pi * (log(2 * 900 / re) - 0.63);

Qc = 2.5 * Sc;
Qi = 2.4 * pi * (di^2 - dc^2) / 4;
Qs = 3.45 * pi * (dp^2 - di^2) / 4;
Qj = 1.7 * pi * (De^2 - dp^2) / 4;

re = De / 2;
rc = dc / 2;
pa = 1 / 2 / log(re / rc) - 1 / (re^2 / rc^2 - 1);
QA3 = Qc + pa * (Qi + (Qs + Qj) / (1 + l1));
TA3 = T1 + T3 * (1 + l1);

ktf = 3 * (log(2 * 900 / re) - 0.63) / log(900 / re + ((900 / re)^2 - 1)^0.5);
rb = 3 * re;
pb = 1 / 2 / log(rb / re) - 1 / (rb^2 / re^2 - 1);
TB3 = ktf * ps * (1 + l1) / 2 / pi * log(rb / re);
QB3 = (1 - pa) * (Qi + (Qs + Qj) / (1 + l1)) + pb * pi * pi * (rb^2 - re^2) * 10^(-6) / ktf / ps / ds / (1 + l1);

rc = 3 * rb;
pc = 1 / 2 / log(rc / rb) - 1 / (rc^2 / rb^2 - 1);
TC3 = ktf * ps * (1 + l1) / 2 / pi * log(rc / rb);
QC3 = (1 - pb) * pi * (rb^2 - re^2) * 10^(-6) / ktf / ps / ds / (1 + l1) + pc * pi * (rc^2 - rb^2) * 10^(-6) / ktf / ps / ds / (1 + l1);

renv = re * exp(2 * pi * T4 / ktf / ps);
penv = 1 / 2 / log(renv / rc) - 1 / (renv^2 / rc^2 - 1);
TD3 = T4 * (1 + l1) - (TB3 + TC3);
QD3 = (1 - pc) * pi * (rc^2 - rb^2) * 10^(-6) / ktf / ps / ds / (1 + l1) + penv * pi * (renv^2 - rc^2) * 10^(-6) / ktf / ps / ds / (1 + l1);

c = 2.5 / 18 / log(di / dc) * 10^(-9); 
Wd = 2 * pi * f * c * U0^2 * 0.004;  
hd3 = Wd * (0.5 * T1 + T3 + T4);

for i = 1:9
    I = (result.branch(LV370(i,1),14)^2+ result.branch(LV370(i,1),15)^2)^0.5...
         * 1e3 / (result.bus(mpc.branch(LV370(i,1),1),8)) / 3^0.5 / 0.48;
    R370 = R0 * (1 + ys + yp);
    Wc = I^2 * R370; 
    hd(i,1) = Wc * TD3;
    hb(i,1) = hd(i,1) + Wc * TC3;
    hs(i,1) = hb(i,1) + Wc * TB3;
    hc(i,1) = hs(i,1) + Wc * TA3;
    HC(i,1) = hc(i,1) + Tas_year(8760) + hd3;
end

%% cable 550

Sc = 400;  
dc = 22.6; 
R0 = 0.0778 * 10^(-3);  
Ti = 2;
l1 = 0.029491750817908;
U0 = 480/3^0.5; 
f = 60;  

di = dc + 2 * Ti; 
dp = di + 0.2 * 2;  
Tm = 0.035 * dp + 1;  
De = dp + 2 * Tm; 

Rm = R0 * (1 + 0.00403 * (90 - 20));  
xs2 = 8 * pi * f * 10^(-7) / Rm;
ys = xs2^2 / (192 + 0.8 * xs2^2);
xp2 = 8 * pi * f * 10^(-7) / Rm;  
s = De; 
yp = xp2^2 / (192 + 0.8 * xp2^2) * (dc / s)^2 * (0.312 * (dc / s)^2 + 1.18 / (xp2^2 / (192 + 0.8 * xp2^2) + 0.27));

T1 = 3.5 / 2 / pi * log(di / dc);
T3 = 6 / 2 / pi * log(De / dp);
re = De / 2;
T4 = 1.5 * ps / pi * (log(2 * 900 / re) - 0.63);

Qc = 2.5 * Sc;
Qi = 2.4 * pi * (di^2 - dc^2) / 4;
Qs = 3.45 * pi * (dp^2 - di^2) / 4;
Qj = 1.7 * pi * (De^2 - dp^2) / 4;

re = De / 2;
rc = dc / 2;
pa = 1 / 2 / log(re / rc) - 1 / (re^2 / rc^2 - 1);
QA5 = Qc + pa * (Qi + (Qs + Qj) / (1 + l1));
TA5 = T1 + T3 * (1 + l1);

ktf = 3 * (log(2 * 900 / re) - 0.63) / log(900 / re + ((900 / re)^2 - 1)^0.5);
rb = 3 * re;
pb = 1 / 2 / log(rb / re) - 1 / (rb^2 / re^2 - 1);
TB5 = ktf * ps * (1 + l1) / 2 / pi * log(rb / re);
QB5 = (1 - pa) * (Qi + (Qs + Qj) / (1 + l1)) + pb * pi * pi * (rb^2 - re^2) * 10^(-6) / ktf / ps / ds / (1 + l1);

rc = 3 * rb;
pc = 1 / 2 / log(rc / rb) - 1 / (rc^2 / rb^2 - 1);
TC5 = ktf * ps * (1 + l1) / 2 / pi * log(rc / rb);
QC5 = (1 - pb) * pi * (rb^2 - re^2) * 10^(-6) / ktf / ps / ds / (1 + l1) + pc * pi * (rc^2 - rb^2) * 10^(-6) / ktf / ps / ds / (1 + l1);

renv = re * exp(2 * pi * T4 / ktf / ps);
penv = 1 / 2 / log(renv / rc) - 1 / (renv^2 / rc^2 - 1);
TD5 = T4 * (1 + l1) - (TB5 + TC5);
QD5 = (1 - pc) * pi * (rc^2 - rb^2) * 10^(-6) / ktf / ps / ds / (1 + l1) + penv * pi * (renv^2 - rc^2) * 10^(-6) / ktf / ps / ds / (1 + l1);

c = 2.5 / 18 / log(di / dc) * 10^(-9);  % 单位长度电缆电容
Wd = 2 * pi * f * c * U0^2 * 0.004;  % 介质损耗
hd5 = Wd * (0.5 * T1 + T3 + T4);

for i = 1:14
    I = (result.branch(LV550(i,1),14)^2+ result.branch(LV550(i,1),15)^2)^0.5...
        * 1e3 / (result.bus(mpc.branch(LV550(i,1),1),8)) / 3^0.5 / 0.48;
    R550 = R0 * (1 + ys + yp);
    Wc = I^2 * R550; % W/m
    hb(i+9,1) = hd(i+9,1) + Wc * TC5;
    hs(i+9,1) = hb(i+9,1) + Wc * TB5;
    hc(i+9,1) = hs(i+9,1) + Wc * TA5;
    HC(i+9,1) = hc(i+9,1) + Tas_year(8760) + hd5;
end

for j = 1:8760
    % calculate maximum available PV output
    Tpv = Ta_year(j) + 0.03 * W_year(j);
    PV_ava = W_year(j) / 1000 * (1 - 0.0045 * (Tpv - 25));
    PV1(j,1) = PV_ava;
    mpc.bus(:,3) = Pz(:,j);
    mpc.bus(:,4) = Qz(:,j);
    for i = 1:19
        mpc.gen(i+1,2) = PV_ava * PVC(i);
        PVPG(i,j) = PV_ava * PVC(i);
    end
    % The power flow calculation is executed to judge whether 
    % the equipment capacity constraint and voltage constraint are violated
    mpopt = mpoption('out.all', 0 ,'verbose', 0);
    [result,~] = runpf(mpc,mpopt);
    MVI1(j,1) = max(abs((result.branch(MVline1,16).^2 + result.branch(MVline1,17).^2).^0.5));
    MVI2(j,1) = max(abs((result.branch(MVline2,16).^2 + result.branch(MVline2,17).^2).^0.5));
    
    % calculate transformer temperature
    for i = 1:14
        Ki(i,j) = (result.branch(i+38,16)^2 + result.branch(i+38,17)^2)^0.5 / Sr(i);
        Ki2 = Ki(i,j)^2;
        delta_Tto(i,j+1) = (delta_Tto_R * (Ki2 * R + 1) / (R + 1) - delta_Tto(i,j)) * Cto + delta_Tto(i,j);
        delta_Ths(i,j+1) = (delta_Th_R * Ki2 - delta_Ths(i,j)) * Ch + delta_Ths(i,j);
        Ths(i,j+1) = Ta_year(j,1) + delta_Tto(i,j+1) + delta_Ths(i,j+1);
    end
    % calculate cable temperature
    for i = 1:9
        I = (result.branch(LV370(i,1),14)^2+ result.branch(LV370(i,1),15)^2)^0.5...
            * 1e3 / (result.bus(mpc.branch(LV370(i,1),1),8)) / 3^0.5 / 0.48;
        Ht = cabletem370([hd(i,j);hb(i,j);hs(i,j);hc(i,j)], Tas_year(j), ps_year(j), ds_year(j), I);
        hd(i,j+1) = Ht(1,1);
        hb(i,j+1) = Ht(2,1);
        hs(i,j+1) = Ht(3,1);
        hc(i,j+1) = Ht(4,1);
        HC(i,j+1) = hc(i,j+1) + Tas_year(j) + hd3;
    end
    
    for i = 1:14
        I = (result.branch(LV550(i,1),14)^2+ result.branch(LV550(i,1),15)^2)^0.5...
        * 1e3 / (result.bus(mpc.branch(LV550(i,1),1),8)) / 3^0.5 / 0.48;
        Ht = cabletem550([hd(i+9,j);hb(i+9,j);hs(i+9,j);hc(i+9,j)], Tas_year(j), ps_year(j), ds_year(j), I);
        hd(i+9,j+1) = Ht(1,1);
        hb(i+9,j+1) = Ht(2,1);
        hs(i+9,j+1) = Ht(3,1);
        hc(i+9,j+1) = Ht(4,1);
        HC(i+9,j+1) = hc(i+9,j+1) + Tas_year(j) + hd5;
    end
    
    % If any constraint is violated, calculate the curtailed PV generation
    if max(result.bus(:,8)) > 1.05 || max(Ths(:,j+1)) > 120 || max(HC(:,j+1)) > 90 || MVI1(j,1) > mpc.branch(2,6) || MVI2(j,1) > mpc.branch(1,6)     
        xPVcur = zeros(3,1);
        for k = 1:3
            nPVcur = [0.1 0.01 0.001];
            while true
                xPVcur(k,1) = xPVcur(k,1) + 1;
                mpc.gen(2:20,2) = (1 - nPVcur*xPVcur) * PV_ava * PVC;
                [result,~] = runpf(mpc,mpopt);
                % calculate transformer temperature
                for i = 1:14
                    Ki(i,j) = (result.branch(i+38,16)^2 + result.branch(i+38,17)^2)^0.5 / Sr(i);
                    Ki2 = Ki(i,j)^2;
                    delta_Tto(i,j+1) = (delta_Tto_R * (Ki2 * R + 1) / (R + 1) - delta_Tto(i,j)) * Cto + delta_Tto(i,j);
                    delta_Ths(i,j+1) = (delta_Th_R * Ki2 - delta_Ths(i,j)) * Ch + delta_Ths(i,j);
                    Ths(i,j+1) = Ta_year(j,1) + delta_Tto(i,j+1) + delta_Ths(i,j+1);
                end
                % calculate cable temperature
                for i = 1:9
                    I = (result.branch(LV370(i,1),14)^2+ result.branch(LV370(i,1),15)^2)^0.5...
                        * 1e3 / (result.bus(mpc.branch(LV370(i,1),1),8)) / 3^0.5 / 0.48;
                    Ht = cabletem370([hd(i,j);hb(i,j);hs(i,j);hc(i,j)], Tas_year(j), ps_year(j), ds_year(j), I);
                    hd(i,j+1) = Ht(1,1);
                    hb(i,j+1) = Ht(2,1);
                    hs(i,j+1) = Ht(3,1);
                    hc(i,j+1) = Ht(4,1);
                    HC(i,j+1) = hc(i,j+1) + Tas_year(j) + hd3;
                end
    
                for i = 1:14
                    I = (result.branch(LV550(i,1),14)^2+ result.branch(LV550(i,1),15)^2)^0.5...
                        * 1e3 / (result.bus(mpc.branch(LV550(i,1),1),8)) / 3^0.5 / 0.48;
                    Ht = cabletem550([hd(i+9,j);hb(i+9,j);hs(i+9,j);hc(i+9,j)], Tas_year(j), ps_year(j), ds_year(j), I);
                    hd(i+9,j+1) = Ht(1,1);
                    hb(i+9,j+1) = Ht(2,1);
                    hs(i+9,j+1) = Ht(3,1);
                    hc(i+9,j+1) = Ht(4,1);
                    HC(i+9,j+1) = hc(i+9,j+1) + Tas_year(j) + hd5;
                end
                if max(result.bus(:,8)) < 1.05 && max(Ths(:,j+1)) < 120 && max(HC(:,j+1)) < 90
                    break
                end
            end
            xPVcur(k,1) = xPVcur(k,1) - 1;
        end
    end
    PVPG(:,j) = result.gen(2:20,2);
    for i = 1:length(PVC)
        if Pz(find(mpc.bus(:,1)==mpc.gen(i+1,1)),j) > PVPG(i,j)
            selfcon(i,j) = PVPG(i,j);
        else
            selfcon(i,j) = Pz(find(mpc.bus(:,1)==mpc.gen(i+1,1)),j);
            fed(i,j) = PVPG(i,j) - Pz(find(mpc.bus(:,1)==mpc.gen(i+1,1)),j);
        end
    end
end

PV1sum = sum(PV1) * PVC;
PVcurratio = 1 - sum(PVPG,2) ./ PV1sum;
