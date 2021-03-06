function [PVcurratio,PVPG,selfcon,fed] = f_DTR(PVratio, weather)

load CaseSwiss;
load Pload;

mpc.bus(2:94,13) = 0.95;

trans = zeros(18,1);
k = 1;
for i = 77:94
    for j = 1:109
        if mpc.bus(mpc.branch(j,1),7) == 5 && mpc.branch(j,2) == i
            trans(k,1) = j;
            k = k + 1;
        end
    end
end

delta_t = 60;
R = 8;
delta_Tto_R = 55;
delta_Th_R = 25;
Tao_to = 155;
Tao_h = 5; 
Cto = (1 - exp(-delta_t / Tao_to));
Ch = (1 - exp(-delta_t / Tao_h));

Ir = zeros(18,1); 
for i = 1:18
    Ir(i) = mpc.branch(trans(i),6) / 100;
end
Sr = zeros(18,1); 
for i = 1:18
    Sr(i) = mpc.branch(trans(i),6); % MVA
end

Ta_year = weather(:,2);
W_year = weather(:,1);
Tas_year = weather(:,4); 
k0 = 1.998 * 10^(-7);
a = 2.003 * 10^(-7);
w0 = 0.36;
b = 0.358;
Cv_year = zeros(8760,1);
ds_year = zeros(8760,1);
ps_year = zeros(8760,1);
w = weather(:,3);

for i = 1:8760
    Cv_year(i,1) = 837 * 1360 + 4.19 * 10^6 * w(i,1);
    ds_year(i,1) = k0 + a * exp(-0.5 * (log(w(i,1) / w0) / b)^2);
    ps_year(i,1) = 1 / Cv_year(i,1) / ds_year(i,1);
end

PVC = PVratio * mpc.bus(77:94,3);
Pz = zeros(94,8760);
Qz = zeros(94,8760);
for i = 77:94
    Pz(i,:) = mpc.bus(i,3) * Pload';
    Qz(i,:) = mpc.bus(i,4) * Pload';
end
PVPG = zeros(18,8760);
PV1 = zeros(8760,1);

Ths = zeros(18,8761);
delta_Tto = zeros(18,8761);
delta_Ths = zeros(18,8761);
for i = 1:18
    Ki = (Pz(mpc.branch(trans(i),2),8760)^2 + Qz(mpc.branch(trans(i),2),8760)^2)^0.5 / Sr(i);
    delta_Tto(i,1) = delta_Tto_R * (Ki^2 * R + 1) / (R + 1);
    delta_Ths(i,1) = delta_Th_R * Ki^2;
    Ths(i,1) = Ta_year(8760,1) + delta_Tto(i,1) + delta_Ths(i,1);
end

selfcon = zeros(96,8760);
fed = zeros(96,8760);

hd = zeros(4,8761);
hb = zeros(4,8761);
hs = zeros(4,8761);
hc = zeros(4,8761);
HC = zeros(4,8761); 

ZB = 0.42^2 / 100;

% cable parameter
Sc = 630;  
dc = 28.3; 
R0 = 0.0428 * 10^(-3);  
Ti = 2.4;
l1 = 0.042689413008005;
U0 = 420/3^0.5; 
f = 50;  
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
T4 = 1.5 * ps / pi * (log(2 * 800 / re) - 0.63);

Qc = 2.5 * Sc;
Qi = 2.4 * pi * (di^2 - dc^2) / 4;
Qs = 3.45 * pi * (dp^2 - di^2) / 4;
Qj = 1.7 * pi * (De^2 - dp^2) / 4;

re = De / 2;
rc = dc / 2;
pa = 1 / 2 / log(re / rc) - 1 / (re^2 / rc^2 - 1);
QA3 = Qc + pa * (Qi + (Qs + Qj) / (1 + l1));
TA3 = T1 + T3 * (1 + l1);

ktf = 3 * (log(2 * 800 / re) - 0.63) / log(800 / re + ((800 / re)^2 - 1)^0.5);
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

LV = [2; 12; 16; 18]; 

mpc.bus(:,3) = Pz(:,8760);
mpc.bus(:,4) = Qz(:,8760);
result = runpf(mpc);
for i = 1:4
    I = (result.branch(LV(i,1),14)^2+ result.branch(LV(i,1),15)^2)^0.5...
        * 1e3 / (result.bus(mpc.branch(LV(i,1),1),8)) / 3^0.5 / 0.42;
    Rac = R0 * (1 + ys + yp);
    Wc = I^2 * Rac; % W/m
    hd(i,1) = Wc * TD3;
    hb(i,1) = hd(i,1) + Wc * TC3;
    hs(i,1) = hb(i,1) + Wc * TB3;
    hc(i,1) = hs(i,1) + Wc * TA3;
    HC(i,1) = hc(i,1) + Tas_year(8760) + hd3;
end

for j = 1:8760
    % calculate maximum available PV output
    Tpv = Ta_year(j) + 0.03 * W_year(j);
    PV_ava = W_year(j) / 1000 .* (1 - 0.0045 * (Tpv - 25));
    PV1(j,1) = PV_ava;
    mpc.bus(:,3) = Pz(:,j);
    mpc.bus(:,4) = Qz(:,j);
    for i = 1:18
        mpc.gen(i+1,2) = PV_ava * PVC(i);
        PVPG(i,j) = PV_ava * PVC(i);
    end
    mpopt = mpoption('out.all', 0 ,'verbose', 0);
    [result,~] = runpf(mpc,mpopt);
    
    % calculate transformer temperature
    for i = 1:18
        Ki(i,j) = (result.branch(trans(i),16)^2 + result.branch(trans(i),17)^2)^0.5 / Sr(i);
        Ki2 = Ki(i,j)^2;
        delta_Tto(i,j+1) = (delta_Tto_R * (Ki2 * R + 1) / (R + 1) - delta_Tto(i,j)) * Cto + delta_Tto(i,j);
        delta_Ths(i,j+1) = (delta_Th_R * Ki2 - delta_Ths(i,j)) * Ch + delta_Ths(i,j);
        Ths(i,j+1) = Ta_year(j,1) + delta_Tto(i,j+1) + delta_Ths(i,j+1);
    end
    
    % calculate cable temperature
    for i = 1:4
        I = (result.branch(LV(i,1),14)^2+ result.branch(LV(i,1),15)^2)^0.5...
            * 1e3 / (result.bus(mpc.branch(LV(i,1),1),8)) / 3^0.5 / 0.42;
        Ht = cabletem758([hd(i,j);hb(i,j);hs(i,j);hc(i,j)], Tas_year(j), ps_year(j), ds_year(j), I);
        hd(i,j+1) = Ht(1,1);
        hb(i,j+1) = Ht(2,1);
        hs(i,j+1) = Ht(3,1);
        hc(i,j+1) = Ht(4,1);
        HC(i,j+1) = hc(i,j+1) + Tas_year(j) + hd3;
    end
    
    % If any constraint is violated, the optimal power flow is executed
    deltaI = result.branch(1:91,16).^2 + result.branch(1:91,17).^2 - result.branch(1:91,6).^2;
    if max(result.bus(:,8)) > 1.05 || max(Ths(:,j+1)) > 120 || max(HC(:,j+1)) > 90 || max(deltaI) > 0
        xPVcur = zeros(3,1);
        for k = 1:3
            nPVcur = [0.1 0.01 0.001];
            while true
                xPVcur(k,1) = xPVcur(k,1) + 1;
                mpc.gen(2:19,2) = (1 - nPVcur*xPVcur) * PV_ava * PVC;
                [result,~] = runpf(mpc,mpopt);
                % calculate transformer temperature
                for i = 1:18
                    Ki(i,j) = (result.branch(trans(i),16)^2 + result.branch(trans(i),17)^2)^0.5 / Sr(i);
                    Ki2 = Ki(i,j)^2;
                    delta_Tto(i,j+1) = (delta_Tto_R * (Ki2 * R + 1) / (R + 1) - delta_Tto(i,j)) * Cto + delta_Tto(i,j);
                    delta_Ths(i,j+1) = (delta_Th_R * Ki2 - delta_Ths(i,j)) * Ch + delta_Ths(i,j);
                    Ths(i,j+1) = Ta_year(j,1) + delta_Tto(i,j+1) + delta_Ths(i,j+1);
                end
    
                % calculate cable temperature
                for i = 1:4
                    I = (result.branch(LV(i,1),14)^2+ result.branch(LV(i,1),15)^2)^0.5...
                        * 1e3 / (result.bus(mpc.branch(LV(i,1),1),8)) / 3^0.5 / 0.42;
                    Ht = cabletem758([hd(i,j);hb(i,j);hs(i,j);hc(i,j)], Tas_year(j), ps_year(j), ds_year(j), I);
                    hd(i,j+1) = Ht(1,1);
                    hb(i,j+1) = Ht(2,1);
                    hs(i,j+1) = Ht(3,1);
                    hc(i,j+1) = Ht(4,1);
                    HC(i,j+1) = hc(i,j+1) + Tas_year(j) + hd3;
                end
                if max(result.bus(:,8)) < 1.05 && max(Ths(:,j+1)) < 120 && max(HC(:,j+1))
                    break
                end
            end
            xPVcur(k,1) = xPVcur(k,1) - 1;
        end
    end
    
    PVPG(:,j) = result.gen(2:19,2);

    for i = 1:18
        if Pz(mpc.branch(trans(i),2),j) > PVPG(i,j)
            selfcon(i,j) = PVPG(i,j);
        else
            selfcon(i,j) = Pz(mpc.branch(trans(i),2),j);
            fed(i,j) = PVPG(i,j) - Pz(mpc.branch(trans(i),2),j);
        end
    end
end

PV1sum = sum(PV1) * PVC;
PVcurratio = 1 - sum(PVPG,2) ./ PV1sum;

