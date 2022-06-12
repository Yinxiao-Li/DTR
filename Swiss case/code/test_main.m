weather = xlsread('weather_2020.xlsx');
load CaseSwiss;
load Pload;

Ta_year = weather(1:8760,2);
W_year = weather(1:8760,1);
Tpv = Ta_year + 0.03 * W_year;
PV_ava = W_year / 1000 .* (1 - 0.0045 * (Tpv - 25));

CostInv0 = zeros(600,1);
Subsidy0 = zeros(600,1);
selfcon0 = zeros(600,1);
PVR = zeros(600,1);
for k = 1:600
    PVR(k,1) = k/100;
    CostInv0(k,1) = PVR(k,1) * 2313000;
    Subsidy0(k,1) = 0.3 * CostInv0(k,1);
    selfcon0(k,1) = sum(min([PV_ava * PVR(k,1) - Pload, zeros(8760,1)]')' + Pload);
end
fed0 = sum(PV_ava) * PVR - selfcon0;
year = 1:25;
Revenue0_1year = 223.4 * selfcon0 + 120 * fed0;
Revenue0 = sum(1 ./ (1 + 0.05).^year) * Revenue0_1year;
CostOM0 = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv0;
NPV0 = Revenue0 + Subsidy0 - CostInv0 - CostOM0;
PVratio = find(NPV0==max(NPV0))/100;

% calculate the PV curtailment ratio with STR
[PVcurSTR,PVSTR,selfSTR,fedSTR] = f_STR(PVratio, weather);
PVratioSTR = PVratio;
if max(PVcurSTR) > 0.03
    xSTR = zeros(4,1);
    for k = 1:4
        nSTR = [0.1 0.01 0.001 0.0001];
        while true
            xSTR(k,1) = xSTR(k,1) + 1;
            PVratioSTR = (1 - nSTR*xSTR) * PVratio;
            [PVcurSTR,PVSTR,selfSTR,fedSTR] = f_STR(PVratioSTR, weather);
            if max(PVcurSTR) < 0.03
                break
            end
        end
        xSTR(k,1) = xSTR(k,1) - 1;
    end
    PVratioSTR = (1 - nSTR*xSTR - 0.0001) * PVratio;
end

% calculate the PV curtailment ratio with DTR
[PVcurDTR,PVDTR,selfDTR,fedDTR] = f_DTR(PVratio, weather);
PVratioDTR = PVratio;
if max(PVcurDTR) > 0.03
    xDTR = zeros(4,1);
    for k = 1:4
        nDTR = [0.1 0.01 0.001 0.0001];
        while true
            xDTR(k,1) = xDTR(k,1) + 1;
            PVratioDTR = (1 - nDTR*xDTR) * PVratio;
            [PVcurDTR,PVDTR,selfDTR,fedDTR] = f_DTR(PVratioDTR, weather);
            if max(PVcurDTR) < 0.03
                break
            end
        end
        xDTR(k,1) = xDTR(k,1) - 1;
    end
    PVratioDTR = (1 - nDTR*xDTR - 0.0001) * PVratio;
end

load CaseSwiss;
Revenue_1year = 223.4 * [sum(sum(selfSTR)) sum(sum(selfDTR))] + ...
                120 * [sum(sum(fedSTR)) sum(sum(fedDTR))];
CostInv = [PVratioSTR PVratioDTR] * sum(mpc.bus(:,3)) * 2313000;
Subsidy = 0.3 * CostInv;
year = 1:25;
Revenue = sum(1 ./ (1 + 0.05).^year) * Revenue_1year;
CostOM = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv;
NPV = Revenue + Subsidy - CostInv - CostOM;
