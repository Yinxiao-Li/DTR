weather = xlsread('weather_2020.xlsx');
load Pload;

%% Swiss policy

Ta_year = weather(1:8760,2);
W_year = weather(1:8760,1);
Tpv = Ta_year + 0.03 * W_year;
PV_ava = W_year / 1000 .* (1 - 0.0045 * (Tpv - 25));

CostInv0 = zeros(500,1);
Subsidy0 = zeros(500,1);
selfcon0 = zeros(500,1);
PVR = zeros(500,1);
for k = 1:500
    PVR(k,1) = k/100;
    CostInv0(k,1) = PVR(k,1) * 2675000;
    Subsidy0(k,1) = 0.3 * CostInv0(k,1);
    selfcon0(k,1) = sum(min([PV_ava * PVR(k,1) - Pload, zeros(8760,1)]')' + Pload);
end
fed0 = sum(PV_ava) * PVR - selfcon0;
year = 1:25;
Revenue0_1year = 124.2 * selfcon0 + 124.2*120/223.4 * fed0;
Revenue0 = sum(1 ./ (1 + 0.05).^year) * Revenue0_1year;
CostOM0 = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv0;
NPV0 = Revenue0 + Subsidy0 - CostInv0 - CostOM0;
PVratioSwiss = find(NPV0==max(NPV0))/100;

% calculate the PV curtailment ratio with STR
[PVcurSTRSwiss,PVSTRSwiss,selfSTRSwiss,fedSTRSwiss] = f_STR(PVratioSwiss, weather);
PVratioSTRSwiss = PVratioSwiss;

if max(PVcurSTRSwiss) > 0.03
    xSTRSwiss = zeros(4,1);
    for k = 1:4
        nSTR = [0.1 0.01 0.001 0.0001];
        while true
            xSTRSwiss(k,1) = xSTRSwiss(k,1) + 1;
            PVratioSTRSwiss = (1 - nSTR*xSTRSwiss) * PVratioSwiss;
            [PVcurSTRSwiss,PVSTRSwiss,selfSTRSwiss,fedSTRSwiss] = f_STR(PVratioSTRSwiss, weather);
            if max(PVcurSTRSwiss) < 0.03
                break
            end
        end
        xSTRSwiss(k,1) = xSTRSwiss(k,1) - 1;
    end
    PVratioSTRSwiss = (1 - nSTR*xSTRSwiss - 0.0001) * PVratioSwiss;
end

% calculate the PV curtailment ratio with DTR
[PVcurDTRSwiss,PVDTRSwiss,selfDTRSwiss,fedDTRSwiss] = f_DTR(PVratioSwiss, weather);
PVratioDTRSwiss = PVratioSwiss;

if max(PVcurDTRSwiss) > 0.03
    xDTRSwiss = zeros(4,1);
    for k = 1:4
        nDTR = [0.1 0.01 0.001 0.0001];
        while true
            xDTRSwiss(k,1) = xDTRSwiss(k,1) + 1;
            PVratioDTRSwiss = (1 - nDTR*xDTRSwiss) * PVratioSwiss;
            [PVcurDTRSwiss,PVDTRSwiss,selfDTRSwiss,fedDTRSwiss] = f_DTR(PVratioDTRSwiss, weather);
            if max(PVcurDTRSwiss) < 0.03
                break
            end
        end
        xDTRSwiss(k,1) = xDTRSwiss(k,1) - 1;
    end
    PVratioDTRSwiss = (1 - nDTR*xDTRSwiss - 0.0001) * PVratioSwiss;
end

load CaseTexan;
mpc = mpc1;
Revenue_1year = 124.2 * [sum(sum(selfSTRSwiss)) sum(sum(selfDTRSwiss))] + ...
                124.2*120/223.4 * [sum(sum(fedSTRSwiss)) sum(sum(fedDTRSwiss))];
CostInv = [PVratioSTRSwiss PVratioDTRSwiss] * sum(mpc.bus(:,3)) * 2675000;
Subsidy = 0.3 * CostInv;
year = 1:25;
Revenue = sum(1 ./ (1 + 0.05).^year) * Revenue_1year;
CostOM = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv;
NPVSwiss = Revenue + Subsidy - CostInv - CostOM;


%% Chinese policy

CostInv0 = zeros(500,1);
selfcon0 = zeros(500,1);
PVR = zeros(500,1);
for k = 1:500
    PVR(k,1) = k/100;
    CostInv0(k,1) = PVR(k,1) * 2675000;
    selfcon0(k,1) = sum(min([PV_ava * PVR(k,1) - Pload, zeros(8760,1)]')' + Pload);
end
fed0 = sum(PV_ava) * PVR - selfcon0;
year = 1:25;
Revenue0_1year = 124.2 * selfcon0 + 22 * fed0;
Subsidy0_1year = 80/565.3*124.2 * sum(PV_ava) * PVR;

Revenue0 = sum(1 ./ (1 + 0.05).^year) * Revenue0_1year;
Subsidy0 = sum(1 ./ (1 + 0.05).^year) * Subsidy0_1year;
CostOM0 = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv0;
NPV0 = Revenue0 + Subsidy0 - CostInv0 - CostOM0;
PVratioChinese = find(NPV0==max(NPV0))/100;

[PVcurSTRChinese,PVSTRChinese,selfSTRChinese,fedSTRChinese] = f_STR(PVratioChinese, weather);
PVratioSTRChinese = PVratioChinese;

if max(PVcurSTRChinese) > 0.03
    xSTRChinese = zeros(4,1);
    for k = 1:4
        nSTR = [0.1 0.01 0.001 0.0001];
        while true
            xSTRChinese(k,1) = xSTRChinese(k,1) + 1;
            PVratioSTRChinese = (1 - nSTR*xSTRChinese) * PVratioChinese;
            [PVcurSTRChinese,PVSTRChinese,selfSTRChinese,fedSTRChinese] = f_STR(PVratioSTRChinese, weather);
            if max(PVcurSTRChinese) < 0.03
                break
            end
        end
        xSTRChinese(k,1) = xSTRChinese(k,1) - 1;
    end
    PVratioSTRChinese = (1 - nSTR*xSTRChinese - 0.0001) * PVratioChinese;
end

% calculate the PV curtailment ratio with DTR
[PVcurDTRChinese,PVDTRChinese,selfDTRChinese,fedDTRChinese] = f_DTR(PVratioChinese, weather);
PVratioDTRChinese = PVratioChinese;

if max(PVcurDTRChinese) > 0.03
    xDTRChinese = zeros(4,1);
    for k = 1:4
        nDTR = [0.1 0.01 0.001 0.0001];
        while true
            xDTRChinese(k,1) = xDTRChinese(k,1) + 1;
            PVratioDTRChinese = (1 - nDTR*xDTRChinese) * PVratioChinese;
            [PVcurDTRChinese,PVDTRChinese,selfDTRChinese,fedDTRChinese] = f_DTR(PVratioDTRChinese, weather);
            if max(PVcurDTRChinese) < 0.03
                break
            end
        end
        xDTRChinese(k,1) = xDTRChinese(k,1) - 1;
    end
    PVratioDTRChinese = (1 - nDTR*xDTRChinese - 0.0001) * PVratioChinese;
end

Revenue_1year = 124.2 * [sum(sum(selfSTRChinese)) sum(sum(selfDTRChinese))] ...
                + 22 * [sum(sum(fedSTRChinese)) sum(sum(fedDTRChinese))];
Subsidy_1year = 80/565.3*124.2 * [sum(sum(PVSTRChinese)) sum(sum(PVDTRChinese))];
CostInv = [PVratioSTRChinese PVratioDTRChinese] * sum(mpc.bus(:,3)) * 2675000;

year = 1:25;
Revenue = sum(1 ./ (1 + 0.05).^year) * Revenue_1year;
CostOM = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv;
Subsidy = sum(1 ./ (1 + 0.05).^year) * Subsidy_1year;
NPVChinese = Revenue + Subsidy - CostInv - CostOM;
