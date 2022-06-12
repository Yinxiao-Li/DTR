weather = xlsread('weather_2020.xlsx');
load CaseSwiss;
load Pload;

Ta_year = weather(1:8760,2);
W_year = weather(1:8760,1);
Tpv = Ta_year + 0.03 * W_year;
PV_ava = W_year / 1000 .* (1 - 0.0045 * (Tpv - 25));

%% Texan policy

PVR = (0.01:0.01:5); % MW
CostInv0 = PVR * 2313000;
Subsidy0 = 0.26 * CostInv0 + PVR / 0.005 * 2500*2313000/2675000;
Ploadyear = sum(Pload);
PVG = sum(PV_ava) * PVR;
Revenue0_1year = min([PVG * 97/124.2*223.4 - Ploadyear * 223.4; ones(1,length(PVR))]) + Ploadyear * 223.4;
year = 1:25;
Revenue0 = sum(1 ./ (1 + 0.05).^year) * Revenue0_1year;
CostOM0 = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv0;
NPV0 = Revenue0 + Subsidy0 - CostInv0 - CostOM0;

PVratioTexan = find(NPV0 == max(NPV0)) / 100;

% calculate the PV curtailment ratio with STR
[PVcurSTRTexan,PVSTRTexan,selfSTRTexan,fedSTRTexan] = f_STR(PVratioTexan, weather);
PVratioSTRTexan = PVratioTexan;
if max(PVcurSTRTexan) > 0.03
    xSTRTexan = zeros(4,1);
    for k = 1:4
        nSTR = [0.1 0.01 0.001 0.0001];
        while true
            xSTRTexan(k,1) = xSTRTexan(k,1) + 1;
            PVratioSTRTexan = (1 - nSTR*xSTRTexan) * PVratioTexan;
            [PVcurSTRTexan,PVSTRTexan,selfSTRTexan,fedSTRTexan] = f_STR(PVratioSTRTexan, weather);
            if max(PVcurSTRTexan) < 0.03
                break
            end
        end
        xSTRTexan(k,1) = xSTRTexan(k,1) - 1;
    end
    PVratioSTRTexan = (1 - nSTR*xSTRTexan - 0.0001) * PVratioTexan;
end

% calculate the PV curtailment ratio with DTR
[PVcurDTRTexan,PVDTRTexan,selfDTRTexan,fedDTRTexan] = f_DTR(PVratioTexan, weather);
PVratioDTRTexan = PVratioTexan;
if max(PVcurDTRTexan) > 0.03
    xDTRTexan = zeros(4,1);
    for k = 1:4
        nDTR = [0.1 0.01 0.001 0.0001];
        while true
            xDTRTexan(k,1) = xDTRTexan(k,1) + 1;
            PVratioDTRTexan = (1 - nDTR*xDTRTexan) * PVratioTexan;
            [PVcurDTRTexan,PVDTRTexan,selfDTRTexan,fedDTRTexan] = f_DTR(PVratioDTRTexan, weather);
            if max(PVcurDTRTexan) < 0.03
                break
            end
        end
        xDTRTexan(k,1) = xDTRTexan(k,1) - 1;
    end
    PVratioDTRTexan = (1 - nDTR*xDTRTexan - 0.0001) * PVratioTexan;
end

load CaseSwiss;
load Pload;
Pz = zeros(94,8760);
for j = 77:94
    Pz(j,:) = mpc.bus(j,3) * Pload';
end
Pload18 = sum(Pz(77:94,:),2);

Revenue_1year = zeros(1,2);
Revenue_1year(1,1) = sum(min([sum(PVSTRTexan,2) * 97/124.2*223.4 - sum(Pload18,2) * 223.4, zeros(18,1)],[],2) + sum(Pload18,2) * 223.4);
Revenue_1year(1,2) = sum(min([sum(PVDTRTexan,2) * 97/124.2*223.4 - sum(Pload18,2) * 223.4, zeros(18,1)],[],2) + sum(Pload18,2) * 223.4);

CostInv = [PVratioSTRTexan PVratioDTRTexan] * sum(mpc.bus(:,3)) * 2313000;
Subsidy = 0.26 * CostInv+ [PVratioSTRTexan PVratioDTRTexan] * sum(mpc.bus(:,3)) / 0.005 * 2500*2313000/2675000;;

year = 1:25;
Revenue = sum(1 ./ (1 + 0.05).^year) * Revenue_1year;
CostOM = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv;
NPVTexan = Revenue + Subsidy - CostInv - CostOM;


%% Chinese policy
CostInv0 = zeros(500,1);
selfcon0 = zeros(500,1);
PVR = zeros(500,1);
for k = 1:500
    PVR(k,1) = k/100;
    CostInv0(k,1) = PVR(k,1) * 2313000;
    selfcon0(k,1) = sum(min([PV_ava * PVR(k,1) - Pload, zeros(8760,1)]')' + Pload);
end
fed0 = sum(PV_ava) * PVR - selfcon0;
year = 1:25;
Revenue0_1year = 223.4 * selfcon0 + 37 * fed0;
Subsidy0_1year = 80/565.3*223.4 * sum(PV_ava) * PVR;

Revenue0 = sum(1 ./ (1 + 0.05).^year) * Revenue0_1year;
Subsidy0 = sum(1 ./ (1 + 0.05).^year) * Subsidy0_1year;
CostOM0 = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv0;
NPV0 = Revenue0 + Subsidy0 - CostInv0 - CostOM0;
PVratioChinese = find(NPV0==max(NPV0))/100;

% calculate the PV curtailment ratio with STR
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
        xDTRTexan(k,1) = xDTRTexan(k,1) - 1;
    end
    PVratioDTRChinese = (1 - nDTR*xDTRChinese - 0.0001) * PVratioChinese;
end

Revenue_1year = 223.4 * [sum(sum(selfSTRChinese)) sum(sum(selfDTRChinese))] ...
                + 37 * [sum(sum(fedSTRChinese)) sum(sum(fedDTRChinese))];
Subsidy_1year = 80/565.3*223.4 * [sum(sum(PVSTRChinese)) sum(sum(PVDTRChinese))];
CostInv = [PVratioSTRChinese PVratioDTRChinese] * sum(mpc.bus(:,3)) * 2313000;

year = 1:25;
Revenue = sum(1 ./ (1 + 0.05).^year) * Revenue_1year;
CostOM = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv;
Subsidy = sum(1 ./ (1 + 0.05).^year) * Subsidy_1year;
NPVChinese = Revenue + Subsidy - CostInv - CostOM;

