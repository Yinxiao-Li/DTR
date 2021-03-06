weather = xlsread('weather_2020.xlsx');
load Pload;

Ta_year = weather(1:8760,2);
W_year = weather(1:8760,1);
Tpv = Ta_year + 0.03 * W_year;
PV_ava = W_year / 1000 .* (1 - 0.0045 * (Tpv - 25));

PVR = (0.01:0.01:5); % MW
CostInv0 = PVR * 2675000;
Subsidy0 = 0.26 * CostInv0 + PVR / 0.005 * 2500;
Ploadyear = sum(Pload);
PVG = sum(PV_ava) * PVR;
Revenue0_1year = min([0.97 * PVG * 97 - Ploadyear * 124.2; ones(1,length(PVR))]) + Ploadyear * 124.2;
year = 1:25;
Revenue0 = sum(1 ./ (1 + 0.05).^year) * Revenue0_1year;
CostOM0 = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv0;
NPV0 = Revenue0 + Subsidy0 - CostInv0 - CostOM0;

PVratio = find(NPV0 == max(NPV0)) / 100;

% calculate the PV curtailment ratio with STR
[PVcurSTR,PVSTR] = f_STR(PVratio, weather);
PVratioSTR = PVratio;
if max(PVcurSTR) > 0.03
    xSTR = zeros(4,1);
    for k = 1:4
        nSTR = [0.1 0.01 0.001 0.0001];
        while true
            xSTR(k,1) = xSTR(k,1) + 1;
            PVratioSTR = (1 - nSTR*xSTR) * PVratio;
            [PVcurSTR,PVSTR] = f_STR(PVratioSTR, weather);
            if max(PVcurSTR) < 0.03
                break
            end
        end
        xSTR(k,1) = xSTR(k,1) - 1;
    end
    PVratioSTR = (1 - nSTR*xSTR - 0.0001) * PVratio;
end

% calculate the PV curtailment ratio with DTR
[PVcurDTR,PVDTR] = f_DTR(PVratio, weather);
PVratioDTR = PVratio;
if max(PVcurDTR) > 0.03
    xDTR = zeros(4,1);
    for k = 1:4
        nDTR = [0.1 0.01 0.001 0.0001];
        while true
            xDTR(k,1) = xDTR(k,1) + 1;
            PVratioDTR = (1 - nDTR*xDTR) * PVratio;
            [PVcurDTR,PVDTR] = f_DTR(PVratioDTR, weather);
            if max(PVcurDTR) < 0.03
                break
            end
        end
        xDTR(k,1) = xDTR(k,1) - 1;
    end
    PVratioDTR = (1 - nDTR*xDTR - 0.0001) * PVratio;
end

load CaseTexan;
load Pload;
Pz = zeros(53,8760);
for j = 30:48
    Pz(j,:) = mpc1.bus(j,3) * Pload';
end
Pload19 = sum(Pz(30:48,:),2);

Revenue_1year = zeros(1,2);
Revenue_1year(1,1) = sum(min([sum(PVSTR,2) * 97 - sum(Pload19,2) * 124.2, zeros(19,1)],[],2) + sum(Pload19,2) * 124.2);
Revenue_1year(1,2) = sum(min([sum(PVDTR,2) * 97 - sum(Pload19,2) * 124.2, zeros(19,1)],[],2) + sum(Pload19,2) * 124.2);

CostInv = [PVratioSTR PVratioDTR] * sum(mpc1.bus(:,3)) * 2675000;
Subsidy = 0.26 * CostInv+ [PVratioSTR PVratioDTR] * sum(mpc1.bus(:,3)) / 0.005 * 2500;

year = 1:25;
Revenue = sum(1 ./ (1 + 0.05).^year) * Revenue_1year;
CostOM = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv;
NPV = Revenue + Subsidy - CostInv - CostOM;

