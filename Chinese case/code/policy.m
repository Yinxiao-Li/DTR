weather = xlsread('weather_2020.xlsx');
load CaseChinese;
Ta_year = weather(1:8760,3);
W_year = weather(1:8760,4);
Tpv = Ta_year + 0.03 * W_year;
PV_ava = W_year / 1000 .* (1 - 0.0045 * (Tpv - 25));

pr = 0.01;
Pre_year = weather(:,5);
daypre = zeros(365,1);
springrain = [];
springsunny = [];
summerrain = [];
summersunny = [];
fallrain = [];
fallsunny = [];
winterrain = [];
wintersunny = [];

for i = 1:365
    daypre(i,1) = sum(Pre_year(24*i-23:24*i,1));
    if i >= 60 && i <= 151
        if daypre(i,1) > pr
            springrain(length(springrain)+1,1) = i;
        else
            springsunny(length(springsunny)+1,1) = i;
        end
    elseif i >= 152 && i <= 243
        if daypre(i,1) > pr
            summerrain(length(summerrain)+1,1) = i;
        else
            summersunny(length(summersunny)+1,1) = i;
        end
    elseif i >= 244 && i <= 334
        if daypre(i,1) > pr
            fallrain(length(fallrain)+1,1) = i;
        else
            fallsunny(length(fallsunny)+1,1) = i;
        end
    else
        if daypre(i,1) > pr
            winterrain(length(winterrain)+1,1) = i;
        else
            wintersunny(length(wintersunny)+1,1) = i;
        end
    end
end

demand = [0.194595	0.274202	0.471744	0.503686	0.280098	0.238821	0.347912	0.33317
0.247665	0.274202	0.468796	0.449631	0.277149	0.247665	0.321375	0.300737
0.235872	0.262407	0.433415	0.437346	0.265356	0.232923	0.327273	0.283047
0.241769	0.268305	0.400983	0.395577	0.262407	0.227027	0.324324	0.288944
0.244718	0.265356	0.389189	0.375921	0.271253	0.268305	0.309582	0.283047
0.297789	0.448157	0.445209	0.304668	0.277149	0.300737	0.356757	0.29484
0.448157	0.448157	0.513022	0.289926	0.421622	0.474693	0.33317	0.324324
0.400983	0.306633	0.436364	0.201474	0.448157	0.380343	0.44226	0.380343
0.353808	0.150369	0.342015	0.277641	0.415725	0.285995	0.610319	0.471744
0.297789	0.011793	0.412776	0.287469	0.359705	0.336117	0.610319	0.353808
0.462899	0.038329	0.49828	0.299754	0.468796	0.315479	0.227027	0.232923
0.023587	0.044226	0.465848	0.314496	0.377396	0.085503	0.557248	0.244718
0.206388	0.082549	0.44226	0.346437	0.330221	0.194595	0.507125	0.061917
0.200492	0.026533	0.324324	0.432432	0.315479	0.262407	0.457002	0.109091
0.141524	0.008845	0.33317	0.479115	0.33317	0.206388	0.468796	0.347912
0.250614	0.126782	0.342015	0.570025	0.377396	0.306633	0.504177	0.359705
0.474693	0.356757	0.386241	0.599509	0.580835	0.483538	0.586732	0.436364
0.533661	0.557248	0.49828	0.685504	0.636856	0.507125	0.66634	0.592629
0.619165	0.751843	0.504177	0.685504	0.569042	0.5543	0.816707	0.70172
0.577887	0.772481	0.734152	1	0.489435	0.48059	0.908108	0.79312
0.454054	0.577887	0.613268	0.977887	0.40688	0.400983	0.790172	0.692875
0.342015	0.457002	0.501229	0.781327	0.315479	0.324324	0.689926	0.545455
0.274202	0.386241	0.418673	0.675676	0.309582	0.280098	0.513022	0.398034
0.256511	0.315479	0.356757	0.594595	0.271253	0.244718	0.398034	0.33317
];
Pload = zeros(8760,1);

for i = 1:365
    if any(springrain == i)
        Pload(24*(i-1)+1:24*i,1) = demand(:,1);
    end
    if any(springsunny == i)
        Pload(24*(i-1)+1:24*i,1) = demand(:,2);
    end
    if any(summerrain == i)
        Pload(24*(i-1)+1:24*i,1) = demand(:,3);
    end
    if any(summersunny == i)
        Pload(24*(i-1)+1:24*i,1) = demand(:,4);
    end
    if any(fallsunny == i)
        Pload(24*(i-1)+1:24*i,1) = demand(:,5);
    end
    if any(fallrain == i)
        Pload(24*(i-1)+1:24*i,1) = demand(:,6);
    end
    if any(winterrain == i)
        Pload(24*(i-1)+1:24*i,1) = demand(:,7);
    end
    if any(wintersunny == i)
        Pload(24*(i-1)+1:24*i,1) = demand(:,8);
    end
end

% Texan policy

PVR = (0.01:0.01:5); % MW
CostInv0 = PVR * 5250000;
Subsidy0 = 0.26 * CostInv0 + PVR / 0.005 * 2500*5250000/2675000;
Ploadyear = sum(Pload);
PVG = sum(PV_ava) * PVR;
Revenue0_1year = min([PVG * 97/124.2*565.3 - Ploadyear * 565.3; ones(1,length(PVR))]) + Ploadyear * 565.3;
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

load CaseChinese;
loadbus = (50:97)';

Pz = zeros(97,8760);
for j = loadbus
    Pz(j,:) = mpc.bus(j,3) * Pload';
end
Ploadsum = sum(Pz(loadbus,:),2);

Revenue_1year = zeros(1,2);
Revenue_1year(1,1) = sum(min([sum(PVSTRTexan,2) * 97/124.2*565.3 - sum(Ploadsum,2) * 565.3, zeros(48,1)],[],2) + sum(Ploadsum,2) * 565.3);
Revenue_1year(1,2) = sum(min([sum(PVDTRTexan,2) * 97/124.2*565.3 - sum(Ploadsum,2) * 565.3, zeros(48,1)],[],2) + sum(Ploadsum,2) * 565.3);

CostInv = [PVratioSTRTexan PVratioDTRTexan] * sum(mpc.bus(:,3)) * 5250000;
Subsidy = 0.26 * CostInv+ [PVratioSTRTexan PVratioDTRTexan] * sum(mpc.bus(:,3)) / 0.005 * 2500 * 5250000 / 2675000;

year = 1:25;
Revenue = sum(1 ./ (1 + 0.05).^year) * Revenue_1year;
CostOM = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv;
NPVTexan = Revenue + Subsidy - CostInv - CostOM;


%% Swiss policy
CostInv0 = zeros(500,1);
Subsidy0 = zeros(500,1);
selfcon0 = zeros(500,1);
PVR = zeros(500,1);
for k = 1:500
    PVR(k,1) = k/100;
    CostInv0(k,1) = PVR(k,1) * 5250000;
    Subsidy0(k,1) = 0.3 * CostInv0(k,1);
    selfcon0(k,1) = sum(min([PV_ava * PVR(k,1) - Pload, zeros(8760,1)]')' + Pload);
end
fed0 = sum(PV_ava) * PVR - selfcon0;
year = 1:25;
Revenue0_1year = 565.3 * selfcon0 + 384.4 * fed0;
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

Revenue_1year = 565.3 * [sum(sum(selfSTRSwiss)) sum(sum(selfDTRSwiss))] + ...
                384.4 * [sum(sum(fedSTRSwiss)) sum(sum(fedDTRSwiss))];
CostInv = [PVratioSTRSwiss PVratioDTRSwiss] * sum(mpc.bus(:,3)) * 5250000;
Subsidy = 0.3 * CostInv;
year = 1:25;
Revenue = sum(1 ./ (1 + 0.05).^year) * Revenue_1year;
CostOM = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv;
NPVSwiss = Revenue + Subsidy - CostInv - CostOM;

