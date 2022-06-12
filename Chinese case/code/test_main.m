weather = xlsread('weather_2020.xlsx');
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

CostInv0 = zeros(500,1);
selfcon0 = zeros(500,1);
PVR = zeros(500,1);
for k = 1:500
    PVR(k,1) = k/100;
    CostInv0(k,1) = PVR(k,1) * 5250000;
    selfcon0(k,1) = sum(min([PV_ava * PVR(k,1) - Pload, zeros(8760,1)]')' + Pload);
end
fed0 = sum(PV_ava) * PVR - selfcon0;
year = 1:25;
Revenue0_1year = 565.3 * selfcon0 + 384.4 * fed0;
Subsidy0_1year = 80 * sum(PV_ava) * PVR;

Revenue0 = sum(1 ./ (1 + 0.05).^year) * Revenue0_1year;
Subsidy0 = sum(1 ./ (1 + 0.05).^year) * Subsidy0_1year;
CostOM0 = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv0;
NPV0 = Revenue0 + Subsidy0 - CostInv0 - CostOM0;
PVratio = find(NPV0==max(NPV0))/100;

% calculate the PV curtailment ratio with STR
[PVcurSTR,PVSTR,selfSTR,fedSTR] = f_STR(PVratio, weather);

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

[PVcurDTR,PVDTR,selfDTR,fedDTR] = f_DTR(PVratio, weather);

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

load CaseChinese;
Revenue_1year = 565.3 * [sum(sum(selfSTR)) sum(sum(selfDTR))] ...
                + 384.4 * [sum(sum(fedSTR)) sum(sum(fedDTR))];
Subsidy_1year = 80 * [sum(sum(PVSTR)) sum(sum(PVDTR))];
CostInv = [PVratioSTR PVratioDTR] * sum(mpc.bus(:,3)) * 5250000;

year = 1:25;
Revenue = sum(1 ./ (1 + 0.05).^year) * Revenue_1year;
CostOM = sum(1 ./ (1 + 0.05).^year) * 0.01 * CostInv;
Subsidy = sum(1 ./ (1 + 0.05).^year) * Subsidy_1year;
NPV = Revenue + Subsidy - CostInv - CostOM;
