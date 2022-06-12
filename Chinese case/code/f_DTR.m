function [PVcurratio,PVPG,selfcon,fed] = f_DTR(PVratio, weather)

load CaseChinese;
trans = (49:96)';
lines = (1:48)';
D1 = 13.6 / 1000; 
r1 = 0.4045; 
Sr = mpc.branch(trans,6); 

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

Pz = zeros(97,8760);
Qz = zeros(97,8760);

for i = 1:365
    if any(springrain == i)
        for j = 1:97
            Pz(j,24*(i-1)+1:24*i) = mpc.bus(j,3) * demand(:,1)';
            Qz(j,24*(i-1)+1:24*i) = mpc.bus(j,4) * demand(:,1)';
        end
    end
    if any(springsunny == i)
        for j = 1:97
            Pz(j,24*(i-1)+1:24*i) = mpc.bus(j,3) * demand(:,2)';
            Qz(j,24*(i-1)+1:24*i) = mpc.bus(j,4) * demand(:,2)';
        end
    end
    if any(summerrain == i)
        for j = 1:97
            Pz(j,24*(i-1)+1:24*i) = mpc.bus(j,3) * demand(:,3)';
            Qz(j,24*(i-1)+1:24*i) = mpc.bus(j,4) * demand(:,3)';
        end
    end
    if any(summersunny == i)
        for j = 1:97
            Pz(j,24*(i-1)+1:24*i) = mpc.bus(j,3) * demand(:,4)';
            Qz(j,24*(i-1)+1:24*i) = mpc.bus(j,4) * demand(:,4)';
        end
    end
    if any(fallsunny == i)
        for j = 1:97
            Pz(j,24*(i-1)+1:24*i) = mpc.bus(j,3) * demand(:,5)';
            Qz(j,24*(i-1)+1:24*i) = mpc.bus(j,4) * demand(:,5)';
        end
    end
    if any(fallrain == i)
        for j = 1:97
            Pz(j,24*(i-1)+1:24*i) = mpc.bus(j,3) * demand(:,6)';
            Qz(j,24*(i-1)+1:24*i) = mpc.bus(j,4) * demand(:,6)';
        end
    end
    if any(winterrain == i)
        for j = 1:97
            Pz(j,24*(i-1)+1:24*i) = mpc.bus(j,3) * demand(:,7)';
            Qz(j,24*(i-1)+1:24*i) = mpc.bus(j,4) * demand(:,7)';
        end
    end
    if any(wintersunny == i)
        for j = 1:97
            Pz(j,24*(i-1)+1:24*i) = mpc.bus(j,3) * demand(:,8)';
            Qz(j,24*(i-1)+1:24*i) = mpc.bus(j,4) * demand(:,8)';
        end
    end
end

PVC = PVratio * mpc.bus(50:97,3);

mpc.gen(2:length(PVC)+1,1:21) = zeros(length(PVC),21);
mpc.gen(2:length(PVC)+1,1) = (50:97)';
mpc.gen(:,7) = 100;
mpc.gen(:,8) = 1;
mpc.gencost(2:length(PVC)+1,1:7) = zeros(length(PVC),7);
mpc.gencost(2:length(PVC)+1,1) = 2;
mpc.gencost(2:length(PVC)+1,4) = 3;
mpc.gencost(2:length(PVC)+1,6) = -0.1;

PVPG = zeros(length(PVC),8760);
PV1 = zeros(8760,1);

selfcon = zeros(length(PVC),8760);
fed = zeros(length(PVC),8760);

% transformer parameters
delta_t = 60;
R = 8;
delta_Tto_R = 55;
delta_Th_R = 25;
Tao_to = 155; 
Tao_h = 5;
Cto = (1 - exp(-delta_t / Tao_to));
Ch = (1 - exp(-delta_t / Tao_h));

u_year = weather(:,1);
v_year = weather(:,2);
Ta_year = weather(:,3);
W_year = weather(:,4);

Ths = zeros(length(trans),8761);
delta_Tto = zeros(length(trans),8761);
delta_Ths = zeros(length(trans),8761);
for i = 1:length(trans)
    a = find(mpc.bus(:,1)==mpc.branch(trans(i),2));
    Ki = (Pz(a,8760)^2 + Qz(a,8760)^2)^0.5 / Sr(i);
    delta_Tto(i,1) = delta_Tto_R * (Ki^2 * R + 1) / (R + 1);
    delta_Ths(i,1) = delta_Th_R * Ki^2;
    Ths(i,1) = Ta_year(8760,1) + delta_Tto(i,1) + delta_Ths(i,1);
end

for j = 1:8760
    % calculate maximum available PV output
    Tpv = Ta_year(j) + 0.03 * W_year(j);
    PV_ava = W_year(j) / 1000 .* (1 - 0.0045 * (Tpv - 25));
    PV1(j,1) = PV_ava;
        
    mpc.bus(:,3) = Pz(:,j);
    mpc.bus(:,4) = Qz(:,j);
    for i = 1:length(PVC)
        mpc.gen(i+1,2) = PV_ava * PVC(i);
        PVPG(i,j) = PV_ava * PVC(i);
    end
    mpopt = mpoption('out.all', 0 ,'verbose', 0);
    [result,~] = runpf(mpc,mpopt);
    
    % calculate transformer temperature
    for i = 1:length(trans)
        Ki(i,j) = (result.branch(trans(i),16)^2 + result.branch(trans(i),17)^2)^0.5 / Sr(i);
        Ki2 = Ki(i,j)^2;
        delta_Tto(i,j+1) = (delta_Tto_R * (Ki2 * R + 1) / (R + 1) - delta_Tto(i,j)) * Cto + delta_Tto(i,j);
        delta_Ths(i,j+1) = (delta_Th_R * Ki2 - delta_Ths(i,j)) * Ch + delta_Ths(i,j);
        Ths(i,j+1) = Ta_year(j,1) + delta_Tto(i,j+1) + delta_Ths(i,j+1);
    end
    
    % calculate overhead line capacity
    DLR_I = DLR(D1, u_year(j), v_year(j),   Ta_year(j),0.8,  r1,    W_year(j), 0.8,  80);
    DLR_S = DLR_I * 3^0.5 * 10 / 1000; 
    
    deltaI = result.branch(lines,16).^2 + result.branch(lines,17).^2 - DLR_S.^2;
    % If any constraint is violated, the optimal power flow is executed
    if max(result.bus(:,8)) > 1.1 || max(Ths(:,j+1)) > 120 || max(deltaI) > 0
        mpc.branch(lines,6) = DLR_S;
        xPVcur = zeros(3,1);
        for k = 1:3
            nPVcur = [0.1 0.01 0.001];
            while true
                xPVcur(k,1) = xPVcur(k,1) + 1;
                mpc.gen(2:length(PVC)+1,2) = (1 - nPVcur*xPVcur) * PV_ava * PVC;
                [result,~] = runpf(mpc,mpopt);
                for i = 1:length(trans)
                    Ki(i,j) = (result.branch(trans(i),16)^2 + result.branch(trans(i),17)^2)^0.5 / Sr(i);
                    Ki2 = Ki(i,j)^2;
                    delta_Tto(i,j+1) = (delta_Tto_R * (Ki2 * R + 1) / (R + 1) - delta_Tto(i,j)) * Cto + delta_Tto(i,j);
                    delta_Ths(i,j+1) = (delta_Th_R * Ki2 - delta_Ths(i,j)) * Ch + delta_Ths(i,j);
                    Ths(i,j+1) = Ta_year(j,1) + delta_Tto(i,j+1) + delta_Ths(i,j+1);
                end
                deltaI = result.branch(lines,16).^2 + result.branch(lines,17).^2 - result.branch(lines,6).^2;
                
                if max(result.bus(:,8)) < 1.1 && max(deltaI) < 0 && max(Ths(:,j+1)) < 120
                    break
                end
            end
            xPVcur(k,1) = xPVcur(k,1) - 1;
        end
        
        PVPG(:,j) = result.gen(2:length(PVC)+1,2);
    end
    
    for i = 1:length(trans)
        if Pz(find(mpc.bus(:,1)==mpc.branch(trans(i),2)),j) > PVPG(i,j)
            selfcon(i,j) = PVPG(i,j);
        else
            selfcon(i,j) = Pz(find(mpc.bus(:,1)==mpc.branch(trans(i),2)),j);
            fed(i,j) = PVPG(i,j) - Pz(find(mpc.bus(:,1)==mpc.branch(trans(i),2)),j);
        end
    end
end

PV1sum = sum(PV1) * PVC;
PVcurratio = 1 - sum(PVPG,2) ./ PV1sum;

