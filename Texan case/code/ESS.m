weather = xlsread('weather_2020.xlsx');
load CaseTexan;
load Pload;
PVratio = 3.473288;

Ta_year = weather(1:8760,2);
W_year = weather(1:8760,1);

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

selfcon = zeros(length(PVC),8760);
fed = zeros(length(PVC),8760);
PVPG = zeros(19,8760);
PV1 = zeros(8760,1);

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
    
    % If any constraint is violated, calculate the curtailed PV generation
    deltaI = result.branch(:,16).^2 + result.branch(:,17).^2 - result.branch(:,6).^2;
    if max(result.bus(:,8)) > 1.05 || max(deltaI) > 0

        xPVcur = zeros(3,1);
        for k = 1:3
            nPVcur = [0.1 0.01 0.001];
            while true
                xPVcur(k,1) = xPVcur(k,1) + 1;
                mpc.gen(2:20,2) = (1 - nPVcur*xPVcur) * PV_ava * PVC;
                [result,~] = runpf(mpc,mpopt);
                deltaI = result.branch(:,16).^2 + result.branch(:,17).^2 - result.branch(:,6).^2;
                if max(result.bus(:,8)) < 1.05 && max(deltaI) < 0
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

% calculate the capacity of ESSs
load CaseTexan;

xESS = zeros(3,1);
for k = 1:3
    nESS = [0.1 0.01 0.001];
    while true
        xESS(k,1) = xESS(k,1) + 1;
        ESSratio = nESS*xESS;
        ESSC = ESSratio * mpc1.bus(30:48,3);
        
        Pmax = 0.5 * ESSC;
        nita = 0.86^0.5;
        SOC = zeros(19,8761);
        SOC(:,1) = 0.5 * ESSC;
        Pch = zeros(19,8760);
        Pdis = zeros(19,8760);
        
        for j = 1:8760
            PVexa = PV1(j) * PVC - PVPG(:,j);
            for i = 1:19
                Pch(i,j) = PVexa(i);
                if PVexa(i) > Pmax(i) || Pch(i,j) * nita + SOC(i,j) > ESSC(i)
                    Pch(i,j) = min([Pmax(i); (ESSC(i) - SOC(i,j)) / nita]);
                end
            end

            for i = 1:19
                if PVPG(i,j) < Pz(i+29,j) && SOC(i,j) > 0.1 * ESSC(i)
                    Pdis(i,j) = min([Pmax(i);(SOC(i,j) - 0.1 * ESSC(i)) * nita]);
                end
            end
            SOC(:,j+1) = SOC(:,j) + Pch(:,j) * nita - Pdis(:,j) / nita;
        end
        
        PVcurratio = 1 - sum(PVPG+Pch,2) ./ PV1sum;
        
        if max(PVcurratio) < 0.03
        	break
        end
    end
    xESS(k,1) = xESS(k,1) - 1;
end
ESSratio = nESS*xESS + 0.001;
ESSC = ESSratio * mpc1.bus(30:48,3);

