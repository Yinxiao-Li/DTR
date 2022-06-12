function [PVcurratio,PVPG,selfcon,fed] = f_STR(PVratio, weather)

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

Ta_year = weather(1:8760,2);
W_year = weather(1:8760,1);

Pz = zeros(94,8760);
Qz = zeros(94,8760);
for i = 77:94
    Pz(i,:) = mpc.bus(i,3) * Pload';
    Qz(i,:) = mpc.bus(i,4) * Pload';
end

PVC = PVratio * mpc.bus(77:94,3);

PVPG = zeros(18,8760);
PV1 = zeros(8760,1);

selfcon = zeros(18,8760);
fed = zeros(18,8760);

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
    
    % If any constraint is violated, the optimal power flow is executed
    deltaI = result.branch(:,16).^2 + result.branch(:,17).^2 - result.branch(:,6).^2;
    if max(result.bus(:,8)) > 1.05 || max(deltaI) > 0
        xPVcur = zeros(3,1);
        for k = 1:3
            nPVcur = [0.1 0.01 0.001];
            while true
                xPVcur(k,1) = xPVcur(k,1) + 1;
                mpc.gen(2:19,2) = (1 - nPVcur*xPVcur) * PV_ava * PVC;
                [result,~] = runpf(mpc,mpopt);
                deltaI = result.branch(:,16).^2 + result.branch(:,17).^2 - result.branch(:,6).^2;
                if max(result.bus(:,8)) < 1.05 && max(deltaI) < 0
                    break
                end
            end
            xPVcur(k,1) = xPVcur(k,1) - 1;
        end
        PVPG(:,j) = result.gen(2:19,2);
    end
    
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
