% ��758���µ��¶�
function [Ht] = cabletem758(H0, Ta, ps, ds, I)

Sc = 630;  % ��������
dc = 28.3;  % ����ֱ��
R0 = 0.0428 * 10^(-3);  % ֱ������(��/m)
Ti = 2.4;
l1 = 0.042689413008005;
U0 = 420/3^0.5;  % �Եص�ѹ�����ѹ��
f = 50;  % Ƶ��

%% ���µ��⾶
di = dc + 2 * Ti;  % ��Եֱ��
dp = di + 0.2 * 2;  % ͭ������ֱ��
% dp = dp + 2 * 0.2; %%%%%%%%%%%% �������
Tm = 0.035 * dp + 1;  % �⻤����
De = dp + 2 * Tm;  % �����⾶

%% �������
Rm = R0 * (1 + 0.00403 * (90 - 20));  % ����¶�(90)����
xs2 = 8 * pi * f * 10^(-7) / Rm;
ys = xs2^2 / (192 + 0.8 * xs2^2);
xp2 = 8 * pi * f * 10^(-7) / Rm;  % ����ЧӦ����
s = De;  % �������ļ��
%  �ٽ�ЧӦ����
yp = xp2^2 / (192 + 0.8 * xp2^2) * (dc / s)^2 * (0.312 * (dc / s)^2 + 1.18 / (xp2^2 / (192 + 0.8 * xp2^2) + 0.27));

%% �������
T1 = 3.5 / 2 / pi * log(di / dc);
T3 = 6 / 2 / pi * log(De / dp);
re = De / 2;
T4 = 1.5 * ps / pi * (log(2 * 800 / re) - 0.63);

%% ����
Qc = 2.5 * Sc;
Qi = 2.4 * pi * (di^2 - dc^2) / 4;
Qs = 3.45 * pi * (dp^2 - di^2) / 4;
Qj = 1.7 * pi * (De^2 - dp^2) / 4;

%% A
re = De / 2;
rc = dc / 2;
pa = 1 / 2 / log(re / rc) - 1 / (re^2 / rc^2 - 1);
QA = Qc + pa * (Qi + (Qs + Qj) / (1 + l1));
TA = T1 + T3 * (1 + l1);

%% B
ktf = 3 * (log(2 * 800 / re) - 0.63) / log(800 / re + ((800 / re)^2 - 1)^0.5);
rb = 3 * re;
pb = 1 / 2 / log(rb / re) - 1 / (rb^2 / re^2 - 1);
TB = ktf * ps * (1 + l1) / 2 / pi * log(rb / re);
QB = (1 - pa) * (Qi + (Qs + Qj) / (1 + l1)) + pb * pi * pi * (rb^2 - re^2) * 10^(-6) / ktf / ps / ds / (1 + l1);

%% C
rc = 3 * rb;
pc = 1 / 2 / log(rc / rb) - 1 / (rc^2 / rb^2 - 1);
TC = ktf * ps * (1 + l1) / 2 / pi * log(rc / rb);
QC = (1 - pb) * pi * (rb^2 - re^2) * 10^(-6) / ktf / ps / ds / (1 + l1) + pc * pi * (rc^2 - rb^2) * 10^(-6) / ktf / ps / ds / (1 + l1);

%% D
renv = re * exp(2 * pi * T4 / ktf / ps);
penv = 1 / 2 / log(renv / rc) - 1 / (renv^2 / rc^2 - 1);
TD = T4 * (1 + l1) - (TB + TC);
QD = (1 - pc) * pi * (rc^2 - rb^2) * 10^(-6) / ktf / ps / ds / (1 + l1) + penv * pi * (renv^2 - rc^2) * 10^(-6) / ktf / ps / ds / (1 + l1);

%% ��Ե��ĵ��µ�����
c = 2.5 / 18 / log(di / dc) * 10^(-9);  % ��λ���ȵ��µ���
Wd = 2 * pi * f * c * U0^2 * 0.004;  % �������
hd = Wd * (0.5 * T1 + T3 + T4);

%% ���¶�
t = 60;
min = 60;

% ϵ������
xsjz = [-1/QD/TD-1/QD/TC, 1/QD/TC, 0, 0
        1/QC/TC, -1/QC/TC-1/QC/TB, 1/QC/TB, 0
        0, 1/QB/TB, -1/QB/TB-1/QB/TA, 1/QB/TA
        0, 0, 1/QA/TA, -1/QA/TA];
[v, l] = eig(xsjz);  % l����ֵ v��������

H1 = H0;

for i = 1:min
    R = R0 * (1 + ys + yp) * (1 + 0.00403 * (H1(4,1) + Ta + hd - 20));
    Wc = I^2 * R;
    
    % �¶ȵ���ֵ
    hDz = Wc * TD;
    hBz = hDz + Wc * TC;
    hsz = hBz + Wc * TB;
    hcz = hsz + Wc * TA;
    HZ = [hDz; hBz; hsz; hcz];
    c = v \ (H1 - HZ);
    Ht = c(1) * v(:,1) * exp(l(1,1) * t) + c(2) * v(:,2) * exp(l(2,2) * t) + c(3) * v(:,3) * exp(l(3,3) * t) + c(4) * v(:,4) * exp(l(4,4) * t) + HZ;
    H1 = Ht;
end





