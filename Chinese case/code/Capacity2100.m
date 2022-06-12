weather2025 = xlsread('weather_2025_1.xlsx');
Tamon_max2025 = max([mean(weather2025(181*24+1:212*24,2)) mean(weather2025(212*24+1:243*24,2))]);
Ta_max2025 = max(weather2025(:,3));
    
weather2100 = xlsread('weather_2100_1.xlsx');
Tamon_max2100 = max([mean(weather2100(181*24+1:212*24,2)) mean(weather2100(212*24+1:243*24,2))]);
Ta_max2100 = max(weather2100(:,3));

ktran = 1 - (Tamon_max2100 - Tamon_max2025) / 100;
D1 = 13.6 / 1000;
r1 = 0.4045;
u = 0.61;
v = 0;
Ta = 40;
E = 0.8;
a = 0.8;
Tc = 80;
W = 1000;
kline = DLR(D1, u, v, Ta+Ta_max2100-Ta_max2025,E, r1,  W, a, Tc) / DLR(D1, u, v, Ta,E, r1,  W, a, Tc);


% ktran is the ratio of the transformer nameplate rating in 2100 to the
% transformer nameplate rating in 2025
% kline is the ratio of the overhead line nameplate rating in 2100 to the
% overhead line nameplate rating in 2025