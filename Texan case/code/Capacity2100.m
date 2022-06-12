weather2025 = xlsread('weather_2025_1.xlsx');
Tamon_max2025 = max([mean(weather2025(181*24+1:212*24,2)) mean(weather2025(212*24+1:243*24,2))]);
Ts_max2025 = max(weather2025(:,4));
    
weather2100 = xlsread('weather_2100_1.xlsx');
Tamon_max2100 = max([mean(weather2100(181*24+1:212*24,2)) mean(weather2100(212*24+1:243*24,2))]);
Ts_max2100 = max(weather2100(:,4));

ktran = 1 - (Tamon_max2100 - Tamon_max2025) / 100;
kcable = (90 - Ts_max2100)^0.5 / (90 - Ts_max2025)^0.5;

% ktran is the ratio of the transformer nameplate rating in 2100 to the
% transformer nameplate rating in 2025
% kcable is the ratio of the cable nameplate rating in 2100 to the
% cable nameplate rating in 2025