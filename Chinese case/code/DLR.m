% calculate overhead line capacity
function [Imax]=DLR(D, u, v,   Ta,E,  r,    W,      a,  Tc)

Tf=(Tc + Ta) / 2;
uf=1.458*10^(-6)*(Tf+273)^1.5/(Tf+383.4);
pf=1.029;
kf=2.424*0.01+7.477*10^-5*Tf-4.407*10^-9*Tf^2;
V = (u^2 + v^2)^0.5;
q = abs(atan(u/v));
NRe=D*pf*V/uf;
Kangle=1.194-cos(q)+0.194*cos(2*q)+0.368*sin(2*q);
qs=a*W*D;
R=r/1000;
RTc=R * (1 + 0.00403 * (Tc - 20));

qc1 = Kangle * (1.01+1.35*NRe^0.52) * kf * (Tc - Ta);
qc2 = Kangle * 0.754 * NRe^0.6 * kf * (Tc - Ta);
qc = max([qc1 qc2]);

qr = 17.8 * D * E * ( ((Tc + 273)/100)^4 - ((Ta + 273)/100)^4 );

Imax=((qc+qr-qs)/RTc)^0.5;

end