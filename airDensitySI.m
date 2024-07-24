function [rhow] = airDensity(Td,T,P)

% convert to SI units
Td = (Td-32)*5/9;
T  = (T-32) *5/9 + 273.15;

% determine vapour pressure 
c0 =   0.99999683;
c1 =  -0.90826951*10^-2;
c2 =   0.78736169*10^-4;
c3 =  -0.61117958*10^-6;
c4 =   0.43884187*10^-8;
c5 =  -0.29883885*10^-10;
c6 =   0.21874425*10^-12;
c7 =  -0.17892321*10^-14;
c8 =   0.11112018*10^-16;
c9 =  -0.30994571*10^-19;

es0 = 6.1078;
Pv = c0 + Td*c1 + Td^2*c2 + Td^3*c3 + Td^4*c4 + Td^5*c5 + Td^6*c6 + Td^7*c7 + Td^8*c8 + Td^9*c9;
Pv = es0/(Pv^8);

% find dry air pressure from total pressure
Pd = P - Pv;

% calculate density
rhow = 100*Pd/(287.0531*T) + 100*Pv/(461.4964*T);
% convert from kg/m^3 to slugs/ft^3 
%rhow = rhow /14.5939 * (2.54*12/100)^3;