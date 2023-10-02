clc
HCIS = [0,1,0,0,1,0,1,0];
addpath('C:\Users\XYZW\Documents\Matlab dissertation - new version\Functions')
addpath('C:\Users\XYZW\Documents\Matlab dissertation - new version\Markowitz_problems\risk_contributions_budget')
CE1 = [75,5000,720,50,2500,25,30000,5];
CE2 = [75,5000,1030,350,4500,5,2000,64];
CE3 = [24000,15000,1210,550,500,187,30000,199];
Y = [300,320,125,100,200,102,100,25];
b = [0.22,0.19,0.17,0.13,0.11,0.08,0.06,0.04];
W = 1;MktCap_total = 2500;
carb_intens1 = CE1./Y;carb_intens2 = CE2./Y;carb_intens3 = CE3./Y;
carb_intens12 = (CE1+CE2)./Y;carb_intens123 = (CE1+CE2+CE3)./Y;
mkt_caps = b*MktCap_total;
%%
clc
sigs = [0.22,0.2,0.25,0.18,0.35,0.23,0.13,0.29];
carb_intens_12 = carb_intens1+carb_intens2;
carb_intens_123 = carb_intens1+carb_intens2+carb_intens3;
ro = [1,0.8,0.7,0.6,0.7,0.5,0.7,0.6;0.8,1,0.75,0.65,0.5,0.6,0.5,0.65;...
    0.7,0.75,1,0.8,0.7,0.7,0.7,0.7;0.6,0.65,0.8,1,0.85,0.8,0.75,0.75;...
    0.7,0.5,0.7,0.85,1,0.6,0.8,0.65;0.5,0.6,0.7,0.8,0.6,1,0.5,0.7;0.7,0.5,0.7,0.75,0.8,0.5,1,0.8;...
    0.6,0.65,0.7,0.75,0.65,0.7,0.8,1];
%%
clc
str_func_CI_decarb_R1_intensity = min_risk_decarb(sigs,ro,b,CE1,Y,0.2,[],'type','relative','carbon','WACI');
str_func_CI_decarb_R12_intensity = min_risk_decarb(sigs,ro,b,CE1+CE2,Y,0.2,[]);
str_func_CI_decarb_R123_intensity = min_risk_decarb(sigs,ro,b,CE1+CE2+CE3,Y,0.2,[]);
%%
str_func_CI_decarb_R1_emission = min_risk_decarb(sigs,ro,b,CE1,Y,0.2,[],'type','relative','carbon','emission');
str_func_CI_decarb_R12_emission = min_risk_decarb(sigs,ro,b,CE1+CE2,Y,0.2,[],'type','relative','carbon','emission');
str_func_CI_decarb_R123_emission = min_risk_decarb(sigs,ro,b,CE1+CE2+CE3,Y,0.2,[],'type','relative','carbon','emission');
%%
clc
[str_CI_decarb_S12_order_stat,ret_CI_decarb_S12_order_stat,vol_CI_S12,RC_CI_S12,~,red12] = ...
    min_risk_decarb_order_stat(sigs,ro,b,CE1+CE2,Y,2);
[str_CI_decarb_S123_order_stat,ret_CI_decarb_S123_order_stat,vol_CI_S123,RC_CI_S123,~,red123] = ...
    min_risk_decarb_order_stat(sigs,ro,b,CE1+CE2+CE3,Y,2);
%%
str_func_CI_decarb_R1_intensity_HCIS = ...
    min_risk_decarb(sigs,ro,b,CE1,Y,0.2,HCIS,'type','relative','carbon','WACI');
str_func_CI_decarb_R12_intensity_HCIS = ...
    min_risk_decarb(sigs,ro,b,CE1+CE2,Y,0.2,HCIS,'type','relative','carbon','WACI');
str_func_CI_decarb_R123_intensity_HCIS = ...
    min_risk_decarb(sigs,ro,b,CE1+CE2+CE3,Y,0.2,HCIS,'type','relative','carbon','WACI');
%%
str_func_CI_decarb_R1_emission_HCIS = ...
    min_risk_decarb(sigs,ro,b,CE1,Y,0.2,HCIS,'type','relative','carbon','emission');
str_func_CI_decarb_R12_emission_HCIS = ...
    min_risk_decarb(sigs,ro,b,CE1+CE2,Y,0.2,HCIS,'type','relative','carbon','emission');
str_func_CI_decarb_R123_emission_HCIS = ...
    min_risk_decarb(sigs,ro,b,CE1+CE2+CE3,Y,0.2,HCIS,'type','relative','carbon','emission');