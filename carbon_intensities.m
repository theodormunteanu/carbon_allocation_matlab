CE1 = [75,5000,720,50,2500,25,30000,5];
CE2 = [75,5000,1030,350,4500,5,2000,64];
CE3 = [24000,15000,1210,550,500,187,30000,199];
Y = [300,320,125,100,200,102,100,25];
sigs = [0.22,0.2,0.25,0.18,0.35,0.23,0.13,0.29];
carb_intens1 = CE1./Y;carb_intens2 = CE2./Y;carb_intens3 = CE3./Y;
carb_intens_12 = carb_intens1+carb_intens2;
carb_intens_123 = carb_intens1+carb_intens2+carb_intens3;
ro = [1,0.8,0.7,0.6,0.7,0.5,0.7,0.6;0.8,1,0.75,0.65,0.5,0.6,0.5,0.65;...
    0.7,0.75,1,0.8,0.7,0.7,0.7,0.7;0.6,0.65,0.8,1,0.85,0.8,0.75,0.75;...
    0.7,0.5,0.7,0.85,1,0.6,0.8,0.65;0.5,0.6,0.7,0.8,0.6,1,0.5,0.7;0.7,0.5,0.7,0.75,0.8,0.5,1,0.8;...
    0.6,0.65,0.7,0.75,0.65,0.7,0.8,1];
b = [0.22,0.19,0.17,0.13,0.11,0.08,0.06,0.04];
%%
MDs = [3.56,7.48,6.54,10.23,2.40,2.30,9.12,7.96];
DTs = [103,155,75,796,89,45,320,245];
ModDur = @(x) MDs*x';ModDur_bench = @(x) ModDur(x-b); Dur_TS_bench  = @(x)(x-b); 
DurTS = @(x) DTs*x';
cov_mat = (sigs'.*sigs).*ro;
carb_intens_func = @(x) carb_intens_123*x';
R = 0.5;
sectors = [1,2,1,1,2,1,2,2];pos1 = find(sectors==1);pos2 = find(sectors==2);
risk_AS = @(x)1/2*(x-b)*(x-b)';
risk_MD = @(x)1/2*(sum((MDs(pos1)*(x(pos1)-b(pos1))').^2)+sum((MDs(pos2)*(x(pos2)-b(pos2))').^2));
risk_DTS = @(x) 1/2*(sum((DTs(pos1)*(x(pos1)-b(pos1))').^2)+sum((DTs(pos2)*(x(pos2)-b(pos2))').^2));
bench_MD = ModDur(b);bench_DTS = DurTS(b);CI_bench = carb_intens_func(b);
%%
func1 = risk_AS;
func2 = @(x) 100*risk_AS(x)+50*risk_MD(x);
func3 = @(x) 100*risk_AS(x)+50*risk_MD(x)+risk_DTS(x);
%%
% Problem 1:
% Minimize the benchmark standard distance by reducing the carbon intensity
% with 50%, while keeping the duration and duration times spread measures
Aeq = [MDs;DTs;ones(1,8)];beq = [bench_MD;bench_DTS;1]; 
A = carb_intens_123;b = (1-R)*CI_bench;
sol1 = fmincon(risk_AS,ones(1,8)/8,A,b,Aeq,beq,zeros(1,8),ones(1,8));
%%
% Problem 2:
% Minimize the aggregate function func2;
Aeq = [DTs;ones(1,8)];beq = [bench_DTS;1];
sol2 = fmincon(func2,ones(1,8)/8,A,b,Aeq,beq,zeros(1,8),ones(1,8));
%%
% Minimize the aggregate function func3;
Aeq = ones(1,8);beq = 1;
sol3 = fmincon(func3,ones(1,8)/8,A,b,Aeq,beq,zeros(1,8),ones(1,8));
