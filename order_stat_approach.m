min_DI = @minimize_DI;
min_TE = @minimize_track_error;
order_stat_CA = @order_statistic_carb_appr;
%%
function allocs = minimize_DI(CE1,CE2,CE3,sigs,ro,Y,b,R,HCIS,Mkt_Cap,rets)
%%
% Minimize diversification index when Carbon Emissions are given
% CE1,CE2,CE3 represent the three tiers of carbon Emissions
% HCIS = index of dummy variables for High Carbon Intensity Sectors 
%        [0/1]: 0 for low intensity, 1 for high intensity
% MktCap: Market Capitalizations

% Allocs: represents the vector of percentages for allocations 
%%
   allocs = struct();
   cov_mat = (sigs'.*sigs).*ro;
   mkt_caps = b*Mkt_Cap;
   carb_intens1 = CE1./Y;carb_intens2 = CE2./Y;
    carb_intens3 = CE3./Y;
    carb_intens_12 = carb_intens1+carb_intens2;
    carb_intens_123 = carb_intens1+carb_intens2+carb_intens3;
    C1 = carb_intens1;D1 = (1-R)*b*carb_intens1';
    C12 = carb_intens_12;D12 = (1-R)*b*carb_intens_12';
    C123 = carb_intens_123;D123 = (1-R)*b*carb_intens_123';
    DR = @(x) x*sigs'/sqrt(x*cov_mat*x');
    n = length(sigs);
    HCIS_bench_exposure = b*HCIS';
    str1 = fmincon(DR,ones(1,n)/n,C1,D1,ones(1,n),1,zeros(1,n),ones(1,n));
    str12 = fmincon(DR,ones(1,n)/n,C12,D12,ones(1,n),1,zeros(1,n),ones(1,n));
    str123 = fmincon(DR,ones(1,n)/n,C123,D123,ones(1,n),1,zeros(1,n),ones(1,n));
    allocs.str1 = str1;allocs.str12 = str12;allocs.str123 = str123;
    C_HCIS_12 = [C12;HCIS];C_HCIS_123 = [C123;HCIS];
    D_HCIS_12 = [D12;-HCIS_bench_exposure];D_HCIS_123 = [D123;-HCIS_bench_exposure];
    str12_HCIS = fmincon(DR,ones(1,n)/n,C_HCIS_12,D_HCIS_12,ones(1,n),1,zeros(1,n),ones(1,n));
    str123_HCIS = fmincon(DR,ones(1,n)/n,C_HCIS_123,D_HCIS_123,ones(1,n),1,zeros(1,n),ones(1,n));
    Carb_em_12 = (CE1+CE2)./mkt_caps;D_em_12 = b*(1-R)*((CE1+CE2)./mkt_caps)';
    str12_em = fmincon(DR,ones(1,n)/n,Carb_em_12,D_em_12,ones(1,n),1,zeros(1,n),ones(1,n));
    Carb_em_123 = (CE1+CE2)./mkt_caps;D_em_123 = b*(1-R)*((CE1+CE2+CE3)./mkt_caps)';
    str123_em = fmincon(DR,ones(1,n)/n,Carb_em_123,D_em_123,ones(1,n),1,zeros(1,n),ones(1,n));
    allocs.str12_HCIS = str12_HCIS;
    allocs.str123_HCIS = str123_HCIS;
    allocs.str12_em = str12_em;
    allocs.str123_em = str123_em;
end
%%
function allocs = minimize_track_error(CE1,CE2,CE3,cov_mat,Y,b,R,HCIS,Mkt_Cap,rets)
    % HCIS = binary vector of 0 and 1: 1 for belonging to HCIS and 0 for
    % not belonging to HCIS. 
    % b = structure of the benchmark portfolio
    
    allocs = struct();
    mkt_caps = b*Mkt_Cap;
    carb_intens1 = CE1./Y;carb_intens2 = CE2./Y;
    carb_intens3 = CE3./Yz;
    carb_intens_12 = carb_intens1+carb_intens2;
    carb_intens_123 = carb_intens1+carb_intens2+carb_intens3;
    C1 = carb_intens1;D1 = (1-R)*b*carb_intens1';
    C12 = carb_intens_12;D12 = (1-R)*b*carb_intens_12';
    C123 = carb_intens_123;D123 = (1-R)*b*carb_intens_123';
    func = @(x) (x-b)*cov_mat*(x-b)';
    n = length(b);
    str1 = fmincon(func,ones(1,n)/n,C1,D1,ones(1,n),1,zeros(1,n),ones(1,n));
    str12 = fmincon(func,ones(1,n)/n,C12,D12,ones(1,n),1,zeros(1,n),ones(1,n));
    str123 = fmincon(func,ones(1,n)/n,C123,D123,ones(1,n),1,zeros(1,n),ones(1,n));
    HCIS_bench_exposure = b*HCIS';
    C_HCIS_12 = [C12;HCIS];
    D_HCIS_12 = [D12;-HCIS_bench_exposure];
    str12_HCIS = fmincon(func,ones(1,n)/n,C_HCIS_12,D_HCIS_12,ones(1,n),1,zeros(1,n),ones(1,n));
    Carb_em_12 = (CE1+CE2)./mkt_caps;D_em_12 = b*(1-R)*((CE1+CE2)./mkt_caps)';
    str12_em = fmincon(func,ones(1,n)/n,Carb_em_12,D_em_12,ones(1,n),1,zeros(1,n),ones(1,n));
    Carb_em_123 = (CE1+CE2)./mkt_caps;D_em_123 = b*(1-R)*((CE1+CE2+CE3)./mkt_caps)';
    str123_em = fmincon(func,ones(1,n)/n,Carb_em_123,D_em_123,ones(1,n),1,zeros(1,n),ones(1,n));
    allocs.str1 = str1;allocs.str12 = str12;allocs.str123 = str123;
    allocs.str12_HCIS = str12_HCIS;allocs.str123_em = str123_em;
    allocs.perf1 = str1*rets';
    allocs.perf12 = str12*rets';
    allocs.perf123 = str123*rets';
end
%%
function str = order_statistic_carb_appr(CE1,CE2,CE3,b,Y,R,sigs,ro)
   cov_mat = covariance(sigs,ro);
   CI1 = CE1./Y;CI2 = CE2./Y;CI3 = CE3./Y;
   n = length(sigs);
   sorted_carb_intens = sort(CI1+CI2+CI3);
   pos_k = arrayfun(@(i) find(CI1+CI2+CI3 == sorted_carb_intens(i)),1:k);
   func = @(x) (x-b)*cov_mat*(x-b)';
   z = zeros(k,n);
   for i = 1:k
       z(i,pos_k(i))=1;
   end
   Aeq = [ones(1,n);z];beq = [1;zeros(k,1)];
   C12 = CI1+CI2;D12 = (1-R)*b*C12';
   str = fmincon(func,ones(1,n)/n,C12,D12,Aeq,beq,zeros(1,n),ones(1,n));
end
%%
function alloc = maximize_performance(CE1,CE2,CE3,cov_mat,Y,b,R,HCIS,Mkt_Cap,rets)
   
end