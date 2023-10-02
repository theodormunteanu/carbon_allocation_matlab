function [str,DR,RCs,pos] = MDP_str_constraint(sigs,ro,scores,bench,alpha,CIs,k,R)
   % Maximum diversification portfolio when scores (credit, ESG / equity )
   % are given
   % CIs = Carbon intensities / emissions
   n = length(sigs);
   cov_mat = (sigs'.*sigs).*ro;
   f = @(x) -sqrt(x*cov_mat*x')/(x*sigs');
   if isempty(scores)==0
       pos = arrayfun(@(i) find(scores(i)==sort(scores,'asc')),1:k);
       z = zeros(1,n);z(pos) = ones(1,k);
   end
   if isempty(CIs)==0 && isempty(scores)==0
      A = [z;CIs];b = [alpha;(1-R)*(bench*CIs')];
   elseif isempty(CIs)==1 && isempty(scores)==0
      A = z;b = alpha;
   elseif isempty(CIs)==0 && isempty(scores)==1
      A = CIs;b = (1-R)*bench*CIs';
   else
      A = [];b = [];
   end
   Aeq = ones(1,n);beq = 1;
   x0 = ones(1,n)/n;
   str = fmincon(@(x)-f(x),x0,A,b,Aeq,beq,zeros(1,n),ones(1,n));
   DR = -f(str);
   RCs = str.*(cov_mat*str')'/(str*cov_mat*str');
end