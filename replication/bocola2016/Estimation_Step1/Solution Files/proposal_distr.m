function [obj] = proposal_distr(data,bounds,param,e,initial,gdp_last,policies,Sigma,s)

sim_data   = simul_filter(bounds,param,e,initial,policies,s);

sim_data(1)= (sim_data(1)-gdp_last);

err        = data-sim_data;

obj        = err'*Sigma*err + e(2,:)*e(2,:)';

end

