function payoff = BlackScholes(V,t,sigma,D,T,r,q)
d1 = (log(V/D)+ (r -q + sigma^2 / 2)* (T-t))/(sigma*sqrt(T-t));
d2 = d1 - sigma*sqrt(T-t);
payoff = V * exp(-q*(T-t))* normcdf(d1) - normcdf(d2) * exp(-r*(T-t)) * D;