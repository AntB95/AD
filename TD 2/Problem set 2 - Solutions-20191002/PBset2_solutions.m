%--------------------------------------------------------------------------
% Advanced Derivatives 
% Problem Set 2
% Marc-Aurèle Divernois
%
% This code calls BlackScholes.m
%--------------------------------------------------------------------------
% Exercice 1
%--------------------------------------------------------------------------
clear all;
clc;

% Statements
V_0 = 100;
sigma_0 = 0.4;
D = 50;
T = 10;
r = 0.05;
t=0;
q=0;

% Current value of the firm's equity
E = BlackScholes(V_0,t,sigma_0,D,T,r,q);
strikes = [0.60; 0.80; 1; 1.2]*E;
maturities = [2;5;7;9];

% Find V_star
V_star = NaN(size(strikes,1),size(maturities,1));

for i = 1:length(strikes)
    for j = 1:length(maturities)
        fun = @(V) (BlackScholes(V,maturities(j),sigma_0,D,T,r,q) - strikes(i))^2;
        V_star(i,j) = fminsearch(fun,V_0);
    end
end

% Prices of compound options

a_1 = (log(V_0./V_star) + (r+sigma_0^2/2)*(maturities -t)')...
    ./(sigma_0*sqrt(maturities-t))';
a_2 = a_1 - sigma_0*sqrt(maturities-t)';
b_1 = (log(V_0/D) + (r+sigma_0^2/2)*(T -t))/(sigma_0*sqrt(T-t));
b_2 = b_1 - sigma_0*sqrt(T-t);
rho = sqrt((maturities-t)/(T-t));
mu = zeros(2,1);

C = NaN(size(V_star));

for i = 1:length(strikes)
    for j = 1:length(maturities)
        covmat = [1 rho(j);rho(j) 1];    
        C(i,j) = V_0 * mvncdf([a_1(i,j); b_1],mu,covmat) - ... 
            D * exp(-r*(T-t)) * mvncdf([a_2(i,j); b_2],mu,covmat) - ...
            exp(-r*(maturities(j)-t)) * strikes(i) * normcdf(a_2(i,j));
    end
end

% Black-Scholes implied volatilities of firm's equity

implied_vol = NaN(size(C));

for i = 1:length(strikes)
    for j = 1:length(maturities)
        fun2 = @(sigma) (BlackScholes(E,t,sigma,strikes(i),maturities(j),r,q) - C(i,j))^2;
        [implied_vol(i,j),fval(i,j)] = fminsearch(fun2,0.52);
    end
end

% Volatility surface

x = repmat(strikes,1,4);
y = repmat(maturities,1,4)';

figure
surf(x,y,implied_vol)
title('Volatility surface')
xlabel('Strikes')
ylabel('Maturities')
zlabel('Volatility')

%--------------------------------------------------------------------------
% Exercice 2
%--------------------------------------------------------------------------
clear all;
clc;

% Statements
S = 100;
gamma = -0.08;
lambda = 0.2;
sigma = 0.2;
r = 0.04;
K = [0.8;0.9;1;1.1]*S;
T=[0.02;0.08;0.25;0.5];
t=0;

n=4; % We can enter any number here but it converges very fast...
C = zeros(size(K,1),size(T,1));

% Compute option prices

for i=1:length(K)
    for j=1:length(T)
        for k=1:n
            C_temp = exp(-lambda*T(j))*((lambda*T(j))^(k-1)/(factorial(k-1))) * ...
            BlackScholes(S*(1+gamma)^(k-1),t,sigma,K(i),T(j),r,lambda*gamma);
            C(i,j) = C(i,j) + C_temp;
        end
    end
end

% Compute option volatilities 

implied_vol2 = NaN(size(C));

for i = 1:length(K)
    for j = 1:length(T)
        fun3 = @(sigma) (BlackScholes(S,t,sigma,K(i),T(j),r,0) - C(i,j))^2;
        [implied_vol2(i,j),fval(i,j)] = fminsearch(fun3,sigma);
    end
end

% Volatility surface

x = repmat(K,1,4);
y = repmat(T,1,4)';

figure
surf(x,y,implied_vol2)
title('Volatility surface')
xlabel('Strikes')
ylabel('Maturities')
zlabel('Volatility')

% Plots

figure
subplot(2,2,1)
plot(x(:,1),implied_vol2(:,1))
title('Maturity : 0.02')
xlabel('Strikes')
ylabel('Volatility')
subplot(2,2,2)
plot(x(:,2),implied_vol2(:,2))
title('Maturity : 0.08')
xlabel('Strikes')
ylabel('Volatility')
subplot(2,2,3)
plot(x(:,3),implied_vol2(:,3))
title('Maturity : 0.25')
xlabel('Strikes')
ylabel('Volatility')
subplot(2,2,4)
plot(x(:,4),implied_vol2(:,4))
title('Maturity : 0.5')
xlabel('Strikes')
ylabel('Volatility')

