%Exercice 3
clear all;

% Parameters
S = 100;
gamma = 0.1;
lambda = 1.0;
sigma = 0.2;
r = 0.04;
K = linspace(60,150,1000)';
T = [0.02;0.08];
t = 0;
q = 0;
n = 100;
% Price of the call option 
C = zeros(size(K,1),size(T,1));
% using the formula for a call and an EDS with a jump process 
for i=1:length(K)
    for j=1:length(T)
        temp = 0;
        for k=1:n
            C(i,j) = C(i,j) + exp(-lambda*T(j))*((lambda*T(j))^(k-1)/(factorial(k-1)))*BlackScholes(S*(1 - gamma)^(k-1),t,sigma,K(i),T(j),r,q-lambda*gamma);
        end
    end
end

% Computing C tild the undiscouted price 
C_t = C.*exp(r*(T'-t));
% We compute phi by differenetation of C 
phi = (C_t(1:end-2,:)-2*C_t(2:end-1,:)+C_t(3:end,:))/(K(2)-K(1))^2;

% We plot the result 
plot(K(2:end-1,:),phi)
xlabel('Strike value')
ylabel('Value for phi')
title('Implied probabilty distribution')
legend('T1=0.02','T2=0.08')
