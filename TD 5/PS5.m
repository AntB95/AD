clear all;
clc;
%read data
data = xlsread('SX5E_Impliedvols.xlsx');
% take the date column starting at T = 0 
T = [0 data(1,2:end)];
% take data vol
data_vol = data(2:end,2:end);
%parameters
t=0;
r=0;
q=0;
S0=2772.70;
% take the strike from the data
K = data(2:end,1)*S0;

%creating and initialisation of the price matrix
NK = length(K);
NT = length(T);
%creation of the price matrix 
Cobs = zeros(NK,NT);
%initialisation
C0 = max(S0-K,0)';
Cobs(:,1) = C0;

%implementation of CObs matrix
for i=1:NK
    for j =1:NT-1
        %if we have data in the vol matrix we complete the price
        if (data_vol(i,j) >0)
            Cobs(i,j+1) = BlackScholes(S0,t,data_vol(i,j),K(i),T(j+1),r,q);
        end    
    end
end

% true or false matrix
%1 value in the matrix if the vol is known 0 if not 
tf_matrix = data_vol>0;

%Optimization pb
%The Andreasen-Huge algorithm
for j = 1:NT-1 
    not_zeros = find((data_vol(:,j)>0));
    %delta t and k
    dt = T(j+1)-T(j);
    dK = K(2)-K(1);
    %solving the optimization problem where Cmodel is obtain by the
    %Andreasen huge algo
    opti = @(x) sum(tf_matrix(:,j).*(AH(tf_matrix,j,x,Cobs,dt,dK,K)-Cobs(:,j+1)).^2);
    %inotialization of the optimization pb
    x0 = 0.2*ones(length(NK),1)*S0;
    %optimization to find the volatility
    [X{j},fval] = fminsearch(opti,x0);
    %New value of the call option using Andrasen Huge algo and the
    %volatiliy we obtain using the optimization line 49
    NewC = AH(tf_matrix,j,X{j},Cobs,dt,dK,K);
    Cobs(:,j+1) = NewC;
end

%generate imply volatility
IV = NaN(size(data_vol));

for i = 1:NK
    for j = 2:NT
        %optimize for the vol using the Cobs obtain with AH algo 
        opti2 = @(sigma) (BlackScholes(S0,t,sigma,K(i),T(j),r,q) - Cobs(i,j))^2;
        [IV(i,j-1),fval(i,j)] = fminsearch(opti2,0.5);
    end
end
% Volatility surface plot

x = repmat(K,1,NT-1);
y = repmat(T(2:end),NK,1);

%Plot
figure
surf(x,y,IV)
title('Volatility surface')
%label
xlabel('Strikes')
ylabel('Maturities')
zlabel('Volatility')

% Main step of the program:
% - initilaisation
% - Stand from T0 and construct the 1 line Ci = max(S0-Ki,0)
% - Construct the matrix of observation using BS model and the strike and the
% volatily given in the data for Tj (we have Nk prices)
% - Use AH algo to find C_model at time Tj and the optimization problem
%sum(C_model - C_obs)^2 to find the NK volatility at time Tj
% - Using this volatility and the AH algo find the new price C at time
% T_j+1 and loop on j 
% - With all the price solve the optimization prob between the BS price and
% the C_obs price to find the IV 
% - Plot the IV surface as a function of strike and maturity 

