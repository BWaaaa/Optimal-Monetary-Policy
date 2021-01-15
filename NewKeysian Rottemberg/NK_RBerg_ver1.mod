'Simulation for New Keysian Rottemberg Model'
Author = 'Brian Wang'

var A V Y_til C H w MC R p Y; 
% A: Aggregate Productivity;
% Y_til = Y/Yf;
% C: Aggregate Consumption;
% H: Aggregate Labour;
% w: real wage;
% MC: real marginal cost of intermediate good producer;
% R: Gross Nominal Interest Rate (R = 1+r), which is better for log-linearization for interpretation as percentage change;
% p: Gross Inflation Level: p = 1 + pi, which is better for log-linearization for more fundamential reason as steady pi = 0;
% Y: Aggregate Output (the real output level);

% Y_tilda: defined as logY_tilde = y_tilde = log(Y) - log(Yf), in which Yf is the ideal Aggregate Output under flex-price path, 
% With refer to the flex-path output, please go to document NewKeysian Basic Model/NK_Basics1_ver1.mod. 
% Also, see LN1 3.5.1, y_tilde is the log difference as defined on Page18 between the real output change, 
% and the ideal situation where there is no sticky price distrtion. The monopolistic distortion is maintained throughout.
% Equivilantly we get  Yf = Y/Y_tilde;

varexo e v;

parameters bet sig n m phi rho eps psi s;

bet     = 0.99;     % household discount factor
sig     = 2;        % CRRA, 0.5, 1, 2, 3, sig > 0
phi     = 5;        % Labour Intensity 
rho     = 0.9;      % External Persistency
eps     = 9;        % Dixit Stiglitz elasticity of substitution
psi     = 80;       % Burning waste for price adjustment 
n       = 1.5;      % elasticity of inflation in Talor Rule
m       = 0.5/4;    % elasticity of output gap in Taylor Rule
s       = 0.05;     % Size of the technology shock

% sig (CRRA) can be crucial in deciding the correlation between H(labour) and Y(output). 
% If sig < 1, correlation (H,Y) usually is positive:
% This is explained by the fact that consumer are MORE elastic in intertemporal consumption, so they will produce and consume more with a positive productivity shock;
% If sig > 1, correlation(H,Y) usually is negative:
% This is explained by the fact that consumer is LESS elastic in intertemporal consumption, so they would rather enjoy more leisure with a positive productivity shock;

%  Note: please ensure the chock is mediated by s=0.05. Larger shocks will destroy the convergence because dynare uses stochastic simulation 
%  Intuitively, We will DESTROY the economy if we let the shock to have variance 1. Noted the steady state output is only 1, by such a great shock the economic will never go back to norm. 

model;

    H = Y/A;                                                    % 1
    MC = w/A;                                                   % 2
    w = H^phi * C^sig;                                          % 3
    C^(-sig) = bet * (C(+1)^(-sig) * R/p(+1));                  % 4
    C = Y - (psi/2) * ((p-1)^2) * Y;                            % 5
    (p-1)*p + ((1-eps)/2)*(p-1)^2 = ((1-eps)/psi) + (eps/psi) * MC + bet * (C(+1)^(-sig)/C^(-sig)) * (Y(+1)/Y) * p(+1) * (p(+1)-1);    % 6
    R = (1/bet) * p^n * (Y_til)^m * exp(V);                     % 7
    ((Y/Y_til)/A)^phi * (Y/Y_til)^sig =  (eps-1)/eps * A;       % 8
    A = A(-1)^rho * exp(s*e);                                   % 9
    V = s*v;                                                    % 10
    
end;

%  1 for the demand for labour given the production level, as the firm side FOC. See LN1 eq(20);
%  2 for optimized marginal cost as the firm side FOC. This in turn decides the total production level, so interpreted as supply for good. See LN1 eq(21);
%  3 for labour supply as the household side FOC. See LN1 eq(5),eq(6);
%  4 for Euler Equation on intertemporal choice based on nominal bond. See LN1 eq(5),eq(7);
%  5 for good demand from household (or good market constraint). See LN1 eq(41);
%  6 for Phillip Curve (or money market clearance). See LN1 eq(29);
%  7 for Exponential form of Taylor Rule;
%  8 for output under Flex Price Path, in taylor rule logY_tilda is defined as log(Y)-log(Y_FlexPricePath). 
%  9 for technology process, with reduced size of the shock!
%  10 to show monetary policy shock, with reduced size of the shock!
%  The s is a size adjustment to the productivity shock, which is CRUCIAL to this model to be around 0.05. 

initval;
A   = 1;
V   = 0;
Y   = 1;        %  should expect not 1, as we always have monoplistic competition
Y_til = 1;      %  should expect 1, during steady state inflation is 0, so they have no difference
C   = 1;
H   = 1;
w   = 8/9;
MC  = 8/9;
R   = 1.1;      %  R = 1 + r
p   = 1;        %  p = q + pi 
e   = 0;
v   = 0;
end;

shocks;
var e;
stderr 1;
var v;
stderr 1;
end;

steady;
check;

stoch_simul(periods=1000,irf=100);
%  We will DESTROY the economy if we let the shock to have variance 1. Noted the steady state output is only 1, by such a great shock the economic will never go back to norm. 

dynasave('simudata_NK_RBerg.mat');
