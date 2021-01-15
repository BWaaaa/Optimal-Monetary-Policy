% This model does the work of Professor Richard Dennis's LN1 Page 18 3.5.1,  
% In fact, this model is the Flex-Path model where the burning cost coefficient psi = 0.
% Most of the notation will be consistent with the Lecture Note; 
% for thoes that are not, explanation will be included to make it as readable as possible;

'Simulation for NK Model: without Price Rigidity'
Author = 'Brian Wang'

var A Y C H w MC Rp; 

% A: Aggregate Productivity;
% Y: Aggregate Output;
% C: Aggregate Consumption;
% H: Aggregate Labour;
% w: real wage;
% MC: real marginal cost of intermediate good producer;
% Rp: denote the ratio between Gross Nominal Interest and Gross Inflation, = R/p = (1+r)/(1+pi) ; 

% The reason we do not setup R (Gross Nominal Interest) and p (Gross Inflation) is because 
% If we include both, this basic model will have more unknown than equations. 
% The intuition is that Rp as a whole represents the intertemporal price of consumption, as nominal bonds are the only saving vehicle.
% In this basic model, prices are not effecting prodcution decisions.
% As to consumer decisions, price level itself can be at random in each period. Consumer will adjust comptemporarily, as long as the ratio Rp is given. 
 
varexo e;

% Exogenous shock to aggregate productivity;

parameters bet sig phi rho eps s;

bet     = 0.99;         % household discount factor
sig     = 2;            % CRRA, 0.5, 1, 2, 3, sig > 0
phi     = 5;            % Labour Intensity 
eps     = 9;            % Dixit Stiglitz elasticity of substitution
rho     = 0.9;          % External Persistency
% psi   = 0;            % Burning waste for price adjustment 
s       = 0.05;         % Size-adjustment of the technology shock

% sig (CRRA) can be crucial in deciding the correlation between H(labour) and Y(output). 
% If sig < 1, correlation (H,Y) usually is positive:
% This is explained by the fact that consumer are MORE elastic in intertemporal consumption, so they will produce and consume more with a positive productivity shock;
% If sig > 1, correlation(H,Y) usually is negative:
% This is explained by the fact that consumer is LESS elastic in intertemporal consumption, so they would rather enjoy more leisure with a positive productivity shock;
% In this specific model, when sig = 1, H is a constant process. 

% The s is a size adjustment to the productivity shock, which is CRUCIAL to this model to be around 0.05. 
% We will DESTROY the economy if we let the shock to have variance 1. Noted the steady state output is only 1, by such a great shock the economic will never go back to norm. 

model; 
    H = Y/A;                                                    % 1
    MC = w/A;                                                   % 2
    w = H^phi * C^sig;                                          % 3
    C^(-sig) = bet * C(+1)^(-sig) * Rp;% /p ;                   % 4 
    C = Y;                                                      % 5
    (1-eps) + (eps) * MC = 0;                                   % 6 
    A = A(-1)^rho * exp(s*e);                                   % 7 
end;

%  1 for the demand for labour given the production level, as the firm side FOC. See LN1 eq(20);
%  2 for optimized marginal cost as the firm side FOC. This in turn decides the total production level, so interpreted as supply for good. See LN1 eq(21);
%  3 for labour supply as the household side FOC. See LN1 eq(5),eq(6);
%  4 for Euler Equation on intertemporal choice based on nominal bond. See LN1 eq(5),eq(7);
%  5 for good demand from household (or good market constraint). See LN1 eq(41);
%  6 for 'Phillip Curve' when psi = 0. This is to say, without money stickiness, the relative mark-up for a intermediate good producer is fixed wrt Dixit-Stiglitiz substitution effect. See LN1 3.5.1;
%  7 for technology process, with reduced size of the shock.

initval;
A   = 1;
Y   = 1;        % should expect < 1, as we always have monoplistic competition
C   = 1;
H   = 1;
w   = 8/9;
MC  = 8/9;
Rp  = 1.3;      
e   = 0;
end;

shocks;
var e;
stderr 1;
end;

steady;
check;
stoch_simul(order = 2, periods=1000,irf=100);   % By default, let dynare do the approximateion in second order approximateion
dynasave('simudata_NK_Basics_ver1.mat');
