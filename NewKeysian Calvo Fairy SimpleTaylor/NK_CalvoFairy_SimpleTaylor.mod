'SIMULATION FOR NEW KEYSIAN CALVO-FAIRY MODEL: SIMPLE TAYLOR RULE'
'Brian Wang'
% =========================================================================
% DESCRIPTION:
% This file contains New Keysian Calvo model, following Prof.Dennis's LN1
%
% Simple Taylor Rule: CB sets nominal interest rate with a trade-off 
% between deviation of output from its steady state level and inflation. 
        % In math: i_t = rho + n * pi + m * y_hat_t + vt, in which 
        % i: nominal interest rate
        % rho=-ln(beta): natural rate of interest (at steady state)
        % n: elasticity of inflation 
        % pi: inflation level
        % m: elasticity of output gap
        % y_hat=ln(Y)-ln(Y_stst): log deviation of output from steady state
        % v: exogenous monetary policy shock that evolves following AR(1)
%
% -------------------------------------------------------------------------
% REFERENCES:
% 1.Professor Richard Dennis Lecture Note 1
% 2.Monetary policy, inflation, and the business cycle : an introduction to
%   the new Keynesian framework and its applications / Jordi Galí, Ch.3 
%
% =========================================================================
% MODEL: 

% Claim all varibles:
var A V Y C H w MC R Pi Pi_star; 

% A: Homogenous Productivity;
% V: Monetary Shock
% Y: Aggregate Output Index / Final Good Production;
% C: Aggregate Consumption Index / Final Good Consumption;
% H: Aggregate Labour; 
    % in CFairy model H is the integral of intermediate firm labour demand
    % as they are ex-post hetereogenous. Meanwhile, it represents each
    % homogenous household's labour supply
% w: Real Wage;
% MC: Real Marginal Cost of intermediate good producer;
    % 1/MC is mark-up
% R: Gross Nominal Interest Rate (R = 1+i)
    %  better interpretable as percentage change under log-linearization
% Pi: Price Level in total --> Gross Inflation Pi=(1+pi)
% Pi_star: Firms Optimal Relative Pricing 
    % Pi_star=P_star(t)/P(t), represents the optimal price by firms tho are 
    % chosen by the Fairy, in other words, the monopoly price
% Notice: We need to have prices all in gross ratio because in New Keysian 
% models prices themselves are undetermined, while the ratios are. 

varexo e v;

% e: Aggregate Productivity Shock;
% v: Cost-push shock in monetary policy;


% Callibrate all parameters:
parameters bet sig n m phi rho eps thet s; 

bet     = 0.99;     % Household Discount Factor
sig     = 1;        % CRRA, 1 = log utility
phi     = 5;        % Labour Intensity 
rho     = 0.9;      % External Persistency
eps     = 9;        % Dixit Stiglitz Elasticity of Substitution
thet    = 0.75;     % (1-thet is) The Smile of Gardevoir (霜奶仙的微笑) ~_~
n       = 1.5;      % Elasticity of Inflation in Talor Rule
m       = 0.5/4;    % Elasticity of Output Gap in Taylor Rule
s       = 0.05;     % Size (and direction if negative) of Shocks  

% Notice:
% 1. sig (CRRA) can be crucial in deciding the correlation between H and Y: 
    % If sig<1, correlation (H,Y) usually is positive:
    % Because consumers are MORE elastic in intertemporal consumption, so 
    % they produce and consume more with a positive productivity shock;
    % If sig>1, correlation (H,Y) usually is negative:
    % Because consumers are LESS elastic in intertemporal consumption, so 
    % they would rather enjoy more leisure at a positive productivity shock
% 2. Computers are modest. They only report what they see. Therefore, keep
    % s small (0.05). If I magnify the shock to such a degree that it sees 
    % nothing, it says nothing. Quote from CompSci: Rubbish in, rubbsih out


% Calvo Fairy Model:
model;

    H = Y/A * (thet + (1-thet)*(Pi_star*Pi)^(-eps))/Pi^(-eps);  % 1  -H
    MC = w/A;                                                   % 2  -MC
    log(Pi_star) = bet*thet*log(Pi_star(+1)*Pi(+1))+(1-bet*thet)*(-(log(1/MC)-log(eps/(eps-1)))); % 3 -Pi_Star
    w = H^phi * C^sig;                                          % 4  -w
    C^(-sig) = bet * (C(+1)^(-sig) * R/Pi(+1));                 % 5  -C
    C = Y;                                                      % 6  -Y
    Pi^(1-eps)=(1-thet)*(Pi_star*Pi)^(1-eps)+thet;              % 7  -Pi
    R = (1/bet) * Pi^n * (Y/0.983)^m * exp(V);                   % 8  -R
    A = A(-1)^rho * exp(s*e);                                   % 9  -A
    V = s*v;                                                    % 10 -V

end;
 
%  1 for the Labour Demand as the firm side FOC. - LN1 eq(37) ;
%  2 for optimized marginal cost as the firm side FOC. - LN1 eq(38);
   % This decides the production level, so interpreted as Goods Supply;
%  3 for Fairy's Impact on lucky firms. - LN1 eq(32) (可爱就是正义) ~_~；
%  4 for Labour Supply as the household side FOC. - LN1 eq(5),eq(6);
%  5 for Goods Demand from Euler Equation on nominal bonds - LN1 eq(7);
%  6 for Goods Market Clear-out - LN1 eq(43);
%  7 for Price Movement in CFairy - LN1 eq(29);
%  7 for Exponential form of Taylor Rule
   % in which Long Run output Y_stst = 0.983;  
%  8 for technology process, with reduced size of the shock;
%  9 to show monetary policy shock, with reduced size;


% Solve for Steady State Values:
initval;

A       = 1;
V       = 0;
Y       = 1;        %  should<1, as we always have monoplistic competition
C       = 1;
H       = 1;
w       = 8/9;
MC      = 8/9;
R       = 1.3;      % R=1+i,nominal
Pi      = 1;
Pi_star = 1;        % at st.st optimal price = aggregate price
e       = 0;
v       = 0;

end;

shocks;
var e;
stderr 1;
var v;
stderr 1;
end;

steady;
check;


% Sotchastic simulation to generate IRF:
stoch_simul(periods=140,drop=100,irf=40,order=1);
dynasave('simudata_NK_Rberg_FixTaylor.mat');

% Notes:
% periods: how many period we want the whole simulation take;
% drop: periods we want to drop before reaching steady state;
% irf: periods we want to show impulse response over;
% order: by default, Dynare uses 2nd order appx. It generates a stochastic 
% result and the graphs actually shows average IRF, which may be quirky. 
% To see a smooth graph and the direct effect of shocks, use 1st order appx
% or reduce the size of the shock even smaller to s=0.01.

