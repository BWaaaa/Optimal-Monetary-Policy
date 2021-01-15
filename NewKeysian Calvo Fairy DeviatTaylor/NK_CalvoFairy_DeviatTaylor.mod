'SIMULATION FOR NEW KEYSIAN CALVO-FAIRY MODEL: DEVIATION TAYLOR RULE'
'Brian Wang'
% =========================================================================
% DESCRIPTION:
% This file contains New Keysian Calvo model
%
% Deviation Taylor Rule: CB sets nominal interest rate with a trade-off 
% between deviation of output from Flex-Price Path level and inflation. 
% Flex-Price Path: This is when prices are totally flexible, and therefore 
% no inflation and no space for a CB to manipulate through MP. In fact, the 
% New Keysian model here regress back towards a RBC model. 
        % In math: i_t = rho + n * pi + m * y_til_t + vt, in which 
        % i: nominal interest rate
        % rho=-ln(beta): natural rate of interest (at steady state)
        % n: elasticity of inflation 
        % pi: inflation level
        % m: elasticity of output gap
        % y_til=ln(Y)-ln(Y_fp): log deviation of output from fp path value
        % v: exogenous monetary policy shock that evolves following AR(1)
%
% Notice: One important conclusion from theory is that the Flex-Price Path 
% of Calvo Fairy Model is exactly the same as which of Rottemberg Burning 
% Model. This is because, respectively, the Flex-Price Path model equals 
    % 1.NK Rottemberg: psi=0;
    % 2.Calvo Fairy: thet=0;
% 
% -------------------------------------------------------------------------
% REFERENCES:
% 1.Lecture Note 
% 2.Monetary policy, inflation, and the business cycle : an introduction to
%   the new Keynesian framework and its applications / Jordi Galí, Ch.3 
%
% =========================================================================
% MODEL: 

% Claim all varibles:
var A V Y C H w MC R Pi Pi_star 
Y_fp C_fp H_fp w_fp MC_fp R_fp 
Y_til C_til H_til w_til MC_til R_til; 

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
% X_fp: corresponding varaibles under flex-price path
% X_til: the log difference between X and X_fp, absolute deviation
% Notice: 
% 1. We need to have prices all in gross ratio because in New Keysian 
    % models prices themselves are undetermined, while the ratios are. 
% 2. We do not have Pi_fp because in flex-price path they will all be 1

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
thet    = 0.75;     % (1-thet is) The Smile of the Fairy
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
    % nothing, it says nothing. 


% Calvo Fairy Model:
model;

    H = Y/A * (thet + (1-thet)*(Pi_star*Pi)^(-eps))/Pi^(-eps);  % 1  -H
    MC = w/A;                                                   % 2  -MC
    log(Pi_star) = bet*thet*log(Pi_star(+1)*Pi(+1))+(1-bet*thet)*(-(log(1/MC)-log(eps/(eps-1)))); % 3 -Pi_Star
    w = H^phi * C^sig;                                          % 4  -w
    C^(-sig) = bet * (C(+1)^(-sig) * R/Pi(+1));                 % 5  -C
    C = Y;                                                      % 6  -Y
    Pi^(1-eps)=(1-thet)*(Pi_star*Pi)^(1-eps)+thet;              % 7  -Pi
    R = (1/bet) * Pi^n * Y_til^m * exp(V);                      % 8  -R
    A = A(-1)^rho * exp(s*e);                                   % 9  -A
    V = s*v;                                                    % 10 -V

    H_fp = Y_fp/A;                                  % 11
    MC_fp = w_fp/A;                                 % 12
    MC_fp = (eps-1)/eps;                            % 13
    w_fp = H_fp^phi * C_fp^sig;                     % 14
    C_fp^(-sig) = bet * (C_fp(+1)^(-sig) * R_fp/1); % 15
    C_fp = Y_fp; % 16

    Y_til = Y/Y_fp;         % 21
    H_til = H/H_fp;         % 22
    C_til = C/C_fp;         % 23
    w_til = w/w_fp;         % 24
    MC_til = MC/MC_fp;      % 25
    R_til = R/R_fp;         % 26

end;
 
%  1 for the Labour Demand as the firm side FOC. - LN1 eq(37) ;
%  2 for optimized marginal cost as the firm side FOC. - LN1 eq(38);
   % This decides the production level, so interpreted as Goods Supply;
%  3 for Fairy's Impact on lucky firms. - LN1 eq(32);
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

Y_fp    = 0.983;
C_fp    = 0.983;
H_fp    = 1;
w_fp    = 8/9;
MC_fp   = 8/9;
R_fp    = 1.3;

Y_til   = 1;
C_til   = 1;
H_til   = 1;
w_til   = 1;
MC_til  = 1;
R_til   = 1;

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
stoch_simul(periods=140,drop=100,irf=40,order=1) A V Y_til C_til
H_til w_til MC_til R_til Pi Pi_star;
dynasave('simudata_NK_Rberg_FixTaylor.mat');

% Notes:
% periods: how many period we want the whole simulation take;
% drop: periods we want to drop before reaching steady state;
% irf: periods we want to show impulse response over;
% order: by default, Dynare uses 2nd order appx. It generates a stochastic 
% result and the graphs actually shows average IRF, which may be quirky. 
% To see a smooth graph and the direct effect of shocks, use 1st order appx
% or reduce the size of the shock even smaller to s=0.01.

