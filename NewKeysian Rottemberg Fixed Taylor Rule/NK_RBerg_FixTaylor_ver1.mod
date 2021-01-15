% This model does the work of Professor Richard Dennis's LN1 log-linear New Keysian Model
% with a Taylor Rule pegged on a Fixed Long-run Output (the steady state level, I will approximately take it as 0.98).
% The reason for this simplification on Taylor Rule is to show an independent NK model, avoiding introducing the flex-path model, 
% which is required in a standard Taylor Rule implicitly by Y_tilda.
% Most of the notation will be consistent with the Lecture Note; 
% for thoes that are not, explanation will be included to make it as readable as possible;

'Simulation for New Keynesian Rottemberg Model: with Fixed Long-run Output as Benchmark in Taylor Rule'
Author = 'Brian Wang'

var A V Y C H w MC R p; 

% A: Aggregate Productivity;
% Y: Aggregate Output;
% C: Aggregate Consumption;
% H: Aggregate Labour;
% w: real wage;
% MC: real marginal cost of intermediate good producer;
% R: Gross Nominal Interest Rate (R = 1+r), which is better for log-linearization for interpretation as percentage change;
% p: Gross Inflation Level: p = 1 + pi, which is better for log-linearization for more fundamential reason as steady pi = 0;

varexo e v;

% e: Aggregate Productivity Shock;
% v: Cost-push shock in monetary policy;

parameters bet sig n m phi rho eps psi s;

bet     = 0.99;     % household discount factor
sig     = 0.5;        % CRRA, 0.5, 1, 2, 3, sig > 0
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

%  Computers are modest. They only report what they see. If i magnify the shocl to such a degree, it can see nothing, so it says nothing! Rubbish in, rubbsih out; gold in, gold out.

model;

    H = Y/A;                                                    % 1
    MC = w/A;                                                   % 2
    w = H^phi * C^sig;                                          % 3
    C^(-sig) = bet * (C(+1)^(-sig) * R/p(+1));                  % 4
    C = Y - (psi/2) * ((p-1)^2) * Y;                            % 5
    (p-1)*p + ((1-eps)/2)*(p-1)^2 = ((1-eps)/psi) + (eps/psi) * MC + bet * (C(+1)^(-sig)/C^(-sig)) * (Y(+1)/Y) * p(+1) * (p(+1)-1);    % 6
    R = (1/bet) * p^n * (Y/0.98)^m * exp(V);                   % 7
    A = A(-1)^rho * exp(s*e);                                   % 8
    V = s*v;                                                    % 9

end;
 
%  1 for the demand for labour given the production level, as the firm side FOC. See LN1 eq(20);
%  2 for optimized marginal cost as the firm side FOC. This in turn decides the total production level, so interpreted as supply for good. See LN1 eq(21);
%  3 for labour supply as the household side FOC. See LN1 eq(5),eq(6);
%  4 for Euler Equation on intertemporal choice based on nominal bond. See LN1 eq(5),eq(7);
%  5 for good demand from household (or good market constraint). See LN1 eq(41);
%  6 for Phillip Curve (or money market clearance). See LN1 eq(29);
%  7 for Exponential form of Taylor Rule, in which benchmark output is set to be Long Run output Y_st.st = 0.98;  
%  8 for technology process, with reduced size of the shock;
%  9 to show monetary policy shock, with reduced size;

initval;

A   = 1;
V   = 0;
Y   = 1;        %  should expect not 1, as we always have monoplistic competition
C   = 1;
H   = 1;
w   = 8/9;
MC  = 8/9;
p   = 1;        % 1+pi
R   = 1.3;      % R=1+r
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
dynasave('simudata_NK_Rberg_FixTaylor.mat');
