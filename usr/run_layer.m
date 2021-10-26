clear; close all; clc;
    
%addpath(genpath('/home/gary/Documents/Simulations/'))
%% SET MODEL PARAMETERS

runID = 'linear_19Oct';  % run identifier tag
nop   = 50;       % print output every 'nop' steps
lvplt = 0;        % plot figures live (1) or in background (0)     
svfig = 1;        % save figures to file (1)

% set domain parameters
N     = 200;      % num. grid size
D     = 1e3;      % phys. domain depth [m]

% set physical parameters

mu    = 1e-3;     % pore fluid viscosity (water) [Pa s]
a     = 5e-3;     % grain size of matrix (sandstone) [m]
b     = 100;      % geom. factor for permeability [1]
n     = 3;        % permeability powerlaw [1]
rhol0 = 1000;     % fluid density [kg/m3]
grav  = 9.81;     % gravity [m/s2]
kC    = 1e-8;     % chemical diffusivity [m2/s]  
kT    = 1e-6;     % thermal diffusivity [m2/s]
aT    = 2e-4;     % thermal expansivity [1/K]
gC    = 1;        % chemical density contrast [1/wt]

% set initial condition parameters
finit      = 'layer'; % initial condition: 'linear' or 'layer'
f0         = 0.05;      % background porosity [vol]
f1         = 0.01;      % base porosity [vol]  
df         = 0.01;    % perturbation amplitude [vol]
f_Layer    = 0.1;
LayerDepth = 500;
LayerWidth = 100;
f_Fault    = 0.1;
FaultDepth = 400;
FaultPos   = 500;
FaultAngle = 30;
FaultWidth = 50;

Tinit = 'linear';  % initial condition: 'linear' or 'layer'
T0    = 20;       % top temperature [C]
T1    = 60;      % base temperature [C]
dT    = 0.1;      % perturbation amplitude [C]
Cinit = 'layer';  % initial condition: 'linear' or 'layer'
C0    = 0.5;    % top concentration  [C]
C1    = 0.1;     % base concentration [C]
dC    = 0.1;      % perturbation amplitude [C]
zlay  = 0.5;      % relative depth of layer boundary
smth  = (N/30)^2; % smoothness of random noise

% set model timing parameters
tend  = 1e11;     % model stopping time [s]

% set numerical solver parameters
CFL   = 0.25;     % Courant number to limit time step size
tol   = 1e-5;     % residual tolerance for iterative solver
alpha = 0.99;     % step size for iterative solver
beta  = 0.95;     % damping parameter for iterative solver

% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% run code
addpath ../src
main
