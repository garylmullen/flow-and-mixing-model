clear; close all; clc;
    
%addpath(genpath('/home/gary/Documents/Simulations/'))
%% SET MODEL PARAMETERS

runID = 'demo';  % run identifier tag
nop   = 10;       % print output every 'nop' steps
lvplt = 1;        % plot figures live (1) or in background (0)     
svfig = 1;        % save figures to file (1)

% set domain parameters
N     = 300;      % num. grid size
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
finit      = 'linear'; % initial condition: 'linear' or 'layer'
f0         = 0.05;    % background porosity [vol]
f1         = 0.01;    % base porosity [vol]  
df         = 0.01;    % perturbation amplitude [vol]

Tinit = 'linear';  % initial condition: 'linear' or 'layer'
T0    = 20;        % top temperature [C]
T1    = 100;       % base temperature [C]
dT    = 0.1;       % perturbation amplitude [C]
Cinit = 'linear';  % initial condition: 'linear' or 'layer'
C0    = 0.001;     % top concentration  [C]
C1    = 0.001;     % base concentration [C]
dC    = 0.0;       % perturbation amplitude [C]
zlay  = 0.5;       % relative depth of layer boundary
wlay  = 0.02;      % relative width of layer boundary
xstruct = [D/2,D/2];  % position of structures
zstruct = [400,400];  % onset depth of structures
dstruct = [D,500];  % final depth of structures
wstruct = [D/20,D/5]; % thickness of structures
astruct = [30,90];    % angle of structures to vertical
fstruct = [0.1,0.1];  % porosity of structures (nan = do not set)
Tstruct = [nan,nan];  % temperature of structures (nan = do not set)
Cstruct = [nan,0.05]; % salinity of structures (nan = do not set)
smth    = (N/40)^2;  % smoothness of initial fields

% set model timing parameters
tend  = 1e11;     % model stopping time [s]

% set numerical solver parameters
CFL   = 0.25;     % Courant number to limit time step size
tol   = 1e-7;     % residual tolerance for iterative solver
maxit = 5e3;      % maximum number of iterations
alpha = 0.99;     % step size for iterative solver
beta  = 0.95;     % damping parameter for iterative solver

% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% run code
addpath ../src
main
