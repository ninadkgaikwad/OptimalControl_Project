function ExamplesSetup

% ------------------%
% Turn off warnings %
% ------------------%
warning off

% --------------------------%
% Get the current directory %
% --------------------------%
currdir  = pwd;
snoptdir = strcat([currdir,'/snopt2017/']);
ipoptdir = strcat([currdir,'/ipopt/']);

% ----------------%
% Set up Adigator %
% ----------------%
adigatorSetup(currdir);

% ------------------------------------------%
% Add NLP Solver SNOPT/IPOPT to MATLAB Path %
% ------------------------------------------%
addpath(snoptdir);
addpath(ipoptdir);

% ---------------------------------%
% Add Directory with LGR Functions %
% ---------------------------------%
RadauDirectory = strcat(currdir,'/collocation');
addpath(RadauDirectory);

% -----------------------%
% Turn on warnings again %
% -----------------------%
warning on
