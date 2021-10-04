function adigatorSetup(currentDirectory);

% --------------------------------------------------------%
% This function sets up adigator for use with SNOPT/IPOPT %
% --------------------------------------------------------%

addpath(strcat(currentDirectory,'/adigator/'));
addpath(strcat(currentDirectory,'/adigator/lib/cadaUtils'));
adigatorHomeDirectory = strcat(currentDirectory,'/adigator/');
dirs = {'doc','examples','lib','util'};
for i=1:length(dirs)
  currdir = fullfile(adigatorHomeDirectory,dirs{i});
  addpath(currdir);
end
savepath;
