% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function y = logsumexp_ADiGatorGrd(x)
global ADiGator_logsumexp_ADiGatorGrd
if isempty(ADiGator_logsumexp_ADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_logsumexp_ADiGatorGrd.logsumexp_ADiGatorGrd.Gator1Data;
% ADiGator Start Derivative Computations
cada1f1dx = exp(x.f(:)).*x.dx;
cada1f1 = exp(x.f);
cada1f2dx = cada1f1dx;
cada1f2 = sum(cada1f1);
y.dx = 1./cada1f2.*cada1f2dx;
y.f = log(cada1f2);
%User Line: y = log(sum(exp(x)));
y.dx_size = 32;
y.dx_location = Gator1Data.Index1;
end


function ADiGator_LoadData()
global ADiGator_logsumexp_ADiGatorGrd
ADiGator_logsumexp_ADiGatorGrd = load('logsumexp_ADiGatorGrd.mat');
return
end