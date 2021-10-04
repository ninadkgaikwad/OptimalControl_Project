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

function y = logsumexp_ADiGatorHes(x)
global ADiGator_logsumexp_ADiGatorHes
if isempty(ADiGator_logsumexp_ADiGatorHes); ADiGator_LoadData(); end
Gator1Data = ADiGator_logsumexp_ADiGatorHes.logsumexp_ADiGatorHes.Gator1Data;
Gator2Data = ADiGator_logsumexp_ADiGatorHes.logsumexp_ADiGatorHes.Gator2Data;
% ADiGator Start Derivative Computations
cada2f1dx = x.dx;
cada2f1 = x.f(:);
cada2f2dx = exp(cada2f1(:)).*cada2f1dx;
cada2f2 = exp(cada2f1);
cada1f1dxdx = x.dx(:).*cada2f2dx;
cada1f1dx = cada2f2.*x.dx;
% Deriv 1 Line: cada1f1dx = exp(x.f(:)).*x.dx;
cada1f1 = exp(x.f);
% Deriv 1 Line: cada1f1 = exp(x.f);
cada1f2dxdx = cada1f1dxdx; cada1f2dx = cada1f1dx;
% Deriv 1 Line: cada1f2dx = cada1f1dx;
cada1f2 = sum(cada1f1);
% Deriv 1 Line: cada1f2 = sum(cada1f1);
cada2f1dx = -1./cada1f2.^2.*cada1f2dx;
cada2f1 = 1/cada1f2;
cada2tempdx = cada2f1dx(Gator2Data.Index1);
cada2tf1 = cada1f2dx(Gator2Data.Index2);
cada2td1 = cada2tf1(:).*cada2tempdx;
cada2td1(Gator2Data.Index3) = cada2td1(Gator2Data.Index3) + cada2f1.*cada1f2dxdx;
y.dxdx = cada2td1;
y.dx = cada2f1*cada1f2dx;
% Deriv 1 Line: y.dx = 1./cada1f2.*cada1f2dx;
y.f = log(cada1f2);
% Deriv 1 Line: y.f = log(cada1f2);
%User Line: y = log(sum(exp(x)));
y.dx_size = 32;
% Deriv 1 Line: y.dx_size = 32;
y.dx_location = Gator1Data.Index1;
% Deriv 1 Line: y.dx_location = Gator1Data.Index1;
y.dxdx_size = [y.dx_size,32];
y.dxdx_location = [y.dx_location(Gator2Data.Index4,:), Gator2Data.Index5];
end


function ADiGator_LoadData()
global ADiGator_logsumexp_ADiGatorHes
ADiGator_logsumexp_ADiGatorHes = load('logsumexp_ADiGatorHes.mat');
return
end