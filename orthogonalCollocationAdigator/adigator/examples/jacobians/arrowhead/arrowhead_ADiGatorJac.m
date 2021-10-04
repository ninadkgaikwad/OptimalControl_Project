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

function y = arrowhead_ADiGatorJac(x)
global ADiGator_arrowhead_ADiGatorJac
if isempty(ADiGator_arrowhead_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_arrowhead_ADiGatorJac.arrowhead_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
%User Line: % Distributed under the GNU General Public License version 3.0
N.f = length(x.f);
%User Line: N = length(x);
y.f = zeros(N.f,1);
%User Line: y = zeros(N,1);
cada1f1dx = x.dx(1);
cada1f1 = x.f(1);
cada1f2dx = 2.*cada1f1.^(2-1).*cada1f1dx;
cada1f2 = cada1f1^2;
cada1f3dx = 2.*cada1f2dx;
cada1f3 = 2*cada1f2;
cada1f4dx = 2.*x.f(:).^(2-1).*x.dx;
cada1f4 = x.f.^2;
cada1f5dx = cada1f4dx;
cada1f5 = sum(cada1f4);
cada1td1 = zeros(100,1);
cada1td1(1) = cada1f3dx;
cada1td1 = cada1td1 + cada1f5dx;
cada1f6dx = cada1td1;
cada1f6 = cada1f3 + cada1f5;
y.dx = cada1f6dx;
y.f(1) = cada1f6;
%User Line: y(1) = 2*x(1)^2+sum(x.^2);
cada1f1dx = x.dx(1);
cada1f1 = x.f(1);
cada1f2dx = 2.*cada1f1.^(2-1).*cada1f1dx;
cada1f2 = cada1f1^2;
cada1f3 = 2:N.f;
cada1f4dx = x.dx(Gator1Data.Index1);
cada1f4 = x.f(cada1f3);
cada1f5dx = 2.*cada1f4(:).^(2-1).*cada1f4dx;
cada1f5 = cada1f4.^2;
cada1tempdx = cada1f2dx(Gator1Data.Index2);
cada1td1 = zeros(198,1);
cada1td1(Gator1Data.Index3) = cada1tempdx;
cada1td1(Gator1Data.Index4) = cada1td1(Gator1Data.Index4) + cada1f5dx;
cada1f6dx = cada1td1;
cada1f6 = cada1f2 + cada1f5;
cada1f7 = 2:N.f;
cada1td1 = zeros(298,1);
cada1td1(Gator1Data.Index5) = cada1f6dx;
cada1td1(Gator1Data.Index6) = y.dx(Gator1Data.Index7);
y.dx = cada1td1;
y.f(cada1f7) = cada1f6;
%User Line: y(2:N) = x(1)^2+x(2:N).^2;
y.dx_size = [100,100];
y.dx_location = Gator1Data.Index8;
end


function ADiGator_LoadData()
global ADiGator_arrowhead_ADiGatorJac
ADiGator_arrowhead_ADiGatorJac = load('arrowhead_ADiGatorJac.mat');
return
end