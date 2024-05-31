function g = design_toptgrad1D (g0, gn, kdes, gmax, srmax, dt)

% DESIGN_TOPTGRAD1D Design time optimal 1D gradient using convex optimization
% methods (cvx matlab toolbox). Here no kspace path is predefined. This is
% useful for designing rewinders, etc.
%
% Usage: g = design_toptgrad1D (g0, gn, kdes, gmax, srmax, dt)
%
% Returns
% -------
% g: output gradient in T/m
%
% Expects
% -------
% g0: initial grad amp in T/m
% gn: final grad amp in T/m
% kdes: desired k value in rad/m
% gmax: max grad amp in T/m
% srmax: max slew rate in T/m/s
% dt: dwell time in sec. defaults to 4 us.
%
%
% See also: minTimeGradient
%
%
% Copyright (C) 2010 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Thu May 13 15:30:29 2010
%

if nargin< 6
  dt = 4e-6;
end

gamma = 2.675e8; % gyromagnetic ratio, rad/sec/tesla

max_n = 500; % preset max number of time points.

gsum = kdes/gamma/dt;

%*********************************************************************
% use bisection algorithm to solve the problem
%*********************************************************************

n_bot = 1;
n_top = max_n;


while( n_top - n_bot > 1)
  % try to find a feasible design for given specs
  n = ceil( (n_top + n_bot)/2 );

% Build the matrix for slew rate constraints.
  As = -eye(n) + diag(diag(eye(n-1)),1);
  As = As(1:end-1,:)./dt;
  
  % formulate and solve the feasibility problem
cvx_begin quiet
cvx_precision high
  variable g(n);
     g(1) <= g0;
     g(1) >= g0;
     g(end) <= gn;
     g(end) >= gn;

     g <= gmax;
     g >= -gmax;
     As*g <= srmax;
     As*g >= -srmax;
     sum(g) >= gsum;
     sum(g) <= gsum;
cvx_end

  % bisection
  if strfind(cvx_status,'Solved') % feasible
    fprintf(1,'Problem is feasible for n = %d\n', n);
    n_top = n;
  else % not feasible
    fprintf(1,'Problem is not feasible for n = %d \n', n);
    n_bot = n;
  end
end

% optimal number of time points
n = n_top;
fprintf(1,'\nMinimum number of time points for given specs is %d.\n',n);


% Build the matrix for slew rate constraints.
  As = -eye(n) + diag(diag(eye(n-1)),1);
  As = As(1:end-1,:)./dt;

cvx_begin quiet
cvx_precision high
  variable g(n);
     g(1) <= g0;
     g(1) >= g0;
     g(end) <= gn;
     g(end) >= gn;

     g <= gmax;
     g >= -gmax;
     As*g <= srmax;
     As*g >= -srmax;
     sum(g) >= gsum;
     sum(g) <= gsum;
cvx_end

g = g(:).';

disp('-> Time optimal gradient designed..')
