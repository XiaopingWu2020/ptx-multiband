function kp = place_spoke_symmetric (fox, nspokes, theta0, nc)
% PLACE_SPOKE_SYMMETRIC Symmetrically place spokes on each circle. This means if an
% even number of spokes are desired, there will be no DC spoke.
%
% Usage: kp = place_spoke_symmetric (fox, nspokes, theta0, nc)
%
% Returns
% -------
% kp: spoke postions in kspace. 2-by-nspokes matrix.
%
% Expects
% -------
% fox: nominal field of excitation in cm.
% nspokes: total number of spokes.
% theta0: initial rotation in deg of the placement. defaults to 0.
%
% nc: number of circles. defaults to 1.
%
% See also: define_spoke_placement design_spoke3d
%
%
% Copyright (C) 2010 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Tue Aug  3 12:02:44 2010
%

if nargin< 3
  theta0=0;
end
if nargin< 4
  nc = 1;
end

if nspokes == 1
  kp= [0;0];
  disp('-> Done!')
  return;
end

r = 2e2*pi/fox;

theta0 = deg2rad(theta0);

isodd = true;
nspc = (nspokes - 1)/ nc;
if rem(nspokes,2) == 0
  isodd = false;
  nspc = nspokes/ nc;
end

dtheta = 2*pi/nspc;

nn=1:nspc;

rr=r* (nc:-1:1);

kk=[];
for idx= 1:nc,
  aa = theta0+ dtheta*nn;

  kkk = rr(idx).* exp(1i*aa);
  kk=[kk kkk];
  
  theta0= aa(end);
end

if isodd
  kp = [imag(kk) 0;real(kk) 0];else
  kp = [imag(kk);real(kk)];
end


disp('-> Done!')
