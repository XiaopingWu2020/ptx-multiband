function [grad, rfsub] = design_spoke3d (rf, thk, kp, useHann, maxsr, maxamp, ...
                                         dt)

% DESIGN_SPOKE3D Design 3d gradients corresponding to a 3D spoke trajectory
% provided the slice selective rf pulse, thickness and the spoke positions.
%
% Usage: [grad, rfsub] = design_spoke3d (rf, thk, kp, useHann, maxsr, maxamp, dt)
%
% Returns
% -------
% grad: 3d gradient
% rfsub: sub rf pulse
%
% Expects
% -------
% rf: subpulse used for slice (slab) selection
% thk: thickness in m
% kp: spoke postions in kspace. 2-by-nspokes matrix.
% 
% useHann: flag indicating how the gradient lobe and sub rf will be formed. true
% means using a hann window to filter both the step gradient and input rf. false
% means appending gradient ramps to the step gradient, and put zeros to both
% sides of input rf accordingly. defaults to true.
%
% maxsr: maximum gradient slew rate in T/m/s. defaults to 166.
% 
% maxamp: maximum gradient amplitude in mT/m. defaults to 50.
% 
% dt: dwell time in sec. defaults to 4e-6 s.
% 
% 
%
% See also: find_spoke_position assemble_rf_spoke3d construct_sysmat_spoke3d
% define_spoke_placement place_spoke_symmetric
%
% 
%
% Copyright (C) 2009 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Fri Sep 18 15:07:38 2009
%

% defaults
if nargin < 4
  useHann = true;
end

if nargin < 5
  maxsr = 166; % tesla/m/s
end
if nargin < 6
  maxamp = 50; 
end
maxamp = 1e-3* maxamp; %tesla/m

if nargin < 7
  dt = 4e-6; % s
end

%gamma = constants.mrConstants.GAMMA_H1;
gamma = 2.675e8;


% design spoke gradient and the subpulse with the same length
gd = gspkDesigner(maxamp,maxsr,dt);

if useHann % NOtE: this windowing results in distorted slice profiles
  slen = length(rf);
  hh = hann(slen);
  hh = hh(:).';
  gspoke = zeros(1,slen);
  gspoke(:) = max(abs(gd.design(rf,thk)));

  rfsub = rf.* hh;
  gspoke = gspoke.*hh;
else
  gspoke = gd.design(rf,thk);
  slen = length(gspoke);
  gsmax = max(gspoke);
  
  rfsub = complex(zeros(1,slen));
  rfsub(gspoke==gsmax) = rf;
end


hasDC = true;
if sum(abs(kp(:,end)))~= 0
  hasDC = false;
  kp = [kp [0;0]];
end

nspks = size(kp,2);


% assemble gradient
gz = [];
mysign = 1;

gd1 = gtrapzDesigner(maxamp,0.3*maxsr,dt);

if nspks > 1
  dkval = diff(kp,1,2);
  dkmax = max(abs(dkval(:)));

% design blip gradient  
  blipmax = gd1.design(dkmax/gamma);

  blen = length(blipmax);
  blipmax(slen) = 0;
  blips = [blipmax;blipmax];
  
  gxy = [];

  for ind= 1: nspks-1  
    gz = [gz mysign*gspoke];
    mysign= -1* mysign;
    
    sf = diag(dkval(:,ind)./ dkmax);
    
    gxy = [gxy sf*blips];  
  end
else
  gxy = [0;0];
end

gssr = gd1.design(0.5*dt*sum(-1*mysign*gspoke));

if hasDC
  gz = [gz, mysign*gspoke, gssr];else
  gz = [gz, -gssr];
end

gxy(:,size(gz,2)) = 0;

if nspks > 1
  gxy = circshift(gxy,[0, round(slen-blen/2)]);end


grad = [gxy; gz];

disp('-> Design of 3D spoke trajectory done.')
