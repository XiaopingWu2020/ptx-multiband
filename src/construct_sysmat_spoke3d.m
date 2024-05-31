function [sysmat, kp] = construct_sysmat_spoke3d (grad, b1maps, mask, fox, ...
                                                  b0map, dt, poffset)

% CONSTRUCT_SYSMAT_SPOKE3D Construct system matrix for pulse design using 3D
% spoke trajectory.
%
% Usage: [sysmat, kp] = construct_sysmat_spoke3d (grad, b1maps, mask, fox,
% b0map, dt, poffset)
%
% Returns
% -------
% sysmat: nspacepts-by-nspokes*nchs system matrix
%
% Expects
% -------
% grad: 3d gradient. 3-by-timepts matrix.
% b1maps: b1 maps
% mask: spatial mask
% fox: field of excitation in m.
% b0map: b0 map in tesla
% 
% dt: dwell time in s. defaults to 4e-6.
%
% poffset: [offsetx,offsety,offsetz] in mm specifying the offset of FOV with
% respect to grad isocenter. defaults to [0 0 0].
% 
%
% See also: construct_sysmat_spoke3d_phase calc_spoke_position
% assemble_grad_spoke3d
%
%
% Copyright (C) 2009 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Fri Sep 18 15:51:23 2009
%

if nargin < 5|| isempty(b0map)
  b0map = zeros(size(mask));
end
if nargin < 6
  dt = 4e-6;end % s
if nargin < 7
  poffset = [0, 0, 0];
end
  
[b1arr,posarr] = create_array(b1maps,mask,fox,1e-3*poffset);
[kp,slen,rlen] = calcK(grad,dt);

glen = length(grad(1,:));
gObj = gradPulse(ones(1,glen-rlen),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;

phasetrack = slen:slen:glen-rlen;
gamma = 2.675e8;
kernalmat = 1i.* gamma.* dt.*exp(1i* ( posarr*kp + b0map(mask)*kb0(phasetrack) ) );
                      
nchs = size(b1arr,2);
[nspa,nspo] = size(kernalmat);
sysmat = complex(zeros(nspa,nchs*nspo)); % this prealloc results in much faster construct
for idx = 1:nchs,
  iidx0 = (idx-1)*nspo + 1;
  sysmat(:,iidx0:(iidx0+nspo-1)) = diag(b1arr(:,idx)) * kernalmat;
end

disp('-> Construction of sysmat done.')

% -----------------------
function [kp,slen,rlen] = calcK (grad,dt)
% calc spoke positions
zz = find(grad(3,:)==0);
slen = zz(2); % lobe length
rlen = zz(end)-zz(end-2); % rewinder length

nspks = length(find(grad(3,:)==0))/ 2 - 1;

gradObj = gradPulse(grad,dt,'spoke3d','Tx');
ktraj = gradObj.KspaceTrajectory;

kp = zeros(2, nspks);
for idx= 1: nspks,
  kp(:,idx) = ktraj(1:2, round((idx-1)*slen+ slen/2));
end
