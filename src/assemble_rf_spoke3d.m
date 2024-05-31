function rf = assemble_rf_spoke3d (rfsub, wts, shifts)

% ASSEMBLE_RF_SPOKE3D Assemble RF for slice selective excitation with spoke
% trajectory. A shift of subpulse can be specified for certain channels to
% compensate for dB0 gradients thru plane.
%
% Usage: rf = assemble_rf_spoke3d (rfsub, wts)
%
% Returns
% -------
% rf: 
%
% Expects
% -------
% rfsub: vector for the subpulse
% 
% wts: weights, nspokes x nchs matrix
% 
% shifts: a vector containing number of elements in the sub rf that will be
% shifted for individual channels. positives shift the first sub pulse to the
% right, and negatives to the left. defaults to a zero vector indicating no
% shifts for all channels.
%
%
% See also: design_spoke3d construct_sysmat_spoke3d assemble_grad_spoke3d
%
%
% Copyright (C) 2009 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Sat Sep 19 21:22:20 2009
%

rfsub = rfsub(:).'; % make it a row vector.
[nspks,nchs] = size(wts);

if nargin < 3
  shifts = zeros(nchs,1);
end

slen = length(rfsub);
rf = complex(zeros(nchs, slen*nspks));

for iCh=1: nchs,
  mysign = 1;
  for jSpk=1:nspks,
    rfsub1= circshift(rfsub,[0,mysign*shifts(iCh)]);
    rf(iCh, (jSpk-1)*slen+1:jSpk*slen) = wts(jSpk,iCh).* rfsub1;
    mysign = -1*mysign;
  end
end
