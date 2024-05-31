function [rfnew, rfphase] = modulate_rfphase (rfold, gss, r, dt, isSymmetric)

% MODULATE_RFPHASE Modulate the phases of slice selective rf pulses that have
% been designed for a particular slice in the way that the excited slice will be
% shifted along the slice selective gradient direction by a certain amount of
% distance.
%
% Usage: [rfnew, rfphase] = modulate_rfphase (rfold, gss, r, dt, isSymmetric)
%
% Returns
% -------
% rfnew: nchs x ntimepts matrix containing new multi-ch rf pulses.
% rfphase: rf phases appended to the old rf.
%
% Expects
% -------
% rfold: nchs x ntimepts matrix containing old multi-ch rf pulses.
% 
% gss: a scalar or a 1 x ntimepts vector for the slice selective gradient in Tesla/m
% 
% r: desired shift in cm
% 
% dt: dwell time in us. defaults to 4 us
%
% isSymmetric: a flag indicating if the phase accumulation is symmetric
% about the middle of the pulse. if true, the phase will go from negative
% to positive with 0 phase at the middle; if false, the phase starts from
% 0. defaults to false.
%
%
% See also: calc_freqoffset design_rf_gauss characterize_gradient
%
%
% Copyright (C) 2010 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Mon Oct  4 12:10:45 2010
%

if nargin< 4
  dt = 4;
end

if nargin<5
    isSymmetric=false;
end

dt = 1e-3*dt; % ms
gss = tpm2gpcm(gss);
gamma = 4.26; % for proton, in kHz/G

[nchs,nt] = size(rfold);
if isscalar(gss)
    gss(1:nt)=gss;
end

rfphase = calcRfPhase(gss(1:nt),r,gamma,dt); % exclude rewinder

if isSymmetric
    rfphase= rfphase-mean(rfphase);
end

rfPha = repmat(rfphase,[nchs 1]);

rfnew = rfold.* exp(1i.*rfPha);

disp('-> done.')


%% ====================
%%
%%  Local functions
%%
%% ====================

function rfphase = calcRfPhase (grad, r, gamma, dt)
foffset = -gamma.* r.* grad; % khz. NOTE: the minus sign may not be appropriate.

woffset = 2*pi*foffset; % rad/ms

rfphase = cumsum(woffset)* dt; % rad

