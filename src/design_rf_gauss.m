function [rf, bw] = design_rf_gauss (t, ta, cutoff, dt)

% DESIGN_RF_GAUSS Design a Gaussian rf pulse for slice select excitation.
% temporal resolution is 4 us.
%
%
% Usage: [rf, bw] = design_rf_gauss (t, ta, cutoff, dt)
%
% Returns
% -------
% rf: 1 x nt vector containing a gaussian pulse in tesla.
% bw: bandwidth of the pulse in khz estimated by FT of pulse.
%
% Expects
% -------
% t: pulse duration in ms (not including the ramp and rewinder of the gradient).
% 
% ta: tip angle in deg when on resonance. defaults to 30 deg.
% 
% cutoff: a value showing where to cut the ideal gaussian curve for the rf
% waveform. defaults to 0.01, meaning that points that have a value less than 1
% percent of maximum are cutoff. Basically its good to keep this value less than
% 0.01.
% 
% dt: dwell time in sec. defaults to 4e-6.
%
%
% See also: design_rf_slr design_rf_sinc design_gss save_rf_varian
%
%
% Copyright (C) 2008 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Mon Oct 13 20:24:53 2008
%

if nargin < 3
  cutoff = 0.01; % 1 percent
end
if nargin < 2
  ta = 30; % tip angle, deg
end

if nargin < 4
  dt = 4e-6;
end

gamma1 = 2.675e8; % for proton, in rad/s/tesla

ta = deg2rad(ta);
t= 1e-3*t; % s
nt = round(t/dt);

% gaussian pulse
beta = abs(log(cutoff));
mygauss = getGauss(linspace(-1,1,nt),beta);

rfmag = ta ./ (gamma1.*dt.* sum(mygauss)); % tesla
rf = rfmag* mygauss;

if nargout > 1
  nz = 2^20- length(rf); % num of zeros padded before FT.
  fftrf= abs(fftshift(fft([rf zeros(1,nz)])));

  bw = 1/dt/1000 * length(find((fftrf>=0.5*max(fftrf)))) / length(fftrf); % khz
end


%% ====================
%%
%%  Local functions
%%
%% ====================

function y = getGauss(x, b)
% y = exp(-b*x^2), where -1<= x <=1

y = exp(-b*x.^2);
