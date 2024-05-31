function [rf, bw] = design_rf_sinc (t, ta, nlobes, dt)

% DESIGN_RFGRAD_SS Design a SINC rf pulse for slice select excitation. temporal
% resolution is 4 us. Hanning filtering is applied.
%
% Usage: [rf, bw] = design_rf_sinc (t, ta, nlobes, dt)
%
% Returns
% -------
% rf: 1 x nt vector containing a sinc pulse in tesla.
% bw: bandwidth of the pulse in khz estimated by FT of pulse. 
%
% Expects
% -------
% t: pulse duration in ms (not including the ramp and rewinder of the gradient).
% ta: tip angle in deg when on resonance. defaults to 30 deg.
% 
% nlobes: num of side lobes of the sinc waveform. this number, along with the
% main lobe which accounts for 2 side lobes, determines the BWTP of the pulse,
% which is the total effective number of side lobes, 2*nlobes + 2. defaults to 3
% corresponding to a BWTP of 8.
%
% dt: dwell time in sec. defaults to 4e-6.
%
%
% See also: hanningFilter design_rf_slr design_rf_gauss design_gss save_rf_varian
%
%
% Copyright (C) 2008 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Fri Oct  3 15:43:10 2008
%

if nargin < 3
  nlobes = 3;
end
if nargin < 2
  ta = 30; % tip angle, deg
end
if nargin < 4
  dt = 4e-6; % s
end

gamma1 = 2.675e8; % for proton, in rad/s/tesla

t = t*1e-3; % s
nt = round(t/dt);
ta = deg2rad(ta);
nz = nlobes + 1;

% sinc pulse
mysinc = sinc(linspace(-nz, nz, nt));
rfmag = ta ./ (gamma1.*dt.* sum(mysinc)); % tesla
rf = rfmag.* hanningFilter(mysinc); % tesla

if nargout > 1
  nz = 2^20- length(rf); % num of zeros padded before FT.
  fftrf= abs(fftshift(fft([rf zeros(1,nz)])));

  bw = 1/dt/1000 * length(find((fftrf>=0.5*max(fftrf)))) / length(fftrf); % khz
end
