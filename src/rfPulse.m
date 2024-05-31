%%% rfPulse.m --- A class for RF pulses
%% 
%% Filename: rfPulse.m
%% Description: 
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
%% Maintainer: 
%% Copyright (C) 2009 CMRR at UMN
%% Created: Thu Sep  3 00:50:54 2009 (CDT)
%% Version: 
%% Last-Updated: Fri Aug 27 10:17:36 2010 (CDT)
%%           By: Xiaoping Wu
%%     Update #: 36
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Commentary: 
%% 
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Change log:
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Code:

classdef rfPulse < pulse
    properties (Dependent = true, SetAccess = private)
      BandWidth % hz
    end % properties

% -----------------------

    methods 
        function obj = rfPulse (rf,dt,name)
        %
          if nargin < 3
            name = '';
          end
          
          if nargin < 2
            dt = 4e-6; % sec
          end
          
          if nargin < 1
            rf = [];
          end
          
          args{1}= rf;
          args{2}= dt;
          args{3}= name;
          
          obj = obj@pulse(args{:});
                    
        end
        % -----------------------
        
        function bw = get.BandWidth (obj)
          bw = obj.calcBW;
        end

        function plot (obj, varargin)
          time = obj.createTimeArray;

          amp = abs(obj.Pulse);
          pha = unwrap(angle(obj.Pulse));
          
          subplot(1,2,1), plot(1e3*time',amp',varargin{:});
          subplot(1,2,2), plot(1e3*time',pha',varargin{:});
          
        end
        
        function plot_amp (obj, varargin)
          time = obj.createTimeArray;
          
          amp = abs(obj.Pulse);
          
          plot(1e3*time',amp',varargin{:});
          xlabel('Time (ms)')
          ylabel('Amp (a.u.)')
          title('RF Waveform(s)')
          
        end
        
        function write(obj)
        % write rf shape to text file(s) for a specified system, e.g. Varian.
                    
          write_rfs_varian(obj.Pulse,obj.Name);
        end
                
    end % methods
    
    methods (Access = protected)
      function bw = calcBW (obj)
        % calc bandwidth in hz        
          npuls = obj.NumOfPulses;
          dt = obj.TimeStep;
          
          fftRf = abs(fftshift(fft(obj.Pulse, 2^20, 2)));          
  
          bw = zeros(npuls,1);
          for ip= 1: npuls,
            bw(ip) = 1/ dt* length(find((fftRf(ip,:)>=0.5*max(fftRf(ip,:))))) ...
                     / length(fftRf(ip,:)); % hz
          end
        end
    end % methods

end % classdef


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% rfPulse.m ends here
