%%% gspkDesigner.m --- A class for design of slice selective gradients
%without rewinder, corresponding to a spoke in kspace
%% 
%% Filename: gspkDesigner.m
%% Description: 
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
%% Maintainer: 
%% Copyright (C) 2009 CMRR at UMN
%% Created: Fri Sep  4 10:57:52 2009 (CDT)
%% Version: 
%% Last-Updated: Thu Nov 18 15:06:34 2010 (CST)
%%           By: Xiaoping Wu
%%     Update #: 53
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

classdef gspkDesigner < gtrapzDesigner
    properties 
      %GAMMA = constants.mrConstants.GAMMA_H1;
      GAMMA = 2.675e8
    end % properties
   
 % -----------------------

    methods 
        function obj = gspkDesigner (maxamp, maxsr, dt)
        % Constructor
        % Usage: obj = gspkDesigner (maxamp, maxsr, dt)
        %
        % maxamp: max gradient amp available in tesla/m
        % maxsr: max slew rate available in tesla/m/s
        % dt: dwell time in sec
        
        
        % Pre init.
        % anything not using output (obj)
          if nargin < 3
            dt = 4e-6;
          end
          
          if nargin < 2
            maxsr = 166;
          end
          
          if nargin < 1
            maxamp = 50e-3;
          end
          
          args{1}= maxamp;
          args{2}= maxsr;
          args{3}= dt;

        % compvalue = classname.staticMethod();

        
        % Object init.
        % call super class, if applicable, before accessing object

          obj = obj@gtrapzDesigner(args{:});

        
        % Post init.
        % anything including accessing object

        % obj.classMethod();          
        end

        function disp(obj)
           disp@gtrapzDesigner(obj)
           fprintf(1,'Gamma: %f\n',obj.GAMMA);
        end
        
        function grad = design (obj, rf, thk)
        % Usage: grad = design (obj, rf, thk)
        % rf: the RF pulse, like a sinc.
        % thk: slice (slab) thickness in m.
 
          if ~isa(rf,'rfPulse')  
            rf = rfPulse(rf, obj.DwellTime, '');
          end
          
          amp = obj.calcAmp(rf.BandWidth, thk);
          
          plateu = amp* ones(1, rf.NumOfTimePoints);
          rampu = obj.createRampUp(amp);
          rampd = obj.createRampDown(amp);
          
          grad = [rampu, plateu, rampd];
          
        end

    end % methods
    
    methods (Access = protected)
      
      function amp = calcAmp (obj, bw, thk)
      % calc amp of grad plateu
        amp = 2*pi*bw ./ thk./ obj.GAMMA; % tesla/m
      end
    end % methods
end % classdef


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% gspkDesigner.m ends here
