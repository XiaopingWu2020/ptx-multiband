%%% gtrapzDesigner.m --- A class for design of trapzoidal gradient waveforms,
%including triangle and blip gradients.
%% 
%% Filename: gtrapzDesigner.m
%% Description: 
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
%% Maintainer: 
%% Copyright (C) 2009 CMRR at UMN
%% Created: Fri Sep  4 22:59:05 2009 (CDT)
%% Version: 
%% Last-Updated: Tue Sep  8 20:21:19 2009 (CDT)
%%           By: Xiaoping Wu
%%     Update #: 33
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

classdef gtrapzDesigner < gradDesigner
    properties 
      
    end % properties
 
    methods 
        function obj = gtrapzDesigner (maxamp, maxsr, dt)
        % Constructor
        % Usage: obj = gtrapzDesigner (maxamp, maxsr, dt)
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

          obj = obj@gradDesigner(args{:});

        
        % Post init.
        % anything including accessing object
        
        end
        
        function grad = design (obj, gsurf)
        % Purpose: design the gradient object of a blip or trap gradient.
        % Usage: output = obj.design (gsurf)
        % gsurf: specified gradient surface.

          gsign = sign(gsurf);
          gsurf = abs(gsurf);
          if canUseTriangle(obj, gsurf)
            grad = createTriangle(obj, gsurf);
          else
            grad = createTrapzoid(obj, gsurf);
          end
        
          grad = gsign* grad;
        end
        
    end % methods
    
    methods (Access = protected)
      
      function flag = canUseTriangle (obj, gsurf)
        
        maxsr = obj.MaxSlewRate;
        maxamp = obj.MaxGradAmp;
        gsurfmax = maxamp.^2 / maxsr;
        
        flag = gsurfmax > gsurf;
      end
      
      function grad = createTriangle (obj, gsurf)
        amp = obj.calcAmp(gsurf);
        grad = [obj.createRampUp(amp) obj.createRampDown(amp)];
      end

      function grad = createTrapzoid (obj, gsurf)
        maxamp = obj.MaxGradAmp;
        maxsr = obj.MaxSlewRate;
        dt = obj.DwellTime;  
        
        time = (gsurf - maxamp.^2/ maxsr)./ maxamp;
        npts = round(time/dt);       
        
        gplat = maxamp* ones(1,npts);
        rampu = obj.createRampUp(maxamp);
        rampd = obj.createRampDown(maxamp);
        
        grad = [rampu gplat rampd];
          
      end
            
      function ramp = createRampUp (obj, amp)
        maxsr = obj.MaxSlewRate;
        dt = obj.DwellTime;
        
        time = amp/ maxsr;
        npts = floor(time/dt);
        
        ramp = (0:npts)* dt* maxsr;
      end
      
      function ramp = createRampDown (obj, amp)
        ramp = fliplr(createRampUp(obj, amp));
      end

      function amp = calcAmp (obj, gsurf)
      % calc amp of grad plateu
        amp = (gsurf* obj.MaxSlewRate).^0.5;
      end
      
    end % methods
end % classdef




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% gtrapzDesigner.m ends here
