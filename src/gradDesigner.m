%%% gradDesigner.m --- A class for design of gradients
%% 
%% Filename: gradDesigner.m
%% Description: 
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
%% Maintainer: 
%% Copyright (C) 2009 CMRR at UMN
%% Created: Fri Sep  4 10:48:55 2009 (CDT)
%% Version: 
%% Last-Updated: Fri Sep  4 15:01:05 2009 (CDT)
%%           By: Xiaoping Wu
%%     Update #: 9
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

classdef gradDesigner
    properties 
      MaxGradAmp = 50e-3 % tesla/meter
      MaxSlewRate = 166 % tesla/meter/sec
      DwellTime = 4e-6 % sec
    end % properties
 
    methods 
        function obj = gradDesigner(maxamp,maxsr,dt)
        %
         if nargin >0
          obj.MaxGradAmp = maxamp;
          obj.MaxSlewRate = maxsr;
          obj.DwellTime = dt;
         end
        end
        
        function disp(obj)
          fprintf(1,['Max grad amp:%d mT/m\nMax slew rate:%d T/m/s\nDwell ' ...
                     'time:%d us\n'],1e3*obj.MaxGradAmp,obj.MaxSlewRate, 1e6* ...
                  obj.DwellTime);
        end
        
    end % methods
end % classdef



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% gradDesigner.m ends here
