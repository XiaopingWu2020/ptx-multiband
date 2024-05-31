%%% MBpTXpulseDesigner.m --- design multi-band ptx spoke RF pulses w/ total
%% rf power, global SAR or peak local SAR constraints.
%%
%%
%%
%% Filename: MBpTXadvPulseDesigner.m
%% Description:
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu>
%% Maintainer:
%% Copyright (C) 2012 CMRR at UMN
%% Created: Tue Aug 14 10:59:05 2012 (CDT)
%% Version:
%% Last-Updated: Mon Apr 15 16:57:45 2013 (CDT)
%%           By: Xiaoping Wu
%%     Update #: 1357
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

classdef MBpTXadvPulseDesigner < handle
    properties
        B1Map = []
        Mask = []
        Target = []
        FieldOfExcit = [] % m
        NumOfSpokes = 1
        SubRFType = 'gauss'
        Thickness = 0.005 % m
        SubRFDuration = 0.001 % s
        NominalFOX = 1e-2*(10:2:36) % m
        ConvergenceTolerance = 1e-5
        Lambda = [1e-8 1e-7 5e-7 1e-6 2e-6 5e-6 8e-6 1e-5 2e-5 5e-5]
        NominalFlipAngle = 10 % deg
        B0Map = []
        MaxGradSlewRate = 160 % T/m/s
        MaxGradAmplitude = 1e-3*50 % T/m
        ReadOutOffset = [0 0 0] % m
        MB = 2
        SARMatrix= []
        OptimalKspace = {[0;0]}
        RMSEd= 0.1
        GRasterTime= 10e-6
        OverSampleFactor= 5
        Tag= 'ptxMB'
        Orient= 'tranA2P'
        RegularizationType= 'tikhonov'
        USFactor= 8
        MaxNumOfRWIters= 100
        DesignStrategy= 'perb'
        PredictionIsNeeded= true
        OptimalX0IsNeeded= false
    end % properties
    
    properties (SetAccess = protected)
        DwellTime = 2e-6 % s
        NominalTheta = 0:20:160 % deg
        SubRF = []
        Gss= []
        Gr= []
        Gz= []
        Gxy= {}
        PhaseTrack= []
        TargetVector= {}
        BandPattern = {} % in slice index
        SliceLocation= []
        Weight= {}
        NumOfChannels= 0
        RF= {}
        iRFs= {}
        Gradient= {}
        IsOptimal = false
        IsInitialized= false
        IsCalculated= false
        OptimalGBlip= {}
        ResidualError = {}
        AMatrix= []
        RMatrix= {}
        SMatrix= {}
        mVectors= {}
        mVector= []
        mfullVector= []
        AfullMatrix= []
        X= {}
        NI= {}
        RHO= {}
        ETA= {}
        RMSE= {}
        MinMaxSAR= {}
        MaxSAR= {}
        MinSWA= {}
        W= {}
        G= {}
        NumOfCGIters= 20
        ToleranceMLS= 1e-4
        Mxy= {}
        MinRfEnergy= {}
    end % properties
    
    
    methods
        function obj = MBpTXadvPulseDesigner (b1map, mask, fox)
            % Constructor
            % Usage: obj = classname (b1map, mask, fox)
            % b1map:
            % mask:
            % fox: in cm
            
            % Pre init.
            % anything not using output (obj)
            
            if nargin == 0
                disp(['-> The object is constructed, but still requires necessary ' ...
                    'properties...'])
                return
            end
            
            if ndims(mask)~=3
                disp('-> 2D mask, ignored. Please give a 3D mask.')
                return
            end
            
            % compvalue = classname.staticMethod();
            
            
            % Object init.
            % call super class, if applicable, before accessing object
            % obj = obj@superClass(args{:});
            
            
            % Post init.
            % anything including accessing object
            % obj.classMethod();
            % obj.Property = compvalue;
            b1map(isnan(b1map))= 0;
            obj.B1Map = b1map;
            
            obj.Mask = mask;
            obj.FieldOfExcit = fox;
            
            obj.NumOfChannels= size(obj.B1Map,4);
            obj.Target= obj.Mask;
            
            
            disp('-> Object constructed. Please specify other needed properties...')
        end
        
        %% setors
        
        function obj= set.GRasterTime(obj, dt)
            
            if ~isequal(obj.GRasterTime,dt)
                obj.GRasterTime= dt;
                
                obj.setInitStatus(false);
                obj.setCalcStatus(false);
            end
            
        end
        
        function obj= set.OverSampleFactor(obj, osf)
            
            osf= round(osf);
            if ~isequal(obj.OverSampleFactor,osf)
                obj.OverSampleFactor= osf;
                
                obj.setInitStatus(false);
                obj.setCalcStatus(false);
            end
            
        end
        
        function obj = set.MB(obj, mb)
            
            if mb<2
                mb= 2;
                disp('-> MB must >=2, set MB=2.')
            end
            
            if obj.MB~= mb
                if ~obj.validateMB(mb);
                    disp('-> Inappropriate MB, ignored...')
                    return;
                end
                
                obj.MB= mb;
                
                obj.setInitStatus(false);
                obj.setOptimStatus(false);
                obj.setCalcStatus(false);
            end
            
        end
        
        
        function obj = set.NumOfSpokes (obj, nspokes)
            % Purpose:
            % Usage: output = obj.funcname (a,b,c)
            % a:
            % b:
            % c:
            
            if nspokes < 1
                nspokes = 1;
                disp('-> NumOfSpokes must be >= 1 and is set to 1.')
            end
            
            if ~isequal(obj.NumOfSpokes,nspokes)
                obj.NumOfSpokes = nspokes;
                obj.setOptimStatus(false);
                obj.setInitStatus(false);
                obj.setCalcStatus(false);
            end
            
        end
        
        
        function obj = set.Target (obj, targ)
            if ~isequal(obj.Target, targ)
                obj.Target = targ;
                obj.constructTargetVector();
                
                obj.setOptimStatus(false);
                obj.setInitStatus(false);
                obj.setCalcStatus(false);
            end
        end

       function obj = set.NominalFlipAngle (obj, fa)
            if ~isequal(obj.NominalFlipAngle, fa)
                obj.NominalFlipAngle = fa;
                obj.constructTargetVector();
                
                obj.setInitStatus(false);
                obj.setCalcStatus(false);
            end
        end

        function obj = set.NominalFOX (obj, myfox)
            % Purpose:
            % Usage: output = obj.funcname (a,b,c)
            % a:
            % b:
            % c:
            
            if ~isequal(obj.NominalFOX, myfox)
                obj.NominalFOX = myfox;
                obj.setOptimStatus(false);
                obj.setCalcStatus(false);
            end
            
        end
        
        function obj= set.Lambda(obj, lamb)
            %
            if ~isequal(lamb,obj.Lambda)
                obj.Lambda= lamb;
                obj.setCalcStatus(false);
            end
            
        end
        
        function obj= set.OptimalX0IsNeeded(obj, flag)
            
            if ~isequal(flag, obj.OptimalX0IsNeeded)
                obj.OptimalX0IsNeeded= flag;
                obj.setCalcStatus(false);
            end
            
        end
        
        
        function obj= set.Orient(obj, myorient)
            
            if ~(strcmpi(myorient,'tranA2P')||strcmpi(myorient,'tranR2L')|| ...
                    strcmpi(myorient,'sagA2P')||strcmpi(myorient,'corR2L'))
                disp(['-> Orient must be tranA2P, tranR2L, sagA2P ' ...
                    'or corR2L. set to current orient'])
                return;
            end
            
            obj.Orient= myorient;
            
        end
        
        
        function obj= set.DesignStrategy(obj, mystrategy)
            
            if ~(strcmpi(mystrategy,'perm')||strcmpi(mystrategy,'perb'))
                disp('-> DesignStrategy must be perm or perb. set to current type')
                mystrategy = obj.DesignStrategy;
            end
            
            oldtype= obj.DesignStrategy;
            obj.DesignStrategy = mystrategy;
            
            if ~strcmpi(oldtype,obj.DesignStrategy)
                obj.setCalcStatus(false);
            end
            
        end
        
        
        function obj= set.RegularizationType(obj, regtype)
            %
            if ~(strcmpi(regtype,'tikhonov')||strcmpi(regtype,'rfpower')|| ...
                    strcmpi(regtype,'gsar')||strcmpi(regtype,'lsar')|| ...
                    strcmpi(regtype,'gsarNoRFshape')||strcmpi(regtype,'lsarNoRFshape'))
                disp(['-> RegularizationType must be tikhonov, rfpower, gsar, gsarNoRFshape ' ...
                    'lsar or lsarNoRFshape. set to current type'])
                regtype = obj.RegularizationType;
            end
            
            oldtype= obj.RegularizationType;
            obj.RegularizationType = regtype;
            
            if ~strcmpi(oldtype,obj.RegularizationType)
                obj.setCalcStatus(false);
                obj.setInitStatus(false);
            end
            
        end
        
        function obj = set.SubRFType (obj, rftype)
            % Purpose:
            % Usage: output = obj.funcname (a,b,c)
            % a:
            % b:
            % c:
            
            if ~(strcmpi(rftype,'gauss')||strcmpi(rftype,'sinc10')||strcmpi(rftype,'sinc8')||...
                    strcmpi(rftype,'sinc6')||strcmpi(rftype,'sinc4'))
                disp('-> SubRFType must be sinc4, sinc6, sinc8, sinc10 or gauss. set to gauss')
                rftype = 'gauss';
            end
            
            oldtype= obj.SubRFType;
            obj.SubRFType = rftype;
            
            if ~strcmpi(oldtype,obj.SubRFType)
                obj.setInitStatus(false);
                obj.setCalcStatus(false);
            end
            
        end
        
        function obj= set.Thickness(obj, thk)
            
            if obj.Thickness~=thk
                obj.Thickness= thk;
                obj.setInitStatus(false);
                obj.setCalcStatus(false);
            end
            
        end
        
        
        function obj = set.B0Map (obj, b0map)
            % Purpose:
            % Usage: output = obj.funcname (a,b,c)
            % a:
            % b:
            % c:
            
            if isempty(b0map)
                b0map= zeros(size(obj.Mask));
            end
            
            if ~isequal(obj.B0Map, b0map)
                obj.B0Map = b0map;
                obj.setOptimStatus(false);
                obj.setCalcStatus(false);
            end
            
        end
        
        
        
        %% public methods
        
        function writeALL(obj)
        % a temporary writer that writes entire waveforms for only one MB excitation and
        % corresponding ref scans.
        if length(obj.BandPattern)>1
           disp('->Only ONE MB pulse can be written. ignored...')
           return
        end
        % write MB pulses whose pulse id has to be 0
                irf=1;
                [rf,grad]= obj.resamplePulse(irf);
                obj.writePulse(rf,grad,irf);
                
        % write SB pulses for ref scans whose pulse id's increment from 1 to MB factor.
                for ib= 1:obj.MB
                    [rf,grad]= obj.resamplePulseI(irf,ib);
                    write_pTXRFPulse(rf,grad,[obj.Tag,obj.Orient,'Sli',num2str(ib)], ...
                                     obj.Orient, ib, [], ...
                                     obj.OverSampleFactor);
                end
                      
        end
        
        function write(obj)
            % write MB rf and grad to machine readable text files
            for irf= 1:length(obj.BandPattern),
                [rf,grad]= obj.resamplePulse(irf);
                obj.writePulse(rf,grad,irf);
            end
            
        end
        
        function writeIP(obj)
            % write individual rf and grad to machine readable text files for ref scans
            for irf= 1:length(obj.BandPattern),
                for ib= 1:obj.MB
                    [rf,grad]= obj.resamplePulseI(irf,ib);
                    obj.writePulseI(rf,grad,obj.BandPattern{irf}(ib));
                end
            end
            
        end
        
        
        function design (obj)
            %
            if ~obj.validateMB(obj.MB)
                disp('-> Inappropriate MB, ignored...')
                return;
            end
            
            if ~obj.IsInitialized
                obj.init();end
            
            if isequal(obj.DesignStrategy,'perb')
                obj.designPERB();
            else % 'perm'
                obj.designPERM();
            end
            
            
            disp('-> Pulse calculations done...')
            
        end
        
        function rmag= calcRelativeVolt(obj)
            % calc volts used on the scanner for summed and individual reference
            % scans to retain constant flip angles. this is done by calculating max
            % amps of individual sub pulses wrt that of the summed pulse.
            
            rmag={};
            
            for irf= 1:length(obj.BandPattern),
                myrf={};
                myrf{1}.rfs= obj.RF{irf};
                for ib= 1:obj.MB
                    myrf{ib+1}.rfs= obj.iRFs{irf}{ib};
                end
                magr= compare_pTXRf(myrf,1);
                rmag{irf}=magr;
            end
        end
        
        
    end % methods
    
    %% protected methods
    methods (Access = protected)
        
        function designPERB(obj)
            
            if ~obj.IsOptimal
                disp('-> optimizing spokes k-space placement...')
                obj.findOptimalKspace();
            end
            
            disp('-> designing gradients...')
            obj.assembleGradient();
            
            disp('-> calculating RF weights...')
            obj.calculateWeights();
            
            disp('-> creating composite rf...')
            obj.composeRF();
            
        end
        
        function designPERM(obj)
            
            if ~obj.IsOptimal
                obj.findOptimalKspacePERM();
            end
            
            obj.assembleGradient();
            
            obj.calculateWeightsPERM();
            
            obj.composeRF();
            
        end
        
        %% - initializors
        
        function init(obj)
            obj.calcDwellTime();
            obj.designSubRF();
            obj.designGss();
            obj.calcPhasetrack();
            obj.createMBPattern();
            obj.calcSliLoc();
            obj.assembleGz();
            obj.createMvectors();
            obj.constructSMatrix();
            
            obj.setInitStatus(true);
        end
        
        
        function calcDwellTime(obj)
            obj.DwellTime= obj.GRasterTime/ obj.OverSampleFactor;
        end
        
        
        function createMvectors(obj)
            %
            for irf=1:length(obj.BandPattern)
                obj.mVectors{irf}= cell2mat(obj.TargetVector(obj.BandPattern{irf})).';
            end
            
        end
        
        function designSubRF (obj)
            % Purpose:
            rftype = obj.SubRFType;
            tp = 1e3* obj.SubRFDuration;
            dt = obj.DwellTime; % s
            switch rftype
                case 'gauss'
                    subrf = design_rf_gauss(tp, 1, 0.01, dt);
                case 'sinc10'
                    subrf = design_rf_sinc(tp, 1,4, dt);
                case 'sinc8'
                    subrf= design_rf_sinc(tp,1,3,dt);
                case 'sinc6'
                    subrf= design_rf_sinc(tp,1,2,dt);
                case 'sinc4'
                    subrf= design_rf_sinc(tp,1,1,dt);
                otherwise
                    subrf= [];
            end
            
            obj.SubRF= subrf;
            
        end
        
        function calcPhasetrack(obj)
            nspokes= obj.NumOfSpokes;
            lenGss= length(obj.Gss);
            len= nspokes* lenGss;
            phatrac= lenGss:lenGss:len;
            
            obj.PhaseTrack= phatrac;
            
        end
        
        function createMBPattern(obj)
            
            ns= size(obj.Mask,3);
            mb= obj.MB;
            
            nrfs= round(ns./mb);
            bandpat= cell(1,nrfs);
            bp= 1:nrfs:ns;
            for idx=1:nrfs
                bandpat{idx}= bp+ idx- 1;
            end
            
            obj.BandPattern= bandpat;
            
        end
        
        function calcSliLoc(obj)
            
            foxz= obj.FieldOfExcit(3);
            ns= size(obj.Mask,3);
            sliceshift= obj.ReadOutOffset(3);
            
            obj.SliceLocation= linspace(-foxz/2,foxz/2,ns) + sliceshift;
            
        end
        
        function designGss(obj)
            
            subrf= obj.SubRF;
            thk= obj.Thickness;
            maxsr= obj.MaxGradSlewRate;
            maxamp= obj.MaxGradAmplitude;
            dt= obj.DwellTime;
            gamma= 2.675e8;
            
            gd= gspkDesigner(maxamp,maxsr,dt);
            gss= gd.design(subrf,thk);
            gr = design_toptgrad1D(0,0,0.5*gamma*dt*sum(gss),maxamp,maxsr,dt);
            
            obj.Gss= gss;
            obj.Gr= gr;
            
        end
        
        function assembleGz(obj)
            %
            gspoke= obj.Gss;
            gz=[];
            for ind=1:obj.NumOfSpokes,
                gz = [gz gspoke];
                gspoke= -gspoke;
            end
            
            obj.Gz= [gz sign(sum(gspoke))*obj.Gr];
            
        end
        
        
        %%
        
        function findOptimalKspacePERM(obj)
            
            nrfs= length(obj.BandPattern);
            kp= cell(1,nrfs);
            rho= cell(1,nrfs);
            
            for irf=1:nrfs
                [ikp,irho]= obj.findOptimalKspaceMBperm(irf);
                kp{irf}= ikp;
                rho{irf}= irho;
            end
            
            obj.OptimalKspace= kp;
            obj.ResidualError= rho;
            
            obj.setOptimStatus(true);
            
        end
        
        function [kp0,myrho]= findOptimalKspaceMBperm (obj,irf)
            % Purpose:
            nspokes= obj.NumOfSpokes;
            if nspokes==1
                kp0= [0;0];
                myrho= [];
                return;
            end
            
            obj.mVector= obj.mVectors{irf};
            obj.mfullVector= obj.mVector;
            
            myfox = obj.NominalFOX;
            mytheta = obj.NominalTheta;
            myrho = zeros(length(myfox),length(mytheta));
            mykp = cell(length(myfox),length(mytheta));
            for ifox= 1:length(myfox),
                myfox1 = myfox(ifox);
                for jtheta= 1:length(mytheta)
                    kp= place_spoke_symmetric(1e2*myfox1,nspokes, ...
                        mytheta(jtheta),1);
                    
                    obj.constructSystemMatrixMBperm(kp,irf);
                    
                    obj.AfullMatrix= obj.AMatrix;
                    
                    [~,~,rho]= obj.solveCGMLS(1);
                    
                    myrho(ifox,jtheta) = rho;
                    mykp{ifox,jtheta} = kp;
                end
            end
            
            [ifox,jtheta] = find(myrho==min(myrho(:)));
            kp0= mykp{ifox(1),jtheta(1)};
            
        end
        
        %
        function findOptimalKspace(obj)
            
            nrfs= length(obj.BandPattern);
            kp= cell(1,nrfs);
            rho= cell(1,nrfs);
            
            for irf=1:nrfs
                [ikp,irho]= obj.findOptimalKspaceMB(irf);
                kp{irf}= ikp;
                rho{irf}= irho;
            end
            
            obj.OptimalKspace= kp;
            obj.ResidualError= rho;
            
            obj.setOptimStatus(true);
            
        end
        
        function [kp0,myrho]= findOptimalKspaceMB (obj,irf)
            % Purpose:
            nspokes= obj.NumOfSpokes;
            if nspokes==1
                kp0= [0;0];
                myrho= [];
                return;
            end
            
            obj.mVector= obj.mVectors{irf};
            obj.mfullVector= obj.mVector;
            
            myfox = obj.NominalFOX;
            mytheta = obj.NominalTheta;
            myrho = zeros(length(myfox),length(mytheta));
            mykp = cell(length(myfox),length(mytheta));
            for ifox= 1:length(myfox),
                myfox1 = myfox(ifox);
                for jtheta= 1:length(mytheta)
                    kp= place_spoke_symmetric(1e2*myfox1,nspokes, ...
                        mytheta(jtheta),1);
                    
                    obj.constructSystemMatrixMB(kp,irf);
                    
                    obj.AfullMatrix= obj.AMatrix;
                    
                    [~,~,rho]= obj.solveCGMLS(1);
                    
                    myrho(ifox,jtheta) = rho;
                    mykp{ifox,jtheta} = kp;
                end
            end
            
            [ifox,jtheta] = find(myrho==min(myrho(:)));
            kp0= mykp{ifox(1),jtheta(1)};
            
        end
        
        
        function reshapeWeight(obj)
            
            if isequal(obj.DesignStrategy,'perb')
                for idx=1:length(obj.Weight),
                    wts= obj.Weight{idx};
                    obj.Weight{idx}= reshape(wts(:),[],obj.MB);
                end
            else % 'perm'
                for idx=1:length(obj.Weight),
                    wts= obj.Weight{idx};
                    obj.Weight{idx}= repmat(wts(:),[1 obj.MB]);
                end
            end
            
        end
        
        
        function weightRF(obj)
            
            gss= obj.Gss;
            rfl= complex(zeros(size(gss)));
            rfl(gss==max(abs(gss)))= obj.SubRF;
            
            obj.iRFs= {};
            for irf=1:length(obj.Weight),
                if isempty(obj.Weight{irf})
                    obj.iRFs{irf}= [];
                    break;
                end
                
                for ib=1:obj.MB,
                    iwt= reshape(obj.Weight{irf}(:,ib),[],obj.NumOfChannels);
                    obj.iRFs{irf}{ib}= assemble_rf_spoke3d(rfl,iwt);
                end
            end
            
        end
        
        
        function predictMxy (obj)
            
            dt= obj.DwellTime;
            poffset= obj.ReadOutOffset;
            fox= obj.FieldOfExcit(1:2);
            
            obj.Mxy= {};
            for irf=1:length(obj.BandPattern),
                if isempty(obj.iRFs{irf})
                    break;
                end
                
                iRF= obj.iRFs{irf};
                gxy= obj.Gxy{irf};
                
                soi= obj.BandPattern{irf};
                for ib=1:obj.MB,
                    b1map= obj.B1Map(:,:,soi(ib),:);
                    mask= obj.Mask(:,:,soi(ib));
                    b0map= obj.B0Map(:,:,soi(ib));
                    rf= iRF{ib};
                    
                    mxypat= run_bloch_sim(rf,gxy,b1map,mask,fox,b0map,0,[],dt, ...
                        1e3*poffset);
                    
                    obj.Mxy{irf}{ib}= mxypat;
                end
            end
            
        end
        
        
        function phaseRF(obj)
            
            for irf=1:length(obj.iRFs),
                if isempty(obj.iRFs{irf})
                    break;
                end
                
                for ib=1:obj.MB,
                    obj.iRFs{irf}{ib}=modulate_rfphase(obj.iRFs{irf}{ib},...
                        obj.Gz,...
                        1e2*obj.SliceLocation(obj.BandPattern{irf}(ib)),...
                        1e6*obj.DwellTime);
                end
            end
            
        end
        
        
        function assembleRF(obj)
            
            obj.RF= {};
            for irf=1:length(obj.iRFs),
                if isempty(obj.iRFs{irf})
                    break;
                end
                
                rf=0;
                for ib=1:obj.MB,
                    rf= rf+obj.iRFs{irf}{ib};
                end
                
                obj.RF{irf}= rf;
            end
            
        end
        
        
        function B= buildCMatrix(obj,irf,ib)
            %
            
            rf= obj.iRFs{irf}{ib};
            rf= rf(:);
            rf= reshape(rf,[],obj.NumOfSpokes);
            rf= rf(obj.Gss==max(obj.Gss),:);
            rf= rf(1:obj.USFactor:end,:);
            rfb=[];
            for idx=1:obj.NumOfSpokes
                rfb= blkdiag(rfb,rf(:,idx));
            end
            
            B= [];
            for idx=1:obj.NumOfChannels,
                B= blkdiag(B,rfb);
            end
            
        end
        
        
        function buildSCMatrix(obj,irf)
            % for regs with rf shape incorporated.
            
            iSC= [];
            switch obj.RegularizationType
                case 'rfpower'
                    for ib=1:obj.MB,
                        iSC= [iSC obj.buildCMatrix(irf,ib)];
                    end
                    obj.RMatrix{1}= sparse(iSC/max(abs(iSC(:))));
                    
                case 'gsar'
                    for ib=1:obj.MB,
                        iSC= [iSC obj.SMatrix{1}*obj.buildCMatrix(irf,ib)];
                    end
                    obj.RMatrix{1}= sparse(iSC/max(abs(iSC(:))));
                    
                otherwise % 'lsar'
                    for idx=1:length(obj.SMatrix)
                        iSC=[];
                        for ib=1:obj.MB,
                            iSC= [iSC obj.SMatrix{idx}*obj.buildCMatrix(irf,ib)];
                        end
                        obj.RMatrix{idx}= sparse(iSC/max(abs(iSC(:))));
                    end
            end
            
        end
        
        
        function constructSMatrix(obj)
            
            regtype= obj.RegularizationType;
            if isequal(regtype,'tikhonov')||isequal(regtype,'rfpower')
                return;
            end
            
            if isempty(obj.SARMatrix)
                error('-> please specify obj.SARMatrix for sar constraint...')
            end
            
            if isequal(obj.DesignStrategy,'perb')&&(isequal(regtype,'gsar')||isequal(regtype,'lsar'))
                
                obj.SMatrix= construct_Svop(obj.SARMatrix,obj.NumOfSpokes* ...
                    length(obj.SubRF(1:obj.USFactor:end)));
                
            else
                obj.SMatrix= construct_Svop(obj.SARMatrix,obj.NumOfSpokes);
            end
            
        end
        
        
        function generatePhasedRf(obj,irf)
            
            rfl= complex(zeros(size(obj.Gss)));
            rfl(obj.Gss==max(obj.Gss))= obj.SubRF;
            myrf= repmat(rfl,[1 obj.NumOfSpokes]);
            
            bp= obj.BandPattern{irf};
            mb= length(bp);
            for ib=1:mb
                obj.iRFs{irf}{ib}=modulate_rfphase(myrf,obj.Gz, ...
                    1e2*obj.SliceLocation(bp(ib)), ...
                    1e6*obj.DwellTime);
            end
            
        end
        
        
        function assembleGxy(obj)
            %
            npts= obj.NumOfSpokes;
            if npts==1
                gxy=[0;0];
                gxy(:,length(obj.Gz))= 0;
                
                for irf=1:length(obj.BandPattern)
                    obj.Gxy{irf}= gxy;
                end
                
                return;
            end
            
            if ~mod(npts,2) % true for even number
                npts= npts+1;
            end
            
            mylen= length(obj.Gss);
            
            gblip= obj.OptimalGBlip;
            
            for irf=1:length(gblip)
                
                grad= gblip{irf};
                gxy= [];
                lof= 0;
                for idx=1:npts-1
                    ig= grad{idx};
                    glen= size(ig,2);
                    ig(:,mylen+floor(glen/2))= 0;
                    
                    gxy= [gxy, circshift(ig,[0 mylen-ceil(glen/2)-lof])];
                    
                    lof= floor(glen/2);
                end
                
                gxy(:,length(obj.Gz))= 0;
                
                obj.Gxy{irf}= gxy;
                
            end
            
        end
        
        
        function assembleGradient(obj)
            
            obj.designGBlip();
            obj.assembleGxy();
            
            for irf= 1:length(obj.BandPattern)
                obj.Gradient{irf}= [obj.Gxy{irf};obj.Gz];
            end
            
        end
        
        function composeRF(obj)
            
            obj.reshapeWeight();
            obj.weightRF();
            
            if obj.PredictionIsNeeded
                obj.predictMxy();
            end
            
            obj.phaseRF();
            obj.assembleRF();
            
        end
        
        
        function constructTargetVector (obj)
            % Purpose:
            mask = obj.Mask;
            targ = obj.Target;
            fa0= obj.NominalFlipAngle;
            
            ns= size(mask,3);
            targv= cell(1,ns);
            for ind=1:ns
                itargv= construct_target_vector(targ(:,:,ind),mask(:,:,ind),fa0);
                targv{ind}= itargv(:).';
            end
            
            obj.TargetVector= targv;
            
        end
        
        
        function createM(obj, irf)
            
            m= obj.mVectors{irf};
            mfull= m;
            mfull(length(m)+ size(obj.RMatrix{1},1)*length(obj.RMatrix))= 0;
            
            obj.mfullVector= sparse(mfull);
            obj.mVector= m;
            
        end
        
        function buildRMatrixSimple (obj)
            % regularization matrix construction for perm
            obj.RMatrix= {};
            regtype= obj.RegularizationType;
            if isequal(regtype,'tikhonov')||isequal(regtype, 'rfpower')
                [~,nn]= size(obj.AMatrix);
                obj.RMatrix{1}= eye(nn);
                
            else% 'lsar' or 'gsar'
                
                for idx=1:length(obj.SMatrix)
                    obj.RMatrix{idx}= sparse(obj.SMatrix{idx}/max(abs(obj.SMatrix{idx}(:))));
                end
                
            end
            
        end
        
        
        function buildRMatrix (obj,irf)
            % regularization matrix construction for perb
            obj.RMatrix= {};
            regtype= obj.RegularizationType;
            
            if isequal(regtype,'tikhonov')
                
                [~,nn]= size(obj.AMatrix);
                obj.RMatrix{1}= sparse(eye(nn));
                
            elseif isequal(regtype,'gsarNoRFshape')
                iSC=[];
                for ib=1:obj.MB,
                    iSC= blkdiag(iSC, obj.SMatrix{1});
                end
                obj.RMatrix{1}= sparse(iSC/max(abs(iSC(:))));
                
            elseif isequal(regtype,'lsarNoRFshape')
                for idx=1:length(obj.SMatrix)
                    iSC=[];
                    for ib=1:obj.MB,
                        iSC= blkdiag(iSC, obj.SMatrix{idx});
                    end
                    obj.RMatrix{idx}= sparse(iSC/max(abs(iSC(:))));
                end
            else
                
                obj.generatePhasedRf(irf);
                obj.buildSCMatrix(irf);
                
            end
            
        end
        
        
        function constructSystemMatrixMBfull(obj,idx)
            % construct Afull
            obj.AfullMatrix= [obj.AMatrix; sqrt(obj.Lambda(idx)).*obj.RMatrix{1}];
        end
        
        
        function constructSystemMatrixMBperm (obj,kp,irf)
            %
            bp= obj.BandPattern{irf};
            mb = length(bp);
            
            sysmat= [];
            for ib= 1:mb,
                sysmat= [sysmat; obj.constructSystemMatrix(kp,bp(ib))];
            end
            
            obj.AMatrix= sysmat;
            
        end
        
        function constructSystemMatrixMB (obj,kp,irf)
            %
            bp= obj.BandPattern{irf};
            mb = length(bp);
            
            sysmat= [];
            for ib= 1:mb,
                sysmat= blkdiag(sysmat, obj.constructSystemMatrix(kp,bp(ib)));
            end
            
            obj.AMatrix= sysmat;
            
        end
        
        function sysmat= constructSystemMatrix (obj,kp,soi)
            
            fox= obj.FieldOfExcit;
            phasetrack= obj.PhaseTrack;
            dt= obj.DwellTime;
            poffset= obj.ReadOutOffset;
            
            b1maps= obj.B1Map(:,:,soi,:);
            mask= obj.Mask(:,:,soi);
            b0map= obj.B0Map(:,:,soi);
            
            [b1arr,posarr] = create_array(b1maps,mask,fox,poffset);
            
            gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
            kb0 = gObj.KspaceTrajectory;
            
            gamma= 2.675e8;
            kernalmat = 1i*gamma*sum(obj.SubRF)*dt*...
                exp(1i* ( posarr*kp + b0map(mask)*kb0(phasetrack) ) );
            
            nchs = obj.NumOfChannels;
            sysmat= cell(1,nchs);
            parfor idx = 1:nchs,
                sysmat{idx} = diag(b1arr(:,idx)) * kernalmat;
            end
            sysmat= cell2mat(sysmat);
            
        end
        
        
        
        function stat= validateMB(obj,mb)
            
            ns= size(obj.Mask,3);
            stat= ~mod(ns,mb);
            
        end
        
        
        function setOptimStatus(obj,stat)
            obj.IsOptimal= stat;
        end
        
        function setInitStatus(obj,stat)
            obj.IsInitialized= stat;
        end
        
        function setCalcStatus(obj, stat)
            obj.IsCalculated= stat;
        end
        
        
        function designGBlip(obj)
            
            if obj.NumOfSpokes==1
                obj.OptimalGBlip= {[0;0]};
                return;
            end
            
            kp= obj.OptimalKspace;
            gblip= cell(size(kp));
            for idx=1:length(gblip)
                gblip{idx}= obj.designGBlipMB(kp{idx});
            end
            
            obj.OptimalGBlip= gblip;
            
        end
        
        function grad= designGBlipMB(obj,kp)
            
            npts = size(kp,2);
            grad={};
            
            maxamp = obj.MaxGradAmplitude;
            maxsr = obj.MaxGradSlewRate;
            dt = obj.DwellTime;
            
            if ~mod(npts,2) % even number
                kp = [kp, zeros(2,1)];
            end
            
            dkval = diff(kp,1,2);
            for ind= 1: length(dkval(1,:))
                dkmax = max(abs(dkval(:,ind)));
                blipmax = design_toptgrad1D(0,0,dkmax,maxamp,maxsr,dt);
                blips = [blipmax;blipmax];
                sf = diag(dkval(:,ind)./ dkmax);
                
                grad{ind} = sf*blips;
            end
            
        end
        %
        
        function calculateWeights(obj)
            
            if ( isequal(obj.RegularizationType,'lsar')||...
                    isequal(obj.RegularizationType,'lsarNoRFshape') )
                obj.calcWeightsLocalSARreg();
            else
                obj.calcWeights();
            end
            
        end
        
        
        function calculateWeightsPERM(obj)
            
            if ~isequal(obj.RegularizationType,'lsar')
                obj.calcWeightsPERM();
            else % 'lsar'
                obj.calcWeightsLocalSARregPERM();
            end
            
        end
        
        
        function calcWeightsLocalSARregPERM(obj)
            %
            for irf=1:length(obj.BandPattern),
                obj.constructSystemMatrixMBperm(obj.OptimalKspace{irf},irf);
                
                obj.buildRMatrixSimple();
                
                obj.createM(irf);
                
                obj.solveLSAR(irf);
            end
        end
        
        function calcWeightsLocalSARreg(obj)
            %
            for irf=1:length(obj.BandPattern),
                obj.constructSystemMatrixMB(obj.OptimalKspace{irf},irf);
                
                obj.buildRMatrix(irf);
                
                obj.createM(irf);
                
                obj.solveLSAR(irf);
            end
        end
        
        
        function calcWeightsPERM (obj)
            % Purpose:
            if ~obj.IsCalculated
                
                obj.reset();
                
                for irf=1:length(obj.BandPattern),
                    obj.constructSystemMatrixMBperm(obj.OptimalKspace{irf},irf);
                    
                    obj.buildRMatrixSimple();
                    
                    obj.createM(irf);
                    
                    obj.solve(irf);
                end
                
                obj.setCalcStatus(true);
                
            end
            
            for irf=1:length(obj.BandPattern)
                obj.pickOptimalSolution(irf);
            end
            
        end
        
        function calcWeights (obj)
            % Purpose:
            if ~obj.IsCalculated
                
                obj.reset();
                
                for irf=1:length(obj.BandPattern),
                    obj.constructSystemMatrixMB(obj.OptimalKspace{irf},irf);
                    
                    obj.buildRMatrix(irf);

                    obj.createM(irf);
                        
                    disp('-> solving the inverse problem...')
                    obj.solve(irf);
                end
                
                obj.setCalcStatus(true);
                
            end
            
            for irf=1:length(obj.BandPattern)
                obj.pickOptimalSolution(irf);
            end
            
        end
        
        
        function solveLSAR(obj, irf)
            
            [xopt,minSARmax,xx0,SARmax,SWAmin,myRMSE,myW,myG,myRHO,myETA]=obj.solveIRWMLS();
            
            obj.Weight{irf}= xopt;
            obj.MinMaxSAR{irf}= minSARmax;
            obj.X{irf}= xx0;
            obj.MaxSAR{irf}= SARmax;
            obj.MinSWA{irf}= SWAmin;
            obj.RMSE{irf}= myRMSE;
            obj.W{irf}= myW;
            obj.G{irf}= myG;
            obj.RHO{irf}= myRHO;
            obj.ETA{irf}= myETA;
            
        end
        
        
        function [xopt,minSARmax,xx0,SARmax,SWAmin,myRMSE,W,G,RHO,ETA]=solveIRWMLS(obj)
            
            S= obj.RMatrix;
            A= obj.AMatrix;
            lamb= obj.Lambda;
            rmse0= obj.RMSEd;
            
            nitsMax= obj.MaxNumOfRWIters;
            dw= 0.1;
            rswa0= 1e-4;
            
            % init weights
            
            ns= length(S);
            w= 1/ns*ones(ns,1);
            
            [~,nn]= size(A);
            [mm1,~]= size(S{1});
            
            % calc for several lamb
            
            nl= length(lamb);
            RHO= zeros(nl,1);
            ETA=RHO;
            XX= complex(zeros(nn,nl));
            
            % iterations begin here
            myRMSE= zeros(nl,nitsMax);
            myrmse= zeros(nl,1);
            nits= 0;
            g= zeros(ns,1);
            G= zeros(ns,nitsMax);
            W= G;
            SWAmin= zeros(1,nitsMax);
            xx0= complex(zeros(nn,nitsMax));
            SARmax= zeros(1,nitsMax);
            
            rswa= inf;
            swamin_pre= 1e-22;
            
            As= complex(zeros(mm1*ns,nn));
            while (rswa>= rswa0) && (nits< nitsMax)
                nits= nits+ 1;
                
                for idx=1:ns,
                    As(((idx-1)*mm1+1):idx*mm1,:)= sqrt(w(idx)).* S{idx};
                end
                
                for jdx=1:nl
                    
                    % construct Afull
                    obj.AfullMatrix= [A; sqrt(lamb(jdx))*As];
                    
                    [x,~,rho,eta,rmse]= obj.solveCGMLS(1);
                    XX(:,jdx)= x;
                    RHO(jdx)= rho;
                    ETA(jdx)= eta;
                    myrmse(jdx)= rmse;
                end
                
                % pick solutions that satisfy performance
                X0= XX(:,myrmse<=rmse0);
                if isempty(X0)
                    disp('-> no solution found to satisfy performance. break...')
                    break;
                end
                
                % calc sum_j(w_j* |S_j*x|^2) for those picked solutions
                Swa= w(1).*S{1}'*S{1};
                for idx=2:ns
                    Swa= Swa+ w(idx).*S{idx}'*S{idx};
                end
                
                [~,nx0]= size(X0);
                swa= zeros(nx0,1);
                for idx=1:nx0
                    ix= X0(:,idx);
                    swa(idx)= real(ix'*Swa*ix);
                end
                
                % find solution that minimizes sum_j(w_j* |S_j*x|^2)
                swamin= min(swa);
                x0=X0(:,swa==swamin);
                
                % update weightings using gradient decent
                for idx=1:ns
                    g(idx)= real(x0'*S{idx}'*S{idx}*x0);
                end
                
                SARmax(nits)= max(g);
                
                g=g/norm(g);
                w= w+ dw*norm(w)*g;
                w= w/sum(w);
                
                G(:,nits)= g;
                W(:,nits)= w;
                SWAmin(nits)= swamin;
                xx0(:,nits)= x0;
                myRMSE(:,nits)= myrmse;
                
                rswa= abs(swamin- swamin_pre)./swamin_pre;
                swamin_pre= swamin;
                
                fprintf('-> %d iteration(s) done...\n', nits)
            end
            % iterations end
            
            G= G(:,1:nits);
            W= W(:,1:nits);
            SWAmin= SWAmin(1:nits);
            SARmax= SARmax(1:nits);
            xx0= xx0(:,1:nits);
            myRMSE= myRMSE(:,1:nits);
            
            minSARmax= min(SARmax);
            xopt= xx0(:,SARmax==minSARmax);
            
        end
        
        
        function solve (obj,irf)
            %
            for jdx=1:length(obj.Lambda)
                obj.constructSystemMatrixMBfull(jdx);
                if obj.NumOfSpokes>1|| ~obj.OptimalX0IsNeeded
                    [x,ni,rho,eta,rmse]= obj.solveCGMLS(1);
                    
                else % single spoke and need optimal x0
                    disp('-> trying different x0 values...')
                    [x,ni,rho,eta,rmse]= obj.solveCGMLS(100);
                end
                
                obj.X{irf}(:,jdx)= x;
                obj.NI{irf}(jdx)= ni;
                obj.RHO{irf}(jdx)= rho;
                obj.ETA{irf}(jdx)= eta;
                obj.RMSE{irf}(jdx)= rmse;
            end
            
        end
        
        function [x, ni, rho, eta, rmse] = solveCGMLS (obj,ntrials)
            % solve magnitude least squares tikhonov regularization problem with CG.
            
            A= obj.AfullMatrix;
            b= obj.mfullVector;
            
            tol= obj.ToleranceMLS;
            ncgs= obj.NumOfCGIters;
            
            n= size(A,2);
            if ntrials==1
                x0= zeros(n,ntrials);
            else
                x0= rand(n,ntrials)+ 1i*rand(n,ntrials);
            end
            
            xs= x0;
            xs(:)=0;
            maxmins= zeros(1,ntrials);
            rhos=maxmins;
            etas=maxmins;
            nis= maxmins;
            
            disp('-> solving cgmls...')
            parfor itrial= 1:ntrials
                [ix,ini,irho,ieta,im]= MBpTXadvPulseDesigner.solve_cgmls(A,b,tol, ...
                    x0(:,itrial), ...
                    ncgs);
                
                if any(isnan(im))
                    maxmins(itrial)= inf;
                else
                    maxmins(itrial)= (max(abs(im))-min(abs(im)))./ mean(abs(im));
                end
                
                xs(:,itrial)= ix;
                rhos(itrial)= irho;
                etas(itrial)= ieta;
                nis(itrial)= ini;
            end
            
            % find the solution that gives least max min diff.
            myidx= find(maxmins== min(maxmins));
            x= xs(:,myidx(1));
            ni= nis(myidx(1));
            rho= rhos(myidx(1));
            eta= etas(myidx(1));
            
            %
            if obj.NominalFlipAngle<=30
                rmse= sqrt(sum( ( abs(obj.AMatrix*x)-abs(obj.mVector) ).^2) ./ ...
                    length(obj.mVector));
            else
                rmse= obj.simulateRMSE(x);
            end
            
        end
        % -----------------------
        
        
        function rmse= simulateRMSE(obj, wts)
            % currently only works for the case of single MB pulse design.
            irf=1;
            
            % reshape wts
            if isequal(obj.DesignStrategy,'perb')
                wts= reshape(wts(:),[],obj.MB);
            else
                wts= repmat(wts(:),[1 obj.MB]);
            end
            
            % weight rf
            gss= obj.Gss;
            rfl= complex(zeros(size(gss)));
            rfl(gss==max(abs(gss)))= obj.SubRF;
            iRF={};
            for ib=1:obj.MB,
                iwt= reshape(wts(:,ib),[],obj.NumOfChannels);
                iRF{ib}= assemble_rf_spoke3d(rfl,iwt);
            end
            
            % run bloch sim
            dt= obj.DwellTime;
            poffset= obj.ReadOutOffset;
            fox= obj.FieldOfExcit(1:2);
            
            
            gxy= obj.Gxy{irf};
            
            soi= obj.BandPattern{irf};
            mxy=[];
            for ib=1:obj.MB,
                b1map= obj.B1Map(:,:,soi(ib),:);
                mask= obj.Mask(:,:,soi(ib));
                b0map= obj.B0Map(:,:,soi(ib));
                rf= iRF{ib};
                
                mxypat= run_bloch_sim(rf,gxy,b1map,mask,fox,b0map,0,[],dt, ...
                    1e3*poffset);
                
                mxy= [mxy;mxypat(mask)];
            end
            
            rmse= sqrt(sum((abs(mxy)-abs(obj.mVector)).^2) ./ ...
                length(obj.mVector));
            
        end
        
        
        function pickOptimalSolution(obj, irf)
            % pick solutions that satisfy performance
            X0= obj.X{irf}(:,obj.RMSE{irf}<= obj.RMSEd);
            if isempty(X0)
                disp('-> no solution found to satisfy performance. break...')
                obj.Weight{irf}=[];
                return;
            end
            
            RMSE0= obj.RMSE{irf}(obj.RMSE{irf}<= obj.RMSEd);
            rfpwr= zeros(size(RMSE0));
            for idx=1:length(rfpwr)
                rfpwr(idx)= real(X0(:,idx)'*obj.RMatrix{1}'*obj.RMatrix{1}*X0(:,idx));
            end
            
            % find solution that minimizes |R*x|^2 = rf pwr
            rfpwrmin= min(rfpwr);
            jdx0= find(rfpwr==rfpwrmin);
            rmse= RMSE0(jdx0);
            
            obj.Weight{irf}=X0(:,jdx0(rmse==min(rmse)));
            obj.MinRfEnergy{irf}= rfpwrmin;
            
        end
        
        
        function writePulse(obj, rf, grad, irf)
            write_pTXRFPulse(rf,grad,[obj.Tag,obj.Orient],obj.Orient, irf-1, ...
                [], obj.OverSampleFactor);
        end
        
        function writePulseI(obj, rf, grad, soi)
            write_pTXRFPulse(rf,grad,[obj.Tag,obj.Orient,'Sli',num2str(soi)], ...
                obj.Orient, soi-1, [], obj.OverSampleFactor);
        end
        
        
        function [rf, grad]= resamplePulse(obj, irf)
            %
            osf= obj.OverSampleFactor;
            
            grad= obj.Gradient{irf}(:,1:osf:end);
            
            if sum(abs(grad(:,end)))~=0
                grad= [grad, [0;0;0]];
            end
            
            
            glen= size(grad,2);
            blen= glen* osf;
            
            rf= obj.RF{irf};
            rf(:,blen)= 0;
            
        end
        
        function [rf, grad]= resamplePulseI(obj, irf,ib)
            %
            osf= obj.OverSampleFactor;
            
            grad= obj.Gradient{irf}(:,1:osf:end);
            
            if sum(abs(grad(:,end)))~=0
                grad= [grad, [0;0;0]];
            end
            
            
            glen= size(grad,2);
            blen= glen* osf;
            
            rf= obj.iRFs{irf}{ib};
            rf(:,blen)= 0;
            
        end
        
        
        function reset(obj)
            obj.X= {};
            obj.NI= {};
            obj.RHO= {};
            obj.ETA= {};
            obj.RMSE= {};
        end
        
        
    end                                 % protected methods
    
    methods (Static=true)
        function [x,ni,rho,eta,m]= solve_cgmls(A,b,tol,x0,ncgs)
            z = exp(1i*angle(A*x0));
            x = MBpTXadvPulseDesigner.cgls_iterate(A, b.*z, ncgs);
            
            costNew= norm(A*x-b.*z);
            costOld = 1e100; % inf
            ni = 1;
            while ((costOld- costNew)/ costOld > tol),
                ni = ni + 1;
                costOld = costNew;
                
                z = exp(1i*angle(A*x));
                x = MBpTXadvPulseDesigner.cgls_iterate(A, b.*z, ncgs);
                costNew = norm(A*x- b.*z);
            end
            
            rho = norm(abs(A*x) - b);
            eta = norm(x);
            m= A(1:end-size(A,2),:)*x;
        end
        
        function x= cgls_iterate(A, b, nCGs)
            %
            xx= cgls(A,b,nCGs);
            x= xx(:,end);
        end
        
    end
    
end % classdef



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MBpTXadvPulseDesigner.m ends here

