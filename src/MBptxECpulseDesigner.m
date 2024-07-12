%%% MBptxECpulseDesigner.m --- design MB ptx pulses with explicit constraints.
%%
%% Filename: MBptxECpulseDesigner.m
%% Description:
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu>
%% Maintainer:
%% Copyright (C) 2012 CMRR at UMN
%% Created: Tue Oct 30 12:44:20 2012 (CDT)
%% Version:
%% Last-Updated: Tue Jul 30 15:46:09 2013 (CDT)
%%           By: Xiaoping Wu
%%     Update #: 55
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


classdef MBptxECpulseDesigner < handle
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
        Constraints = {[1:0.2:5]}
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
        ConstraintType= 'pRFpower'
        USFactor= 8
        MaxNumOfRWIters= 100
        PredictionIsNeeded= true
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
        X= {}
        RMSE= {}
        Mxy= {}
    end % properties
    
    
    methods
        function obj = MBptxECpulseDesigner (b1map, mask, fox)
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
        
        function obj= set.Constraints(obj, constr)
            %
            if ~isequal(constr,obj.Constraints)
                obj.Constraints= constr;
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
        
        
        
        function obj= set.ConstraintType(obj, regtype)
            %
            if ~(strcmpi(regtype,'pRFpower')||strcmpi(regtype,'pRFpowerNOshape')|| ...
                 strcmpi(regtype,'totRFenergy')|| ...
                    strcmpi(regtype,'gsar')||strcmpi(regtype,'lsar')|| ...
                    strcmpi(regtype,'gsarNoRFshape')||strcmpi(regtype,'lsarNoRFshape'))
                disp(['-> ConstraintType must be pRFpower, pRFpowerNOshape, ' ...
                      'or totRFenergy. set to current type.'])
                regtype = obj.ConstraintType;
            end
            
            oldtype= obj.ConstraintType;
            obj.ConstraintType = regtype;
            
            if ~strcmpi(oldtype,obj.ConstraintType)
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
            
            if ~(strcmpi(rftype,'gauss')||strcmpi(rftype,'sinc')||strcmpi(rftype,'sinc8')||...
                    strcmpi(rftype,'sinc6')||strcmpi(rftype,'sinc4'))
                disp('-> SubRFType must be gauss or sinc. set to gauss')
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
            
            obj.designPERB();
            
            
            disp('-> Pulse calculations done...')
            
        end
        
        
    end % methods
    
    %% protected methods
    methods (Access = protected)
        
        function designPERB(obj)
            
            if ~obj.IsOptimal
                obj.findOptimalKspace();
                obj.assembleGradient();
            end
            
            obj.calculateWeights();
            
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
                case 'sinc'
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
            
            obj.SliceLocation= linspace(-foxz/2,foxz/2,ns);
            
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
                    
                    rho= obj.solveCVX(0);
                    
                    myrho(ifox,jtheta) = rho;
                    mykp{ifox,jtheta} = kp;
                end
            end
            
            [ifox,jtheta] = find(myrho==min(myrho(:)));
            kp0= mykp{ifox(1),jtheta(1)};
            
        end
        
        
        function reshapeWeight(obj)
            
            for idx=1:length(obj.Weight),
                wts= obj.Weight{idx};
                obj.Weight{idx}= reshape(wts(:),[],obj.MB);
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
            switch obj.ConstraintType
                case 'pRFpower'
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
            
            regtype= obj.ConstraintType;
            if isequal(regtype,'pRFpower')||isequal(regtype,'pRFpowerNOshape')||...
                  isequal(regtype,'totRFenergy')
                return;
            end
            
            if isempty(obj.SARMatrix)
                error('-> please specify obj.SARMatrix for sar constraint...')
            end
            
            if (isequal(regtype,'gsar')||isequal(regtype,'lsar'))
                
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
            obj.mVector= obj.mVectors{irf};            
        end
        
        
        function buildRMatrix (obj,irf)
            % regularization matrix construction for perb
            obj.RMatrix= {};
            regtype= obj.ConstraintType;
            
            if isequal(regtype,'pRFpowerNOshape')||isequal(regtype,'totRFenergy')
                
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
            
            obj.calcWeights();
        end
        
        
        function calcWeights (obj)
            % Purpose:
            if ~obj.IsCalculated
                
                obj.reset();
                
                for irf=1:length(obj.BandPattern),
                    obj.constructSystemMatrixMB(obj.OptimalKspace{irf},irf);
                    
                    obj.buildRMatrix(irf);
                    
                    obj.createM(irf);
                    
                    obj.solve(irf);
                end
                
                obj.setCalcStatus(false);
                
            end
            
        end
        
        
        function solve (obj,irf)
            %
            for jdx=1:length(obj.Constraints{1})
                x= obj.solveCVX(jdx);
                rmse= obj.calcRMSE(x);
                
                obj.X{irf}(:,jdx)= x;
                obj.RMSE{irf}(jdx)= rmse;
            end
            
            obj.Weight{irf}= obj.X{irf}(:,1);
            
        end
        
        
        function x = solveCVX (obj,idx)
            % solve least squares minimization with explicit constraints using cvx.
            
            A= obj.AMatrix./max(abs(obj.AMatrix(:)));
            n=size(A,2);
            
            if idx==0 % calc rho
                cvx_begin %quiet
                cvx_precision high
                variable x(n) complex;
                minimize(norm(A*x- 1));
                cvx_end
                
                x = norm(A*x - 1);
                return;
                
            end
            
            myconstr= obj.Constraints{1}(idx);
            
            switch obj.ConstraintType
                case 'pRFpower'
                    R= obj.RMatrix{1};
                    
                    cvx_begin %quiet
                    cvx_precision high
                    variable x(n) complex;
                    minimize(norm(A*x- 1));
                    subject to
                    norm(R*x,Inf)<= myconstr;
                    cvx_end
                    
                case 'pRFpowerNOshape'
                    cvx_begin %quiet
                    cvx_precision high
                    variable x(n) complex;
                    minimize(norm(A*x- 1));
                    subject to
                    norm(x,Inf)<= myconstr;
                    cvx_end
                    
                case 'totRFenergy'
                    cvx_begin %quiet
                    cvx_precision high
                    variable x(n) complex;
                    minimize(norm(A*x- 1));
                    subject to
                    norm(x)<= myconstr;
                    cvx_end
                    
                otherwise
                    
            end

            x= x.*sind(obj.NominalFlipAngle)./max(abs(obj.AMatrix(:)));
            
        end
        
        
        function rmse= calcRMSE(obj, x)
            
            % scale x
            
            if obj.NominalFlipAngle<=30
                rmse= sqrt(sum( abs( obj.AMatrix*x- sind(obj.NominalFlipAngle) ...
                                     ).^2 )./size(obj.AMatrix,1));
            else
                rmse= obj.simulateRMSE(x);
            end
            
        end
        % -----------------------
        
        
        function rmse= simulateRMSE(obj, wts)
            % currently only works for the case of single MB pulse design.
            irf=1;
            
            % reshape wts
            
            wts= reshape(wts(:),[],obj.MB);
            
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
            
            rmse= sqrt(sum(abs(mxy- obj.mVector).^2) ./ ...
                length(obj.mVector));
            
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
            obj.RMSE= {};
        end
        
        
    end                                 % protected methods
    
end % classdef



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MBptxECpulseDesigner.m ends here
