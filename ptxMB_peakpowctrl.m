%%% This is a demo script showing how to design pTx multiband RF pulses with explicit peak power constraint. 
%%% The peak power constraint is formulated based on the composite multiband RF waveforms. 
%%% 
%%% If you find this work is helpful, please consider citing the following ISMRM abstract:
%%%
%%% Wu, X., et al. Peak RF power constrained pulse design for multi-band parallel excitation, ISMRM 2013, p4253
%%% 
%%% Created by Xiaoping Wu on 7/11/2024
%%%

clearvars; 
close all
% %
% load pTxCalibData.mat;%fmapMS % load 3D B1+ and B0 mapping and brain masking
% mask= pTxCalibData.mask;
% b1maps= 1e-6* pTxCalibData.b1maps_uTperV;
% b0map= 1e-6* pTxCalibData.b0map_uT;
% sliceLocations= pTxCalibData.sliceLocations_mm;
% fov= pTxCalibData.FOV_mm;
% res= pTxCalibData.sliceThickness_mm;
% x0= pTxCalibData.phases4flipAngleMapping;
% 
% bDrawROI= false; % true for manually drawing brain masks in slices of interest. 
% %%%%%%%
% nskips=2;
% mask00= mask(1:nskips:end,1:nskips:end, :);
% b1maps00= b1maps(1:nskips:end,1:nskips:end,:,:);
% %b0mapMS= b0map(1:nskips:end,1:nskips:end,whichSlices);
% 
% %%
% soi=[12 36];% MB2
% %soi=[10 30 50 70];%MB4
% %soi=[10 20 30 40 50 60 70 80]; %MB8
% %soi=[20 30 40 50 60 70];%MB6
% mask=mask00(:,:,soi);
% b1maps= b1maps00(:,:,soi,:);
% fox= 1e-3.* [fov(1:2) res*(soi(end)-soi(1))];
% save ptxCalib mask b1maps fox

load ptxCalib.mat
%%
MBptxECdesigner = MBptxECpulseDesigner(b1maps,mask,fox);
MBptxECdesigner.B0Map = [];
MBptxECdesigner.GRasterTime = 10e-6;
MBptxECdesigner.OverSampleFactor= 5;
MBptxECdesigner.SubRFType = 'sinc';
%MBptxECdesigner.ReadOutOffset = roffset;
MBptxECdesigner.SubRFDuration = 1e-3;
MBptxECdesigner.Thickness = 6e-3;
MBptxECdesigner.NominalFOX = (0.2:0.05:0.3);
%MBptxECdesigner.Constraints = {[2]};
MBptxECdesigner.NominalFlipAngle=10;% 
%MBptxECdesigner.RMSEd= inf;
MBptxECdesigner.USFactor=8;

%MBptxECdesigner.Target= targ;
%MBptxECdesigner.ConvergenceTolerance= 0.005;
MBptxECdesigner.PredictionIsNeeded=0;
MBptxECdesigner.NumOfSpokes=1;
MBptxECdesigner.MB= 2; %length(soi);

%MBptxECdesigner.design;
%% 1 spoke, MB2
MBptxECdesigner.ConstraintType='pRFpower';
MBptxECdesigner.PredictionIsNeeded=1;
MBptxECdesigner.NumOfSpokes= 1; %2;
constr1= 5; %1:0.5:5;
clear rmse prfpow
for idx=1:length(constr1)
    MBptxECdesigner.Constraints= {constr1(idx)};
    MBptxECdesigner.design;
     
    prfpow(idx)= max(abs(MBptxECdesigner.RF{1}(:)).^2);
    rmse(idx)=MBptxECdesigner.RMSE{1};
%     result.rf = MBptxECdesigner.RF{1};
%     result.cv= MBptxECdesigner.RMSE{1};% for consistency with the plotting script
%     %result.realCV= mycv;
%     
%     resultfile = createfilename(MBptxECdesigner.NumOfSpokes,...
%         coiltype,'result1',MBptxECdesigner.RMSE{1});
%     save(resultfile,'result');
end

%save prfpowctrl1x8MB2 rmse1x8 prfpow1x8
%%
% plot
figure, plot(rmse,prfpow./50,'bo-')
xlabel('RMSE')
ylabel('peak RF power (a.u.)')
title('peak RF power controlled')


% load pRFpowCtrlNOshape_sp1MB4
% MBptxECdesigner.ConstraintType='pRFpowerNOshape';
% MBptxECdesigner.PredictionIsNeeded=0;
% MBptxECdesigner.NumOfSpokes=1;
% %constr2=0:0.1:5;
% for idx=1:length(constr2)
%     MBptxECdesigner.Constraints= {constr2(idx)};
%     MBptxECdesigner.design;
%     
%     result.rf = MBptxECdesigner.RF{1};
%     result.cv= MBptxECdesigner.RMSE{1};% for consistency with the plotting script
%     %result.realCV= mycv;
%     
%     resultfile = createfilename(MBptxECdesigner.NumOfSpokes,...
%         coiltype,'result1',MBptxECdesigner.RMSE{1});
%     save(['./peakpowCtrlNoShape/',resultfile],'result');
% end
% 

%
% load totRFengyCtrl_sp1MB4
% MBptxECdesigner.ConstraintType='totRFenergy';
% MBptxECdesigner.PredictionIsNeeded=0;
% MBptxECdesigner.NumOfSpokes=1;
% %constr3=0:0.1:7;
% for idx=1:length(constr3)
%     MBptxECdesigner.Constraints= {constr3(idx)};
%     MBptxECdesigner.design;
%     
%     result.rf = MBptxECdesigner.RF{1};
%     result.cv= MBptxECdesigner.RMSE{1};% for consistency with the plotting script
%     %result.realCV= mycv;
%     
%     resultfile = createfilename(MBptxECdesigner.NumOfSpokes,...
%         coiltype,'result1',MBptxECdesigner.RMSE{1});
%     save(['./totalpowCtrl/',resultfile],'result');
%     
% end


% %%  2 spokes, MB=4
% load pRFpowCtrl_sp2MB4
% load pRFpowCtrlNOshape_sp2MB4
% load totRFengyCtrl_sp2MB4
% 
% MBptxECdesigner.NumOfSpokes=2;
% MBptxECdesigner.PredictionIsNeeded=0;
% MBptxECdesigner.ConstraintType='pRFpower';
% %constr1=0.4:0.1:5;
% for idx=1:length(constr1)
%     MBptxECdesigner.Constraints= {constr1(idx)};
%     MBptxECdesigner.design;
%     
%     result.rf = MBptxECdesigner.RF{1};
%     result.cv= MBptxECdesigner.RMSE{1};% for consistency with the plotting script
%     %result.realCV= mycv;
%     
%     resultfile = createfilename(MBptxECdesigner.NumOfSpokes,...
%         coiltype,'result1',MBptxECdesigner.RMSE{1});
%     save(resultfile,'result');
% end
% 
% %
% MBptxECdesigner.ConstraintType='pRFpowerNOshape';
% MBptxECdesigner.PredictionIsNeeded=0;
% %constr2=0:0.1:5;
% for idx=1:length(constr2)
%     MBptxECdesigner.Constraints= {constr2(idx)};
%     MBptxECdesigner.design;
%     
%     result.rf = MBptxECdesigner.RF{1};
%     result.cv= MBptxECdesigner.RMSE{1};% for consistency with the plotting script
%     %result.realCV= mycv;
%     
%     resultfile = createfilename(MBptxECdesigner.NumOfSpokes,...
%         coiltype,'result1',MBptxECdesigner.RMSE{1});
%     save(['./peakpowCtrlNoShape/',resultfile],'result');
% end
% 
% 
% %
% MBptxECdesigner.ConstraintType='totRFenergy';
% %constr3=0:0.1:7;
% for idx=1:length(constr3)
%     MBptxECdesigner.Constraints= {constr3(idx)};
%     MBptxECdesigner.design;
%     
%     result.rf = MBptxECdesigner.RF{1};
%     result.cv= MBptxECdesigner.RMSE{1};% for consistency with the plotting script
%     %result.realCV= mycv;
%     
%     resultfile = createfilename(MBptxECdesigner.NumOfSpokes,...
%         coiltype,'result1',MBptxECdesigner.RMSE{1});
%     save(['./totalpowCtrl/',resultfile],'result');
%     
% end
