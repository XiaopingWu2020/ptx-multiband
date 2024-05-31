
load fmapMS % load in b0 and16 channel b1+ mapping at 7T and masks defining the region of interest in the brain

ns=size(maskMS,3);
fox = 1e-3*[256*[1 0.688] 14*(ns-1)];% in plane field of excitation in m

rmse0=0.02; % root mean squared error specifying the desired excitation fidelity

%% = MB4
soi= [2 4 6 8];%[3 7];

mask1= maskMS(:,:,soi);
b1maps1= b1mapsMS(:,:,soi,:);
b0map1= b0mapMS(:,:,soi);

fox1= [fox(1:2) 1e-3*14*(soi(end)-soi(1))];

MBptxAdvDesigner = MBpTXadvPulseDesigner(b1maps1,mask1,fox1);
MBptxAdvDesigner.B0Map = b0map1;
MBptxAdvDesigner.GRasterTime = 10e-6; % sec
MBptxAdvDesigner.OverSampleFactor= 5;
MBptxAdvDesigner.SubRFType = 'sinc';
MBptxAdvDesigner.ReadOutOffset = 1e-3*[-30 0 0]; % mm
MBptxAdvDesigner.SubRFDuration = 1e-3; % sec
MBptxAdvDesigner.Thickness = 6e-3; % m
MBptxAdvDesigner.NominalFOX = (0.1:0.1:0.35);%(0.1:0.05:0.35);
MBptxAdvDesigner.Lambda = 10.^(-4:2:4); % a heuristically determined range of regularization parameters, over which the solutions will be calculated.
MBptxAdvDesigner.NominalFlipAngle=10;% target flip angle in degrees
MBptxAdvDesigner.ConvergenceTolerance= 0.005; % for internal magnitude least squares optimization
MBptxAdvDesigner.RMSEd= rmse0;% desired excitation fidelity. The solution that satisfies this will be returned.
%MBptxAdvDesigner.Target= targ; % excitation target. defaults to unity in roi
MBptxAdvDesigner.DesignStrategy= 'perb'; % "perb": band specific design; "perm": band joint design
MBptxAdvDesigner.RegularizationType= 'tikhonov';

MBptxAdvDesigner.NumOfSpokes=2;
MBptxAdvDesigner.MB= length(soi);

MBptxAdvDesigner.PredictionIsNeeded=true;% otherwise no numerical prediction will be provided. 

MBptxAdvDesigner.design;

% MBptxAdvDesigner.Tag= fn;
%MBptxAdvDesigner.Orient='tranR2L';
% MBptxAdvDesigner.write;

rfs= MBptxAdvDesigner.RF{1};
figure, plot(abs(rfs).')

grads= MBptxAdvDesigner.Gradient{1};
figure, plot(grads.')

% show numerical prediction
figure, myimagesc(asind(abs(cell2mat(MBptxAdvDesigner.Mxy{1})))),
daspect([1 0.688, 1])

