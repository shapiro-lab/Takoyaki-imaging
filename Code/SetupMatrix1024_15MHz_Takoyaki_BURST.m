%% Takoyaki BURST script
%  This script scans multiple focal points simultaneously using the 1024 element 15 MHz 
%  2D matrix array probe from Vermon with the UTA 1024-MUX adapter on a 
%  Vantage 256 system, integrated with BURST scheme. All 4 banks are used 
%  for transmission. For receive apertures, only one random sparse aperture
%  is used for each transmit aperture.
%
%  only support P.nCode = 1, Takoyaki Bmode + BURST
%  Note that the initialization can take a long time. (~30 sec * #
%  of frames) After initialization, please change "High Voltage P2" to your
%  desired high voltage value, and then press the button "Record BURST" to 
%  begin data acquisition. IQ data are accumulated across frames. 
%
%  To run this script, there are two required files:
%  computeTrans_MatrixProbeIncluded.m
%  saved100CompRndApod.mat - also exists in
%  Vantage folder/Example_Scripts/Biomedical/Vantage 128 and
%  256/UTA-1024-MUX/Matrix_3MHz/SparseApertures folder
%
% written by Sunho Lee

clear;

% Parameter settings
P.pathName = '/Volumes/T7/Caltech/Sunho' ; 

P.startDepth_mm = 3.3;   % 
P.endDepth_mm = 10;
P.gaussFiltSigma = 0; % parameter for gaussian filter in CACTUS display    
P.Voltage = 1.6 ;

P.numPreColFrames = 2; % no. of frames before collapse transmits
P.numColFrames = 18; % no. of collapse transmits
P.numPostColFrames = 0; % no. of frames after collapse transmits

P.numSyntheticRcvs = 1; % will only use 1 rnd sparse aperture for 1 TX
P.numHalfCycles =3;

P.txFocus_mm = 7; % focal depth
P.aperpatch = 8;
P.numRays = 32 - P.aperpatch + 1; % no. of Rays (for each dimension)
P.nCode = 1; % 3 = 1 full amplitude  + 2 half amplitude

Display_output = 20; %Display mode only: 20 gives a power of 1.0 for a linear output, 40 gives a power of 0.5 for a square root compression
RcvProfile.AntiAliasCutoff = 15; % 10, 15, 20, 30

P.iV = 1:P.numPreColFrames+P.numPostColFrames+P.numColFrames;

%% Specify system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;   % number of receive channels.
Resource.Parameters.verbose = 2; % 2: warnings and status 1: warnings 0: errors only
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.fakeScanhead = 0; % allow system HW operation with nothing connected
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.Parameters.Connector = 1;

%% Specify Trans structure array.
Trans.name = 'Matrix1024-15';
Trans.units = 'mm';
Trans = computeTrans_MatrixProbeIncluded(Trans);
Trans.HVMux = computeUTAMux1024; % add the HV Mux programming for the UTA 1024-MUX

% this files contains 100 sets of 4 Rnd Apertures that are complementary. 
% CompRndApod is matrix of 100x4x1024. The sum of every 4 matrices (of the 
% second dimention) account for the full aperture. 
load('saved100CompRndApod'); %This file has 100x4x1024

% patch-like aperture
Papod = zeros(P.numRays^2, 1+P.numSyntheticRcvs, 1024);

for i = 1:P.numRays
    for j = 1:P.numRays
        
        center1 = floor(P.aperpatch/2)+i-1;
        center2 = floor(P.aperpatch/2)+j-1;

        M = zeros(32); % initialize 32 x 32 matrix 
        M(center1-floor(P.aperpatch/2)+1:center1+ceil(P.aperpatch/2),...
            center2-floor(P.aperpatch/2)+1:center2+ceil(P.aperpatch/2)) = 1;
    
        Papod((i-1)*P.numRays + j, 1, :) = reshape(M, [1, 1024]);
        Papod((i-1)*P.numRays + j,2:1+P.numSyntheticRcvs,:) = CompRndApod(mod((i-1)*P.numRays + j -1,size(CompRndApod,1))+1,1:P.numSyntheticRcvs,:);

    end
end

Papod2 = Papod;
Papod(:,1,:) = ones(size(Papod,1),1,1024);

% create aperture index for each Papod
for i = 1:size(Papod, 1)
    for j = 1:size(Papod, 2)
        Paper(i,j) = computeMuxAperture(squeeze(Papod(i,j,:))', Trans);
    end
end

% Intermediate Variables
waveLength = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
if strcmp(Trans.units,'mm')
    Trans.ElementPosMm = Trans.ElementPos;
    Trans.ElementPosWL = Trans.ElementPos./waveLength;
else
    Trans.ElementPosMm = Trans.ElementPos.*waveLength;
    Trans.ElementPosWL = Trans.ElementPos;
end

ElementPosx = unique(Trans.ElementPosWL(:,1));
ElementPosy = unique(Trans.ElementPosWL(:,2)); % y direction has a gap

% Convert mm to wavelength
demodFreq = Trans.frequency; % demodulation frequency
P.startDepth =  P.startDepth_mm/waveLength;	% Acquisition depth in wavelengths
P.endDepth =    P.endDepth_mm/waveLength;    % This should preferrably be a multiple of 128 samples
P.txFocus = P.txFocus_mm/waveLength;
%% PData
PData(1).PDelta = [1,1,1]*1/2;

PData(1).Size(1) = ceil(((P.numRays+3)*Trans.spacing)/PData(1).PDelta(1)); % +3 compensates for gaps
PData(1).Size(2) = ceil((P.numRays*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(1), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), P.startDepth];
P.lateral_voxel_size = PData(1).PDelta(1)*100e-6 ;
P.axial_voxel_size = PData(1).PDelta(3)*100e-6 ;

% specify Region structures.
PData(1).Region = repmat(struct('Shape',struct( ...
    'Name','Frustum',...
    'Position',[0,0,P.startDepth],...
    'l1',Trans.spacing,...
    'w1',Trans.spacing,...
    'height',P.endDepth-P.startDepth,...
    'l2',Trans.spacing,...
    'w2',Trans.spacing)),1,P.numRays^2);

for i = 1:P.numRays
    for j = 1:P.numRays
        PData.Region((i-1)*P.numRays + j).Shape.Position(1) = (ElementPosx(i)+ElementPosx(i+P.aperpatch-1))/2;
        PData.Region((i-1)*P.numRays + j).Shape.Position(2) = -(ElementPosy(j)+ElementPosy(j+P.aperpatch-1))/2;
    
        if ismember(j, [2:8:P.numRays])
            PData.Region((i-1)*P.numRays + j).Shape.Position(2) = -(ElementPosy(j)+ElementPosy(j+8-1)-Trans.spacing/2)/2;
            PData.Region((i-1)*P.numRays + j).Shape.w1 = 1.5*Trans.spacing;
            PData.Region((i-1)*P.numRays + j).Shape.w2 = 1.5*Trans.spacing;
        end

        if ismember(j, [8:8:P.numRays])
            PData.Region((i-1)*P.numRays + j).Shape.Position(2) = -(ElementPosy(j)+ElementPosy(j+8-1)+Trans.spacing/2)/2;
            PData.Region((i-1)*P.numRays + j).Shape.w1 = 1.5*Trans.spacing;
            PData.Region((i-1)*P.numRays + j).Shape.w2 = 1.5*Trans.spacing;
        end

    end
end

PData(1).Region = computeRegions(PData(1));

for i = 1:8
    for j = 1:8

        PData.Region(P.numRays^2 + (i-1)*8 + j).Shape.Name = 'Custom';

        center_idx = (j:8:P.numRays)' + P.numRays*((i:8:P.numRays)-1); center_idx = center_idx(:);
        pixels_set = PData.Region(center_idx(1)).PixelsLA;
        for set_n = 2:length(center_idx)
            pixels_set = union(pixels_set, PData.Region(center_idx(set_n)).PixelsLA);
        end

        PData.Region(P.numRays^2 + (i-1)*8 + j).PixelsLA = pixels_set;
        PData.Region(P.numRays^2 + (i-1)*8 + j).numPixels = length(PData.Region(P.numRays^2 + (i-1)*8 + j).PixelsLA);
    
    end
end

PData.Region = PData.Region(end-8^2+1:end);

%% Specify Media.  Use point targets in middle of PData.
Media.MP(1,:) = [0,0,60,1.0];      % single point.

Media.numPoints = size(Media.MP,1);

%% Resources
% Resource.VDAS.el2ChMapDisable = 1;

Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*8^2*P.nCode*P.numSyntheticRcvs ;
Resource.RcvBuffer(1).colsPerFrame = 256;
Resource.RcvBuffer(1).numFrames = P.numPreColFrames+P.numPostColFrames+P.numColFrames;
Resource.ImageBuffer(1).numFrames = P.numPreColFrames+P.numPostColFrames+P.numColFrames;
% Resource.ImageBuffer(2).numFrames = P.numFrames;

Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;

% Encode numbers of frames in less cumbersome way
P.nRcvFrms = Resource.RcvBuffer(1).numFrames;
P.nImgFrms = Resource.ImageBuffer(1).numFrames;

Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.25;
Resource.DisplayWindow(1).Position = [0,580, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0.0,PData(1).Origin(3)];
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).AxesUnits = 'mm';

% Resource.DisplayWindow(2).mode = '2d';
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).Title = '3D Flash Image - YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [430,580, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),PData(1).Origin(3)];
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).AxesUnits = 'mm';

% Resource.DisplayWindow(3).mode = '2d';
Resource.DisplayWindow(3).Type = 'Verasonics';
Resource.DisplayWindow(3).Title = '3D Flash Image - XY plane';
Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(3).Position = [0,40, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
Resource.DisplayWindow(3).Orientation = 'xy';
Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),-PData(1).Origin(2),60];%PData.Region(end).Shape.oPAIntersect];
Resource.DisplayWindow(3).Colormap = gray(256);
Resource.DisplayWindow(3).AxesUnits = 'mm';
% Resource.DisplayWindow(3).mode = '2d';

%% Specify Transmit waveform structure.
TW.type = 'parametric';
% 15.625 MHz center frequency, 67% pulsewidth 1.5 cycle burst
TW.Parameters = [Trans.frequency,.67,P.numHalfCycles,1];  % A, B, C, D
P.TW = TW;

% Specify TPC structure.
TPC(1).hv = P.Voltage;

%% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', P.txFocus, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,1024), ...
    'Delay', zeros(1,256),...
    'peakCutOff', 2,...
    'peakBLMax', 20,...
    'Angle', 0,... 
    'aperture',1,...
    'TXPD',[]), 1,P.numRays^2);

for i=1:P.numRays
    for j = 1:P.numRays
        
        TX((i-1)*P.numRays+j).Apod = squeeze(Papod2((i-1)*P.numRays+j,1,:))';
        TX((i-1)*P.numRays+j).aperture = Paper((i-1)*P.numRays+j,1);
        TX((i-1)*P.numRays+j).Origin = [(ElementPosx(i)+ElementPosx(i+P.aperpatch-1))/2, -(ElementPosy(j)+ElementPosy(j+P.aperpatch-1))/2, 0];
        TX((i-1)*P.numRays+j).Delay = computeTXDelays(TX((i-1)*P.numRays+j));

%         TX((i-1)*P.numRays+j).TXPD = computeTXPD(TX((i-1)*(P.numRays-1)+j),PData(1));   

    end
end

k = P.numRays^2;
Paper_real = zeros(1,64);

for i = 1:8
    for j = 1:8
        
        delay_matrix = reshape(TX(j+25*(i-1)).Delay, [32 32]);
        delay_matrix = repmat(delay_matrix(i:i+7,j:j+7),4,4);
        translated_D = reshape(circshift(delay_matrix, [i-1, j-1]), 1, 1024);

        Ap = ones(32);
        if i > 1; Ap([1:i-1, 32-8+i:32], :) = 0; end
        if j > 1; Ap(:, [1:j-1, 32-8+j:32]) = 0; end

        Ap = reshape(Ap, 1, 1024);
        
        TX(k+j+8*(i-1)) = TX(1);
        TX(k+j+8*(i-1)).Delay = translated_D;
        TX(k+j+8*(i-1)).focus = 0.0;
        TX(k+j+8*(i-1)).Origin = [0.0, 0.0, 0.0];
        TX(k+j+8*(i-1)).Apod = Ap;
        
        Paper_real(k+(i-1)*8) = computeMuxAperture(Ap, Trans);
        TX(k+j+8*(i-1)).aperture = Paper_real(k+(i-1)*8);
        
        if P.nCode > 1
            TX(k+(i-1)*8+j+8^2) = TX(k+j+8*(i-1));
            TX(k+(i-1)*8+j+2*8^2) = TX(k+j+8*(i-1));

            mask = repmat([repmat([1 0], 1, 16) repmat([0 1], 1, 16)], 1, 16);
            TX(k+(i-1)*8+j+8^2).Apod = TX(k+j+8*(i-1)).Apod.*mask;
            TX(k+(i-1)*8+j+2*8^2).Apod = TX(k+j+8*(i-1)).Apod.*(~mask);
        end

    end
end

TX = TX(end-8^2*P.nCode+1:end);
TX(end+1) = struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', 0.0, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,1024), ...
    'Delay', zeros(1,1024),...
    'peakCutOff', 2,...
    'peakBLMax', 20,...
    'Angle', 0,... 
    'aperture',1,...
    'TXPD',[]);
nDummyTX = P.nCode*8^2 + 1;

% showTXPD
%% allows to visualize the TX in 3D (similar to the EventAnalysisTool
% figure;
% 
% for i=129
%     plot3(Trans.ElementPos(find(TX(i).Apod==1),1),Trans.ElementPos(find(TX(i).Apod==1),2),TX(i).Delay(find(TX(i).Apod==1)),'k*');
%     zlim([0 5]);
%     xlim([min(Trans.ElementPos(:,1)) max(Trans.ElementPos(:,1))]);
%     ylim([min(Trans.ElementPos(:,2)) max(Trans.ElementPos(:,2))]);
%     view(gca,-40,49);grid on;
% %     pause(.5);
%     zlabel('Delay (wavelength)'); xlabel('X (mm)'); ylabel('Y (mm)')
% end

%% Specify TGC Waveform structure.
TGC(1).CntrlPts = [100 1010 1023 1023 1023 1023 1023 1023];
TGC(1).rangeMax = 128;
TGC(1).Waveform = computeTGCWaveform(TGC);
maxAcqLength = sqrt(P.endDepth^2 + (32*Trans.spacing)^2 +(35*Trans.spacing)^2 );

%% Receive
Receive = repmat(struct(...
    'Apod', zeros(1,Trans.numelements), ...
    'startDepth', P.startDepth, ...
    'endDepth', 128/(4*2)*ceil(maxAcqLength/(128/(4*2))), ...
    'TGC', 1, ...
    'mode', 0, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'aperture',0, ...
    'sampleMode', 'NS200BW', ...
    'callMediaFunc', 0), 1, P.nRcvFrms*P.numSyntheticRcvs*P.nCode*8^2);

for j = 1:P.nRcvFrms
    k=P.numSyntheticRcvs*8^2*P.nCode*(j-1);
    for n = 1:P.nCode
        kk = P.numSyntheticRcvs*8^2*(n-1);
        for i = 1:8^2
            kkk=P.numSyntheticRcvs*(i-1);
    
            for w=1:P.numSyntheticRcvs
                Receive(w+kkk+kk+k).acqNum = w+kkk;
                Receive(w+kkk+kk+k).framenum = j;
                Receive(w+kkk+kk+k).aperture = Paper(i,w+1);
                if n==1
                    Receive(w+kkk+kk+k).Apod = squeeze(Papod(i,w+1,:))';
                else 
                    Receive(w+kkk+kk+k).Apod = -squeeze(Papod(i,w+1,:))'; % subtract
                    Receive(w+kkk+kk+k).mode = 1; % accumulate
                end
            end
        end
    end
end

%% Recon
% - We need one Recon structure for each frame. 
senscutoff = 0.6;

Recon = repmat(struct('senscutoff', senscutoff, ...
                   'pdatanum', 1, ...
                   'rcvBufFrame', -1, ...
                   'IntBufDest', [1,1],...
                   'ImgBufDest', [1,-1], ...
                   'RINums', 1:8^2*P.numSyntheticRcvs), 1, P.nImgFrms+1) ;

for i = 1:P.nImgFrms+1
    Recon(i).RINums = (i-1)*8^2*P.numSyntheticRcvs + (1:8^2*P.numSyntheticRcvs);
    if i <= P.nImgFrms
        Recon(i).rcvBufFrame = i;
        Recon(i).ImgBufDest = [1,i];
    end
end

%% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...
    'txnum', 1, ...
    'rcvnum', 1, ...
    'scaleFactor', 1, ...
    'regionnum', 1), 1, 8^2*P.numSyntheticRcvs*(P.nImgFrms+1));


ReconInfo(1).mode = 3;

for j = 1:P.nImgFrms+1
    kk = (j-1)*8^2*P.numSyntheticRcvs;
    for w=1:8^2
        k= P.numSyntheticRcvs*(w-1);
        for i = 1:P.numSyntheticRcvs
            ReconInfo(i+k+kk).txnum = w;
            ReconInfo(i+k+kk).rcvnum = i+k;
            ReconInfo(i+k+kk).regionnum = w;
            if i == 1
                ReconInfo(i+k+kk).mode = 3;
            end
            if i == P.numSyntheticRcvs
                ReconInfo(i+k+kk).mode = 5;
            end
        end    
    end
end

ReconInfo(end).mode = 5;

%% Specify Process structure arrays.
pgain=20;

Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...
                         'framenum',-1,...
                         'pdatanum',1,...
                         'srcData','intensity3D',...
                         'pgain', pgain,...
                         'persistMethod','none',...
                         'persistLevel',30,...
                         'interpMethod','4pt',...
                         'compressMethod','power',...
                         'compressFactor',75,...
                         'mappingMethod','full',...
                         'display',1,...
                         'displayWindow',1};
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',1,...
                         'framenum',-1,...
                         'pdatanum',1,...
                         'srcData','intensity3D',...
                         'pgain', pgain,...
                         'persistMethod','none',...
                         'persistLevel',30,...
                         'interpMethod','4pt',...
                         'compressMethod','power',...
                         'compressFactor',75,...
                         'mappingMethod','full',...
                         'display',1,...
                         'displayWindow',2};
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,...
                         'framenum',-1,...
                         'pdatanum',1,...
                         'srcData','intensity3D',...
                         'pgain', pgain,...
                         'persistMethod','none',...
                         'persistLevel',30,...
                         'interpMethod','4pt',...
                         'compressMethod','power',...
                         'compressFactor',75,...
                         'mappingMethod','full',...
                         'display',1,...
                         'displayWindow',3};
Process(4).classname = 'External';
Process(4).method = 'volumetricPlot';
Process(4).Parameters = {'srcbuffer','image',...
                         'srcbufnum',1,...
                         'srcframenum',-1,...
                         'dstbuffer','none'};

nproc_resetStartEvent = 5;
Process(5).classname = 'External';
Process(5).method = 'resetStartEvent';
Process(5).Parameters = {'srcbuffer','none','dstbuffer','none'};

nproc = 6;
nImageSaveProcess = nproc;
for N = 1:P.nImgFrms
    Process(nproc).classname = 'External';                   
    Process(nproc).method = 'saveImageData';
    Process(nproc).Parameters = {'srcbuffer','image',...
        'srcbufnum',1,...
        'srcframenum',N,... 
        'dstbuffer','none'};
    nproc = nproc+1;
end

nImageRecon = nproc;
for N = 1:P.nImgFrms
    Process(nproc).classname = 'Image';
    Process(nproc).method = 'imageDisplay';
    Process(nproc).Parameters = {'imgbufnum',1,...   % no. of buffer to process.
        'framenum',N,...   
        'pdatanum',1,...    % no. of PData structure to use
        'pgain',1,...            % pgain is image processing gain
        'reject',2,...      % reject level
        'persistMethod','simple',...
        'persistLevel',0,...
        'interpMethod','4pt',...
        'grainRemoval','none',...
        'processMethod','none',...
        'averageMethod','none',...
        'compressMethod','power',...
        'compressFactor',40,...
        'mappingMethod','full',...
        'display',1,...      % display image after processing
        'displayWindow',1};
    nproc = nproc+1;
end

% nproc_tic = nproc;
% Process(nproc).classname = 'External';
% Process(nproc).method = 'vstic';
% Process(nproc).Parameters = {'srcbuffer','none'};
% nproc = nproc+1;
% 
% nproc_toc = nproc;
% Process(nproc).classname = 'External';
% Process(nproc).method = 'vstoc';
% Process(nproc).Parameters = {'srcbuffer','none'};
% nproc = nproc+1;

% Save RF Data
nRFSaveProcess = nproc;
for N = 1:P.nImgFrms
    Process(nproc).classname = 'External';                   
    Process(nproc).method = 'saveRFData';
    Process(nproc).Parameters = {'srcbuffer','receive',...
        'srcbufnum',1,...
        'srcframenum',N,... 
        'dstbuffer','none'};
    nproc = nproc+1;
end               

%external function
EF(1).Function = vsv.seq.function.ExFunctionDef('volumetricPlot', @volumetricPlot);
EF(2).Function = vsv.seq.function.ExFunctionDef('resetStartEvent', @resetStartEvent);
EF(3).Function = vsv.seq.function.ExFunctionDef('saveImageData', @saveImageData);
EF(4).Function = vsv.seq.function.ExFunctionDef('setVoltage', @setVoltage);
EF(5).Function = vsv.seq.function.ExFunctionDef('resetVoltage', @resetVoltage);
EF(6).Function = vsv.seq.function.ExFunctionDef('vstic', @vstic);
EF(7).Function = vsv.seq.function.ExFunctionDef('vstoc', @vstoc);
EF(8).Function = vsv.seq.function.ExFunctionDef('saveRFData', @saveRFData);

%% Specify SeqControl structure arrays.
nsc = 1;
nsc_jump = nsc;
SeqControl(nsc).command = 'jump'; % jump back to start
SeqControl(nsc).argument = 1;
nsc = nsc+1;

nsc_acqDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(nsc).argument = 3*P.endDepth_mm/Resource.Parameters.speedOfSound*2*1e3; % 250
nsc = nsc+1;

nsc_frameDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';  % time between frames
SeqControl(nsc).argument = 20000  - P.nCode*(64)*SeqControl(nsc_acqDelay).argument;  % 20 msec
% SeqControl(nsc).argument = 1000000;  % 2s
SeqControl(nsc).condition = 'ignore'; %Recon processing time
nsc = nsc+1;

nsc_TPCDelay = nsc;
SeqControl(nsc).command = 'timeToNextAcq';  % time to change voltage
SeqControl(nsc).argument = 6000;  % 6 msec
nsc = nsc+1;

nsc_setCollapseVoltage = nsc;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).argument = 2;
SeqControl(nsc).condition = 'next';
nsc = nsc+1;

nsc_setImagingVoltage = nsc;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).argument = 1;
SeqControl(nsc).condition = 'next';
nsc = nsc+1;

nsc_returnToMatlab = nsc;
SeqControl(nsc).command = 'returnToMatlab';
nsc = nsc+1;

% nsc_triggerOut = nsc;
% SeqControl(nsc).command = 'triggerOut';
% nsc = nsc+1;

nsc_sync = nsc;
SeqControl(nsc).command = 'sync'; % - Synchronize hardware and software sequencers
% SeqControl(nsc).argument = 100000; % 100 msec timeout for software sequencer (default is 0.5 seconds)
nsc = nsc+1;

% nsc is count of SeqControl objects
%% Specify Event Structure arrays
n = 1;   % start index for Events
lastTTHnsc = 0;
nscStart = nsc;

% standard live imaging
for j = 1:1 % P.nRcvFrms

    v=P.numSyntheticRcvs*8^2*P.nCode*(j-1);
    for w=1:8^2
        k= P.numSyntheticRcvs*(w-1);
        for i = 1:P.numSyntheticRcvs
            for m = 1:P.nCode
                Event(n).info = 'Transmit/Receive.';
                Event(n).tx = w + (m-1)*8^2;   % % Loop over the TX structure.
                Event(n).rcv = i+k+v+(m-1)*8^2*P.numSyntheticRcvs;  % Rcv structure of frame.
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = nsc_acqDelay; %
                n = n+1;
            end
        end
    end
    Event(n-1).seqControl = [nsc_frameDelay]; % nsc_sync nsc

    Event(n).info = 'Return to Matlab';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [nsc_returnToMatlab nsc_sync];
    n = n+1;
end
% change first seqControl to point to last TTH from acquisition loop 
% SeqControl(nscStart).argument = lastTTHnsc;

Event(n).info = 'Jump';
Event(n).tx = 0;         % no TX structure.
Event(n).rcv = 0;        % no Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % call processing function
Event(n).seqControl = nsc_jump; %
n = n+1;

% Acquire time series collapse data
nWideBeam = n;

TTHnsc = nan(1,P.nImgFrms);

for j = 1:P.nImgFrms
    
    v=P.numSyntheticRcvs*8^2*P.nCode*(j-1);
    for w=1:8^2
        k= P.numSyntheticRcvs*(w-1);
        for i = 1:P.numSyntheticRcvs
            for m = 1:P.nCode
                Event(n).info = 'Transmit/Receive.';
                Event(n).tx = w + (m-1)*8^2;   % % Loop over the TX structure.
                Event(n).rcv = i+k+v+(m-1)*8^2*P.numSyntheticRcvs;  % Rcv structure of frame.
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = nsc_acqDelay; %
                n = n+1;
            end
        end
    end
    Event(n-1).seqControl = [nsc_frameDelay nsc]; % modify last acquisition Event's seqControl
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    
    if (j == P.numPreColFrames)
        Event(n).info = 'set voltage';
        Event(n).tx = nDummyTX;
        Event(n).seqControl = [nsc_TPCDelay nsc_setCollapseVoltage];
        n = n+1;
    elseif (j == P.numPreColFrames+P.numColFrames)
        Event(n).info = 'reset voltage';
        Event(n).tx = nDummyTX;
        Event(n).seqControl = [nsc_TPCDelay nsc_setImagingVoltage];
        n = n+1;
    end
end

for i = 1:P.nImgFrms
    
    Event(n).info = ['3D Reconstruct Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = i;
    Event(n).process = nImageRecon+i-1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Save image data'; 
    Event(n).process = i+nImageSaveProcess-1;    % processing
    n = n+1;
    
%     Event(n).info = 'Save RF data'; 
%     Event(n).process = i+nRFSaveProcess-1;    % processing
%     n = n+1;
end

Event(n).info = 'Reset start event';
Event(n).process = nproc_resetStartEvent;
n = n+1;
% 
Event(n).info = 'Return to Matlab';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc_returnToMatlab;
n = n+1;

% Make unspecified event field values zero by default
for i = 1:length(Event)
    if isempty(Event(i).tx); Event(i).tx = 0; end
    if isempty(Event(i).rcv); Event(i).rcv = 0; end
    if isempty(Event(i).recon); Event(i).recon = 0; end
    if isempty(Event(i).process); Event(i).process = 0; end
    if isempty(Event(i).seqControl); Event(i).seqControl = 0; end
end

%% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl;
m = 1;

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB4','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback );

% - Section number change
UI(2).Control = VsSliderControl('LocationCode','UserB1','Label','Z-Section',...
                  'SliderMinMaxVal',[1,P.endDepth,Resource.DisplayWindow(3).ReferencePt(3)],...
                  'SliderStep',[1/(P.endDepth-1) 10/(P.endDepth-1)],'ValueFormat','%3.0f', ...
                  'Callback', @ZSectionChange );

% - Section number change
UI(3).Control = VsSliderControl('LocationCode','UserB2','Label','Y-Section',...
                  'SliderMinMaxVal',[PData(1).Origin(2)-(PData(1).Size(1)-1)*PData(1).PDelta(1), PData(1).Origin(2), Resource.DisplayWindow(1).ReferencePt(2)],...
                  'SliderStep',[1/((PData(1).Size(2)-1)*PData(1).PDelta(2)-1) 10/((PData(1).Size(2)-1)*PData(1).PDelta(2)-1)],'ValueFormat','%3.0f', ...
                  'Callback', @YSectionChange);

% - Section number change
UI(4).Control = VsSliderControl('LocationCode','UserB3','Label','X-Section',...
                  'SliderMinMaxVal',[PData(1).Origin(1), PData(1).Origin(1)+(PData(1).Size(2)-1)*PData(1).PDelta(2), Resource.DisplayWindow(2).ReferencePt(1)],...
                  'SliderStep',[1/((PData(1).Size(1)-1)*PData(1).PDelta(1)-1) 10/((PData(1).Size(1)-1)*PData(1).PDelta(1)-1)],'ValueFormat','%3.0f', ...
                  'Callback', @XSectionChange);

% - compression
CompressionFactor=1.3;
UI(5).Control = VsSliderControl('LocationCode','UserA1','Label','CompressionP',...
                  'SliderMinMaxVal',[1,11,CompressionFactor],...
                  'SliderStep',[1/40,4/40],'ValueFormat','%1.1f', ...
                  'Callback', @CompressionP );

m = 6;
% - Record burst
UI(m).Control = {'UserC7','Style','VsPushButton','Label','Record Burst'};
UI(m).Callback = text2cell('%wideBeam');
m = m+1;

% - Beam/mode change
UI(m).Control = {'UserB8','Style','VsButtonGroup','Title','BURST mode',...
    'NumButtons',2,'Labels',{'BURST','BURST+'}};
UI(m).Callback = text2cell('%BURSTmodeCallback');
m = m + 1;

% - Change start/stop/step
UI(m).Control = {'UserB7','Style','VsSlider','Label','Pre-coll frames',...
    'SliderMinMaxVal',[0,100,P.numPreColFrames],...
    'SliderStep',[0.01,0.1], 'ValueFormat', '%3.0i'};
UI(m).Callback = text2cell('%preBURST');
m = m + 1;

UI(m).Control = {'UserB6','Style','VsSlider','Label','Coll frames',...
    'SliderMinMaxVal',[0,100,P.numColFrames],...
    'SliderStep',[0.01,0.1],'ValueFormat','%3.0i'};
UI(m).Callback = text2cell('%durBURST');
m = m + 1;

UI(m).Control = {'UserB5','Style','VsSlider','Label', 'Post-coll frames',...
    'SliderMinMaxVal',[0,100,P.numPostColFrames],...
    'SliderStep',[0.01,0.1],'ValueFormat','%3.0i'};
UI(m).Callback = text2cell('%postBURST');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
filename = 'SetupMatrix1024_15MHz_Takoyaki_BURST';
save(['MatFiles/',filename]);
VSX
return

%StartSaving
display('Start Saving')
Control = repmat(struct('Command','set&Run','Parameters',[]),1,1);
nStartAcquire = evalin('base','nStartAcquire');
Control(1).Parameters = {'Parameters',1,'startEvent',nStartAcquire};
evalin('base','Resource.Parameters.startEvent = nStartAcquire;');
assignin('base','Control', Control);
%StartSaving

%% **** Callback routines used by UIControls (UI)  ****

%BURSTmodeCallback
    TW = evalin('base', 'TW');
    
    if (UIState == 1)
        TW.Parameters = [18, 0.67, 3, 1]; % BURST = 3 half-pulses
    elseif (UIState == 2)
        TW.Parameters = [18, 0.67, 7, 1]; % BURST+ = 7 half-pulses
    end

    assignin('base', 'TW', TW);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Receive'};
    assignin('base','Control', Control);
    return
%BURSTmodeCallback

%preBURST
    P = evalin('base', 'P');
    numFrames_slider = ceil(get(hObject, 'Value'));
    numFrames_box = str2double(get(hObject,'String'));
    if isnan(numFrames_box)
        numFrames = numFrames_slider; 
    else 
        numFrames = numFrames_box; 
    end 

    iV_original = evalin('base', 'P.numPostColFrames + P.numPreColFrames + P.numColFrames;');
    pre_BURST_original = evalin('base', 'P.numPreColFrames');
    iV_curr = iV_original - pre_BURST_original + numFrames;

    assignin('base', 'numFrames_pre', numFrames);
    evalin('base', 'P.numPreColFrames = numFrames_pre;')
    assignin('base', 'iV_curr', iV_curr); 
    evalin('base', 'P.iV = 1:iV_curr;')
    evalin('base', 'Resource.RcvBuffer(1).numFrames = iV_curr;')
    evalin('base', 'Resource.ImageBuffer(1).numFrames = iV_curr;')
    evalin('base', 'P.nRcvFrms = iV_curr;')
    evalin('base', 'P.nImgFrms = iV_curr;')
    
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'SeqControl'};
    assignin('base','Control', Control);
    return
%preBURST

%durBURST
    P = evalin('base', 'P');
    numFrames_slider = ceil(get(hObject, 'Value'));
    numFrames_box = str2double(get(hObject,'String'));
    if isnan(numFrames_box)
        numFrames = numFrames_slider; 
    else 
        numFrames = numFrames_box; 
    end 
    
    iV_original = evalin('base', 'P.numPostColFrames + P.numPreColFrames + P.numColFrames;');
    dur_BURST_original = evalin('base', 'P.numColFrames');
    iV_curr = iV_original - dur_BURST_original + numFrames; 
    
    assignin('base', 'numFrames_post', numFrames);
    evalin('base', 'P.numColFrames = numFrames_post;')
    assignin('base', 'iV_curr', iV_curr); 
    evalin('base', 'P.iV = 1:iV_curr;')
    evalin('base', 'Resource.RcvBuffer(1).numFrames = iV_curr;')
    evalin('base', 'Resource.ImageBuffer(1).numFrames = iV_curr;')
    evalin('base', 'P.nRcvFrms = iV_curr;')
    evalin('base', 'P.nImgFrms = iV_curr;')
    
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'SeqControl'};
    assignin('base','Control', Control);
    return
%durBURST

%postBURST
    P = evalin('base', 'P');
    numFrames_slider = ceil(get(hObject, 'Value'));
    numFrames_box = str2double(get(hObject,'String'));
    if isnan(numFrames_box)
        numFrames = numFrames_slider; 
    else 
        numFrames = numFrames_box; 
    end 

    iV_original = evalin('base', 'P.numPostColFrames + P.numPreColFrames + P.numColFrames;');
    post_BURST_original = evalin('base', 'dur_BURST_original');
    iV_curr = iV_original - post_BURST_original + numFrames;
    
    assignin('base', 'numFrames_post', numFrames);
    evalin('base', 'P.numPostColFrames = numFrames_post;')
    assignin('base', '1:iV_curr', iV_curr); 
    evalin('base', 'P.iV = 1:iV_curr;')
    evalin('base', 'Resource.RcvBuffer(1).numFrames = iV_curr;')
    evalin('base', 'Resource.ImageBuffer(1).numFrames = iV_curr;')
    evalin('base', 'P.nRcvFrms = iV_curr;')
    evalin('base', 'P.nImgFrms = iV_curr;')
    
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'SeqControl'};
    assignin('base','Control', Control);
    return
%postBURST

%wideBeam
    Resource = evalin('base', 'Resource');
    nWideBeam = evalin('base','nWideBeam');
    P = evalin('base', 'P');
    Trans = evalin('base','Trans');
    hv = evalin('base','TPC.hv');

    % Start event sequence at the wide beam index
    Resource.Parameters.startEvent = nWideBeam;

    % Save parameter struct 
%     unique_ID = evalin('base', 'unique_ID');
%     P.BURSTSaveDir = ['BURST_files_' unique_ID '_' char(datetime('now', 'Format', 'MMddhhmmss'))];
%     evalin('base' , "param_filename = [acq_dir script_type '_params_' unique_ID, '_' , char(datetime('now', 'Format', 'MMddhhmmss'))];");
%     evalin('base', 'save(param_filename);');

    assignin('base','P', P);
    assignin('base','Resource',Resource);
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};
    assignin('base','Control', Control);
    return
%wideBeam



function SensCutoffCallback(~,~,UIValue)
%SensCutoff - Sensitivity cutoff change
    ReconL = evalin('base', 'Recon');
    for i = 1:size(ReconL,2)
        ReconL(i).senscutoff = UIValue;
    end
    assignin('base','Recon',ReconL);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Recon'};
    assignin('base','Control', Control);
end

function YSectionChange(~,~,UIValue)
%YSectionChange
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(1).ReferencePt;
    RefPt(2) = UIValue;
    RS.DisplayWindow(1).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    PDataT = evalin('base','PData');
    P = evalin('base','P');

    PDataT(1).Region(1).Shape.oPAIntersect = RefPt(2);

    PDataT(1).Region = computeRegions(PDataT(1));
    assignin('base','PData',PDataT);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',1,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end

function XSectionChange(~,~,UIValue)
%XSectionChange
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(2).ReferencePt;
    RefPt(1) = UIValue;
    RS.DisplayWindow(2).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    PDataT = evalin('base','PData');
    P = evalin('base','P');

    PDataT(1).Region(1).Shape.oPAIntersect = RefPt(1);

    PDataT(1).Region = computeRegions(PDataT(1));
    assignin('base','PData',PDataT);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',2,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end

function ZSectionChange(~,~,UIValue)
%ZSectionChange
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(3).ReferencePt;
    RefPt(3) = UIValue;
    RS.DisplayWindow(3).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    PDataT = evalin('base','PData');
    P = evalin('base','P');

    PDataT(1).Region(1).Shape.oPAIntersect = RefPt(3);

    PDataT(1).Region = computeRegions(PDataT(1));
    assignin('base','PData',PDataT);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',3,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end

function CompressionP(~,~,UIValue)
%CompressionP
    CompressionFactor = UIValue;
    assignin('base','CompressionFactor',CompressionFactor);
end

%% **** Callback routines used by External function definition (EF) ****
function resetStartEvent
Resource = evalin('base', 'Resource');
Resource.Parameters.startEvent = 1;
assignin('base','Resource',Resource);
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',Resource.Parameters.startEvent};
assignin('base','Control', Control);
end

function volumetricPlot(ImgData)

    persistent handle3D

    PData = evalin('base','PData');
    Trans = evalin('base','Trans');
    P = evalin('base','P');
    v_depth = 1;
    CFactor = evalin('base','CompressionFactor');
    ImgData=ImgData.^(1/CFactor);
    ImgData(:,:,1:v_depth)=0;
    ImgData = flipdim(ImgData,3);

    if isempty(handle3D) || ~ishandle(handle3D)
        figure('Position', [430, 50, 450, 450]);
        handle3D = axes;
    end

    set(handle3D,'NextPlot','replacechildren');
    vol3d('cdata', ImgData,'texture', '2D','Parent', handle3D);
    grid(handle3D, 'on'); colormap(handle3D,'gray');
    set(handle3D,'NextPlot','add');

    xl=PData(1).Size(2);
    yl=PData(1).Size(1);
    zl=PData(1).Size(3);
    dx=PData(1).PDelta(2);
    dy=PData(1).PDelta(1);
    dz=PData(1).PDelta(3);
    plot3(handle3D,Trans.ElementPosWL(:,1)./dx+(xl-1)/2,Trans.ElementPosWL(:,2)./dy+(yl-1)/2,(Trans.ElementPosWL(:,3)+P.endDepth)./dz,'k.');
    view(handle3D,-45,30);
end

function saveImageData(ImgData)

   persistent bloc_count;
      if isempty(bloc_count)
         bloc_count = 1;
      end
   P = evalin('base','P');
   Path = P.pathName;

   if bloc_count == 1 
      mkdir(Path)
      save([Path '\P.mat'],'P');
   end
   eval(['pAM3D_',num2str(bloc_count),' = squeeze(ImgData);']);
   save([Path '\pAM3D_' num2str(bloc_count) '.mat'],['pAM3D_',num2str(bloc_count)]);
   display(['Saved image #',num2str(bloc_count) '!'])
   bloc_count = bloc_count+1;
   if bloc_count == P.nImgFrms + 1
       display('End of the acquisition!')
   end
end

function vstic(varargin)
    tic
end

function vstoc(varargin)
    toc
end

function saveRFData(RFData)

   persistent RF_count;
      if isempty(RF_count)
         RF_count = 1;
      end
   P = evalin('base','P');
   Path = P.pathName;
   
   if RF_count == 1 
      mkdir(Path)
      save([Path '\P.mat'],'P');
   end
   
   eval(['RFData_',num2str(RF_count),' = squeeze(RFData);']);
   save([Path '\RFData_' num2str(RF_count) '.mat'],['RFData_',num2str(RF_count)]);
   display(['Saved RFData #',num2str(RF_count) '!'])
   RF_count = RF_count+1;
   if RF_count == P.numberImagesToSave + 1
       display('End of the acquisition!')
   end
end

function setVoltage
P = evalin('base','P');
% Reset the voltage after collapse
hv1Sldr = findobj('Tag','hv1Sldr');
set(hv1Sldr,'Value',P.collapseVoltage);
hv1Value = findobj('Tag','hv1Value');
set(hv1Value,'String',num2str(P.collapseVoltage,'%.1f'));
feval(get(hv1Sldr,'Callback'), hv1Sldr);
end

%EF#4
function resetVoltage
P = evalin('base','P');
% Reset the voltage after collapse
hv1Sldr = findobj('Tag','hv1Sldr');
set(hv1Sldr,'Value',P.imagingVoltage);
hv1Value = findobj('Tag','hv1Value');
set(hv1Value,'String',num2str(P.imagingVoltage,'%.1f'));
feval(get(hv1Sldr,'Callback'), hv1Sldr);
end