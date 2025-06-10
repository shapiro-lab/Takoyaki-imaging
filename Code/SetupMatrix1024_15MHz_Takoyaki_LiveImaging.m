%% Takoyaki live imaging script
%  This script scans multiple focal points simultaneously using the 1024 element 15 MHz 
%  2D matrix array probe from Vermon with the UTA 1024-MUX adapter on a 
%  Vantage 256 system. All 4 banks are used for transmission. For receive apertures,
%  one complementary set of 4 sparse random apertures are used. The 4
%  receive apertures are synthesized for each transmit aperture.
%
%  P.nCode = 3, Takoyaki AM / P.nCode = 1, Takoyaki Bmode
%  At the beginning, the sequence shows live imaging without saving images.
%  After clicking "save the images", it starts to store the images.
%
%  To run this script, there are two required files:
%  computeTrans_MatrixProbeIncluded.m
%  saved100CompRndApod.mat
%
% written by Sunho Lee
    
% Parameter settings
clear;

P.pathName = 'D:\Sunho\TakoyakiAM' ; 
P.numberImagesToSave = 2;

P.startDepth_mm = 3.3;
P.endDepth_mm = 10;    
P.Voltage = 1.6;

P.numFrames = 2 ;

P.numSyntheticRcvs = 4; % 4 receives for Fully-Sampled Array

P.txFocus_mm = 7; % focal depth
P.aperpatch = 8;
P.numRays = 32 - P.aperpatch + 1; % no. of Rays (for each dimension)
P.nCode = 3; % 3 = 1 full amplitude  + 2 half amplitude

Display_output = 20; %Display mode only: 20 gives a power of 1.0 for a linear output, 40 gives a power of 0.5 for a square root compression
RcvProfile.AntiAliasCutoff = 15; % 10, 15, 20, 30

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

CompRndApod=CompRndApod(1,:,:); % for rcv aperture - will use only one set

% patch-like aperture
Papod = zeros(P.numRays^2, 5, 1024);

for i = 1:P.numRays
    for j = 1:P.numRays
        
        center1 = floor(P.aperpatch/2)+i-1;
        center2 = floor(P.aperpatch/2)+j-1;

        M = zeros(32); % initialize 32 x 32 matrix 
        M(center1-floor(P.aperpatch/2)+1:center1+ceil(P.aperpatch/2),...
            center2-floor(P.aperpatch/2)+1:center2+ceil(P.aperpatch/2)) = 1;
    
        Papod((i-1)*P.numRays + j, 1, :) = reshape(M, [1, 1024]);
        Papod((i-1)*P.numRays + j,2:5,:) = CompRndApod;

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
PData(1).PDelta = [1 1 1]/2;

PData(1).Size(1) = ceil(((P.numRays+3)*Trans.spacing)/PData(1).PDelta(1)); % +3 compensates for gaps
PData(1).Size(2) = ceil(((P.numRays)*Trans.spacing)/PData(1).PDelta(2));
PData(1).Size(3) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(2), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), P.startDepth];
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
            PData.Region((i-1)*P.numRays + j).Shape.w1 = 1.5*Trans.spacing; % to compensate longer jump caused by gaps
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

%% Specify Media. Use point targets in middle of PData.
Media.MP(1,:) = [0,0,60,1.0];      % single point.
Media.numPoints = size(Media.MP,1);

%% Resources
% Resource.VDAS.el2ChMapDisable = 1;

Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*8^2*P.nCode*P.numSyntheticRcvs ;
Resource.RcvBuffer(1).colsPerFrame = 256;
Resource.RcvBuffer(1).numFrames = P.numFrames;
Resource.ImageBuffer(1).numFrames = P.numFrames;
% Resource.ImageBuffer(2).numFrames = P.numFrames;

Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;

Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.25;
Resource.DisplayWindow(1).Position = [0,580, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),-15,PData(1).Origin(3)];
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
Resource.DisplayWindow(2).ReferencePt = [-15,-PData(1).Origin(2),PData(1).Origin(3)];
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
TW.Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D
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

% generate parabolic delay law for each patch
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

% multiple focused points
for i = 1:8
    for j = 1:8
        
        delay_matrix = reshape(TX(j+25*(i-1)).Delay, [32 32]);
        delay_matrix = repmat(delay_matrix(i:i+7,j:j+7),4,4);
        translated_D = reshape(circshift(delay_matrix, [i-1, j-1]), 1, 1024);

        Ap = ones(32);
        % remove truncated unit arrays that do not fully form focused beams
        if i > 1; Ap([1:i-1, 32-8+i:32], :) = 0; end
        if j > 1; Ap(:, [1:j-1, 32-8+j:32]) = 0; end

        Ap = reshape(Ap, 1, 1024);
        
        TX(k+j+8*(i-1)) = TX(1);
        TX(k+j+8*(i-1)).Delay = translated_D;
        TX(k+j+8*(i-1)).focus = 0.0;
        TX(k+j+8*(i-1)).Origin = [0.0, 0.0, 0.0];
        TX(k+j+8*(i-1)).Apod = Ap;
        
        Paper_real(j+(i-1)*8) = computeMuxAperture(Ap, Trans);
        TX(k+j+8*(i-1)).aperture = Paper_real(j+(i-1)*8);
        
        if P.nCode > 1
            TX(k+(i-1)*8+j+8^2) = TX(k+j+8*(i-1));
            TX(k+(i-1)*8+j+2*8^2) = TX(k+j+8*(i-1));

            % checkerboard masking for AM
            mask = repmat([repmat([1 0], 1, 16) repmat([0 1], 1, 16)], 1, 16);
            TX(k+(i-1)*8+j+8^2).Apod = TX(k+j+8*(i-1)).Apod.*mask;
            TX(k+(i-1)*8+j+2*8^2).Apod = TX(k+j+8*(i-1)).Apod.*(~mask);
        end
    end
end

TX = TX(end-8^2*P.nCode+1:end);

% showTXPD
%% allows to visualize the TX in 3D (similar to the EventAnalysisTool)
% f = figure;
% 
% for i=1:64
%     plot3(Trans.ElementPos(find(TX(i).Apod==1),1),Trans.ElementPos(find(TX(i).Apod==1),2),TX(i).Delay(find(TX(i).Apod==1)),'k*');
%     zlim([0 5]);
%     xlim([min(Trans.ElementPos(:,1)) max(Trans.ElementPos(:,1))]);
%     ylim([min(Trans.ElementPos(:,2)) max(Trans.ElementPos(:,2))]);
%     view(gca,-40,49);grid on;
%     zlabel('Delay (wavelength)'); xlabel('X (mm)'); ylabel('Y (mm)')
% 
%     pause(.5);
% 
% end

%% Specify TGC Waveform structure.
TGC(1).CntrlPts = [100 1010 1023 1023 1023 1023 1023 1023];
TGC(1).rangeMax = 128;
TGC(1).Waveform = computeTGCWaveform(TGC);
maxAcqLength = sqrt(P.endDepth^2 + (32*Trans.spacing)^2 + (35*Trans.spacing)^2 );

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
    'callMediaFunc', 0), 1, P.numFrames*P.numSyntheticRcvs*P.nCode*8^2);

for j = 1:P.numFrames
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
                else % for AM
                    Receive(w+kkk+kk+k).Apod = -squeeze(Papod(i,w+1,:))'; % subtract
                    Receive(w+kkk+kk+k).mode = 1; % accumulate
                end
            end
        end
    end
end

%% Recon
senscutoff = 0.6;

Recon(1) = struct('senscutoff', senscutoff, ...
                   'pdatanum', 1, ...
                   'rcvBufFrame', -1, ...
                   'IntBufDest', [1,1],...
                   'ImgBufDest', [1,-1], ...
                   'RINums', 1:8^2*P.numSyntheticRcvs) ;

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...
    'txnum', 1, ...
    'rcvnum', 1, ...
    'scaleFactor', 1, ...
    'regionnum', 1), 1, 8^2*P.numSyntheticRcvs);


ReconInfo(1).mode = 3;
for w=1:8^2
    k= P.numSyntheticRcvs*(w-1);
    for i = 1:P.numSyntheticRcvs
        ReconInfo(i+k).txnum = w;
        ReconInfo(i+k).rcvnum = i+k;
        ReconInfo(i+k).regionnum = w;
        if i == 1
            ReconInfo(i+k).mode = 3;
        end
        if i == P.numSyntheticRcvs
            ReconInfo(i+k).mode = 5;
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
% Save Image Data
Process(5).classname = 'External';  
Process(5).method = 'saveImageData';
Process(5).Parameters = {'srcbuffer','image',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...    % process the most recent frame.
                         'dstbuffer','none'};     

Process(6).classname = 'External';
Process(6).method = 'vstic';
Process(6).Parameters = {'srcbuffer','none'};

Process(7).classname = 'External';
Process(7).method = 'vstoc';
Process(7).Parameters = {'srcbuffer','none'};

% Save RF Data
Process(8).classname = 'External';     
Process(8).method = 'saveRFData';
Process(8).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...      % process the most recent frame.
                         'dstbuffer','none'};                     

%external function
EF(1).Function = vsv.seq.function.ExFunctionDef('volumetricPlot', @volumetricPlot);
EF(2).Function = vsv.seq.function.ExFunctionDef('saveImageData', @saveImageData);
EF(3).Function = vsv.seq.function.ExFunctionDef('vstic', @vstic);
EF(4).Function = vsv.seq.function.ExFunctionDef('vstoc', @vstoc);
EF(5).Function = vsv.seq.function.ExFunctionDef('saveRFData', @saveRFData);

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 250;
SeqControl(3).command = 'timeToNextAcq';  % time between frames
% SeqControl(3).argument = 20000;  % 20 msec
SeqControl(3).argument = 1000000;  % 2s
SeqControl(3).condition = 'ignore'; %Recon processing time
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'triggerOut';
SeqControl(6).command = 'loopCnt'; % - Set loop count. for looping on blocs
SeqControl(6).argument = P.numberImagesToSave/P.numFrames-1;  %
SeqControl(6).condition = 'counter1';
SeqControl(7).command = 'loopTst';  % loop test
SeqControl(7).argument = [];    % set apres
SeqControl(7).condition = 'counter1';

nsc = 8; % nsc is count of SeqControl objects
%% Specify Event Structure arrays
n = 1;   % start index for Events

for j = 1:Resource.RcvBuffer(1).numFrames

    v=P.numSyntheticRcvs*8^2*P.nCode*(j-1);
    for w=1:8^2
        k= P.numSyntheticRcvs*(w-1);
        for i = 1:P.numSyntheticRcvs
            for m = 1:P.nCode
                Event(n).info = 'Transmit/Receive.';
                Event(n).tx = w + (m-1)*8^2;   
                Event(n).rcv = i+k+v+(m-1)*8^2*P.numSyntheticRcvs;  
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = [2,5]; %
                n = n+1;
            end
        end
    end
    Event(n-1).seqControl = [3,5];

    Event(n).info = 'Transfer To Host';
    Event(n).tx = 0;         
    Event(n).rcv = 0;    
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = [nsc,nsc+1]; % set wait time and transfer data
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc + 1;
        SeqControl(nsc).command = 'waitForTransferComplete';
        SeqControl(nsc).argument = nsc-1;
        nsc = nsc + 1;
    n = n+1;
  
    Event(n).info = ['3D Reconstruct Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process XZ - Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process YZ - Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process XY - Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['3D Process  - Frame' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 4;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Return to Matlab';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 4;
    n = n+1;

    Event(n).info = 'HW/SW Sync';
    Event(n).tx = 0;         
    Event(n).rcv = 0;     
    Event(n).recon = 0;    
    Event(n).process = 0;  
    Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'sync';
        nsc = nsc + 1;
    n = n+1;

end

Event(n).info = 'Jump';
Event(n).tx = 0;         % no TX structure.
Event(n).rcv = 0;        % no Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % call processing function
Event(n).seqControl = 1; %
n = n+1;

% Acquire Image Data
nStartAcquire = n;

Event(n).info = 'tic ';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 6;
    Event(n).seqControl = 0;
    n = n+1;

% set loop count
Event(n).info = 'start counter';
Event(n).tx = 0;   
Event(n).rcv = 0;
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = 6;
n = n+1;
SeqControl(7).argument = n;

for j = 1:Resource.RcvBuffer(1).numFrames
    
    v=P.numSyntheticRcvs*8^2*P.nCode*(j-1);
    for w=1:8^2
        k= P.numSyntheticRcvs*(w-1);
        for i = 1:P.numSyntheticRcvs
            for m = 1:P.nCode
                Event(n).info = 'Transmit/Receive.';
                Event(n).tx = w + (m-1)*8^2;   
                Event(n).rcv = i+k+v+(m-1)*8^2*P.numSyntheticRcvs;  
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = [2,5]; %
                n = n+1;
            end
        end
    end
    Event(n-1).seqControl = [3,5];

    Event(n).info = 'Transfer To Host';
    Event(n).tx = 0;        
    Event(n).rcv = 0;   
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = [nsc, nsc+1]; % set wait time and transfer data [nsc,nsc+1]
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc + 1;
        SeqControl(nsc).command = 'waitForTransferComplete';
        SeqControl(nsc).argument = nsc-1;
        nsc = nsc + 1;
    n = n+1;
    
    
    Event(n).info = ['3D Reconstruct Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process XZ - Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process YZ - Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process XY - Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['3D Process  - Frame' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 4;
    Event(n).seqControl = 0;
    n = n+1;
% 
    Event(n).info = 'Return to Matlab';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 4;
    n = n+1;

    Event(n).info = 'HW/SW Sync';
    Event(n).tx = 0;         
    Event(n).rcv = 0;     
    Event(n).recon = 0;    
    Event(n).process = 0;  
    Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'sync';
        nsc = nsc + 1;
    n = n+1;
  
    Event(n).info = 'save Image Data';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 5;    % process
    Event(n).seqControl = 0;
    n = n+1;
    
%     Event(n).info = 'save RF Data'; 
%     Event(n).tx = 0;         % no transmit
%     Event(n).rcv = 0;        % no rcv
%     Event(n).recon = 0;      % no reconstruction
%     Event(n).process = 8;    % process
%     Event(n).seqControl = 0;
%     n = n+1;

    Event(n).info = 'toc ';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 7;
    Event(n).seqControl = 0;
    n = n+1;
end


% Event(n).info = 'save Image Data';
%     Event(n).tx = 0;         % no transmit
%     Event(n).rcv = 0;        % no rcv
%     Event(n).recon = 0;      % no reconstruction
%     Event(n).process = 5;    % process
%     Event(n).seqControl = 0;
%     n = n+1;
    
Event(n).info = 'Loop back until all images are acquired';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % process
Event(n).seqControl = 7;
n = n+1;

%% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl;

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
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
UI(5).Control = VsSliderControl('LocationCode','UserA2','Label','CompressionP',...
                  'SliderMinMaxVal',[1,11,CompressionFactor],...
                  'SliderStep',[1/40,4/40],'ValueFormat','%1.1f', ...
                  'Callback', @CompressionP );

UI(6).Control = {'UserC1','Style','VsPushButton','Label','Save Image'};
UI(6).Callback = text2cell('%StartSaving');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
filename = 'SetupMatrix1024_15MHz_Takoyaki_LiveImaging';
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
   if bloc_count == P.numberImagesToSave + 1
       display('End of the acquisition!')
   end
end

function vstic(varargin)
    tic
end

function vstoc(varargin)
    
    persistent Tend;
    
    if isempty(Tend)
        Tend = toc;
        disp(['Elapsed time is ', num2str(Tend), ' seconds.'])
    else
        Tend(end+1) = toc;
        disp(['Elapsed time is ', num2str(Tend(end)), ' seconds.'])
    end

    P = evalin('base','P');
    if length(Tend) == P.numberImagesToSave
        Path = P.pathName;
        save([Path '\timelog.mat'],'Tend');
    end        

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