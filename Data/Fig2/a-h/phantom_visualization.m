%% import data
% This visualization code is optimized for volume2/Y_axis images
TakoyakiAM = load('./volume2/Y_axis/TakoyakiAM/pAM3D_1.mat').pAM3D_1;
SheetpAM = load('./volume2/Y_axis/SheetpAM/pAM3D_1.mat').pAM3D_1;
xAM3D = load('./volume2/Y_axis/xAM3D/xAM3D_1.mat').xAM3D_1;

%% YZ slice
% Takoyaki AM
f = figure; f.Position(3) = f.Position(4)*5;
tiledlayout('horizontal'); nexttile
imagesc(1540/15.625e3/2*(-85:85), 1540/15.625e3/2*(67.5:202.5), 20*log10(squeeze(TakoyakiAM(:,107,:))'/25/mean(TakoyakiAM(:,107,110:end),'all'))); 
colormap gray; clim([-32 0]); daspect([1 1 1]); title('Takoyaki AM')
xlabel('Y (mm)'); ylabel('Z (mm)')

% Sheet-pAM
nexttile; imagesc(1540/15.625e3/2*(-85:85), 1540/15.625e3/2*(67.5:202.5), 20*log10(squeeze(SheetpAM(:,125,:))'/25/mean(SheetpAM(:,125,110:end),'all'))); 
colormap gray; clim([-32 0]); daspect([1 1 1]); title('Sheet pAM')
xlabel('Y (mm)'); ylabel('Z (mm)')

% 3D xAM
nexttile; imagesc(1540/15.625e3/2*(-103:103), 1540/15.625e3/2*(67.5:202.5), 20*log10(squeeze(xAM3D(:,82,:))'/25/mean(xAM3D(:,82,110:end),'all'))); 
colormap gray; c = colorbar; c.Label.String = 'Intensity (dB)'; clim([-32 0]); daspect([1 1 1]); title('3D xAM')
xlabel('Y (mm)'); ylabel('Z (mm)')

%% XZ slice
% Takoyaki AM
f = figure; f.Position(3) = f.Position(4)*4.5;
tiledlayout('horizontal'); nexttile
imagesc(1540/15.625e3/2*(-76:76), 1540/15.625e3/2*(67.5:202.5), 20*log10(squeeze(TakoyakiAM(66,:,:))'/30/mean(TakoyakiAM(66,:,110:end),'all'))); 
colormap gray; clim([-32 0]); daspect([1 1 1]); title('Takoyaki AM')
xlabel('X (mm)'); ylabel('Z (mm)')

% Sheet-pAM
nexttile; imagesc(1540/15.625e3/2*(-94:94), 1540/15.625e3/2*(67.5:202.5), 20*log10(squeeze(SheetpAM(66,:,:))'/30/mean(SheetpAM(66,:,110:end),'all'))); 
colormap gray; clim([-32 0]); daspect([1 1 1]); title('Sheet pAM')
xlabel('X (mm)'); ylabel('Z (mm)')

% 3D xAM
nexttile; imagesc(1540/15.625e3/2*(-51.5:51.5), 1540/15.625e3/2*(67.5:202.5), 20*log10(squeeze(xAM3D(84,:,:))'/30/mean(xAM3D(84,:,110:end),'all'))); 
colormap gray; c = colorbar; c.Label.String = 'Intensity (dB)'; clim([-32 0]); daspect([1 1 1]); title('3D xAM')
xlabel('X (mm)'); ylabel('Z (mm)')

%% XY slice
% Takoyaki AM
f = figure; f.Position(3) = f.Position(4)*4;
tiledlayout('horizontal'); nexttile
imagesc(1540/15.625e3/2*(-76:76), 1540/15.625e3/2*(-85:85), 20*log10(mean(TakoyakiAM(:,:,78:82),3)/10/mean(TakoyakiAM(157:end,:,78:82),'all'))); 
colormap gray; clim([-22 0]); daspect([1 1 1]); title('Takoyaki AM')
xlabel('X (mm)'); ylabel('Y (mm)')

% Sheet-pAM
nexttile; imagesc(1540/15.625e3/2*(-94:94), 1540/15.625e3/2*(-85:85), 20*log10(mean(SheetpAM(:,:,78:82),3)/10/mean(SheetpAM(157:end,:,78:82),'all'))); 
colormap gray; clim([-22 0]); daspect([1 1 1]); title('Sheet pAM')
xlabel('X (mm)'); ylabel('Y (mm)')

% 3D xAM
nexttile; imagesc(1540/15.625e3/2*(-51.5:51.5), 1540/15.625e3/2*(-103:103), 20*log10(mean(xAM3D(:,:,78:82),3)/10/mean(xAM3D(175:end,:,78:82),'all'))); 
colormap gray; c = colorbar; c.Label.String = 'Intensity (dB)'; clim([-22 0]); daspect([1 1 1]); title('3D xAM')
xlabel('X (mm)'); ylabel('Y (mm)')