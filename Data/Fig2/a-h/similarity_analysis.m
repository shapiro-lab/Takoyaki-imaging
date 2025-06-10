%% Transformation matrix (this part is already done)
% import registered masks from X and Y phantoms and get the transformation
% matrix between the two masks using Medical Registration Estimator. The
% resulting tform matrices were stored in each volume folder

volume_num = 6;

for j = volume_num
    root = ['./Volume', num2str(j), '/X_axis'];

    load(fullfile(root, 'phantom.mat'));
    phantom_X = phantom_temp;

    root = ['./Volume', num2str(j), '/Y_axis'];

    load(fullfile(root, 'phantom.mat'));
    phantom_Y = phantom_temp;

end


%% MS-SSIM analysis

volume_num = 1:6;

for j = volume_num
    load(fullfile(['./Volume', num2str(j), '/tform.mat'])); % transformation matrix

    root = ['./Volume', num2str(j), '/X_axis'];
    load(fullfile(root, 'TakoyakiAM/pAM3D_1.mat'));
    TakoyakiAM_X = pAM3D_1;
    clear pAM3D_1

    load(fullfile(root, 'SheetpAM/pAM3D_1.mat'));
    pAM3D_sheet_X = pAM3D_1;
    clear pAM3D_1

    load(fullfile(root, 'xAM3D/xAM3D_1.mat'));
    xAM3D_X = xAM3D_1;

    load(fullfile(root, 'phantom.mat')); % registered mask
    phantom_X = phantom_temp;

    root = ['./Volume', num2str(j), '/Y_axis'];
    load(fullfile(root, 'TakoyakiAM/pAM3D_1.mat'));
    TakoyakiAM_Y = pAM3D_1;
    clear pAM3D_1

    load(fullfile(root, 'SheetpAM/pAM3D_1.mat'));
    pAM3D_sheet_Y = pAM3D_1;
    clear pAM3D_1

    load(fullfile(root, 'xAM3D/xAM3D_1.mat'));
    xAM3D_Y = xAM3D_1;
    
    % rotate Y axis images based on transformation matrix
    warped_point = imwarp(imwarp(TakoyakiAM_Y, movingReg.tCoarse, "OutputView", imref3d(size(TakoyakiAM_X))), movingReg.dispField);
    warped_sheet = imwarp(imwarp(pAM3D_sheet_Y(:,19:end-18,:), movingReg.tCoarse, "OutputView", imref3d(size(TakoyakiAM_X))), movingReg.dispField);
    temp_xAM = zeros(size(TakoyakiAM_X));
    temp_xAM(:,25:end-25,:) = xAM3D_Y(19:end-18,:,:);
    warped_xAM = imwarp(imwarp(temp_xAM, movingReg.tCoarse, "OutputView", imref3d(size(TakoyakiAM_X))), movingReg.dispField);
    
    % region of interest
    bbox = regionprops3(phantom_X, "BoundingBox").BoundingBox;
    
    % calculate multi-scale SSIM
    MSssim_point(j) = multissim3(log10(TakoyakiAM_X(round(bbox(2)):round(sum(bbox([2 5]))), 25:end-25,round(bbox(3)):round(sum(bbox([3 6]))))), log10(warped_point(round(bbox(2)):round(sum(bbox([2 5]))), 25:end-25,round(bbox(3)):round(sum(bbox([3 6]))))));
    MSssim_sheet(j) = multissim3(log10(pAM3D_sheet_X(round(bbox(2)):round(sum(bbox([2 5]))), 43:end-43,round(bbox(3)):round(sum(bbox([3 6]))))), log10(warped_sheet(round(bbox(2)):round(sum(bbox([2 5]))), 25:end-25,round(bbox(3)):round(sum(bbox([3 6]))))));
    MSssim_xAM(j) = multissim3(log10(xAM3D_X(18+(round(bbox(2)):round(sum(bbox([2 5])))), :,round(bbox(3)):round(sum(bbox([3 6]))))+1), log10(warped_xAM(round(bbox(2)):round(sum(bbox([2 5]))), 25:end-25,round(bbox(3)):round(sum(bbox([3 6]))))+1));

    % figure; imagesc(warped_xAM(round(bbox(2)):round(sum(bbox([2 5]))), 25:end-25,70))

end

%% Visualize bar graph

figure; b = bar(["Tako AM" "Sheet pAM" "3D xAM"], [mean(MSssim_point) mean(MSssim_sheet) mean(MSssim_xAM)], ...
    'BarWidth', 0.5);

b_X = b.XEndPoints;
hold on; 
line(repmat(b_X',  [1 6]), [MSssim_point; MSssim_sheet; MSssim_xAM], 'Color', 0.7*[1 1 1], 'Marker', '.', 'MarkerSize', 10, 'LineStyle', '--');

ylim([0 0.6]); ylabel('MS-SSIM')