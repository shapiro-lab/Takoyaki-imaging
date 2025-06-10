volume_num = 1:6;

for j = volume_num
    root = ['./Volume', num2str(j), '/X_axis'];
    load(fullfile(root, 'TakoyakiAM/pAM3D_1.mat'));
    TakoyakiAM = pAM3D_1;
    clear pAM3D_1

    load(fullfile(root, 'SheetpAM/pAM3D_1.mat'));
    pAM3D_sheet = pAM3D_1;
    clear pAM3D_1

    load(fullfile(root, 'xAM3D/xAM3D_1.mat'));
    load(fullfile(root, 'phantom.mat')); % registered mask
    % volshow(TakoyakiAM, OverlayData = phantom_temp); pause

    vol_CBR_point_X(j) = 20*log10(abs(get_volCBR(TakoyakiAM, phantom_temp, 2)));
    vol_CBR_sheet_X(j) = 20*log10(abs(get_volCBR(pAM3D_sheet, phantom_temp, 2)));
    vol_CBR_xAM_X(j) = 20*log10(abs(get_volCBR(xAM3D_1, phantom_temp, 2)));
    clear phanton_temp

    root = ['./Volume', num2str(j), '/Y_axis'];
    load(fullfile(root, 'TakoyakiAM/pAM3D_1.mat'));
    TakoyakiAM = pAM3D_1;
    clear pAM3D_1s

    load(fullfile(root, 'SheetpAM/pAM3D_1.mat'));
    pAM3D_sheet = pAM3D_1;
    clear pAM3D_1

    load(fullfile(root, 'xAM3D/xAM3D_1.mat'));
    load(fullfile(root, 'phantom.mat')); % registered mask
    % volshow(TakoyakiAM, OverlayData = phantom_temp); pause

    vol_CBR_point_Y(j) = 20*log10(abs(get_volCBR(TakoyakiAM, phantom_temp, 1)));
    vol_CBR_sheet_Y(j) = 20*log10(abs(get_volCBR(pAM3D_sheet, phantom_temp, 1)));
    vol_CBR_xAM_Y(j) = 20*log10(abs(get_volCBR(xAM3D_1, phantom_temp, 1)));
    clear phanton_temp

end

%% Plot bar graph
figure; b = bar(["Tako AM" "Sheet pAM" "3D xAM"],[mean(vol_CBR_point_X) mean(vol_CBR_point_Y) ; mean(vol_CBR_sheet_X) mean(vol_CBR_sheet_Y); mean(vol_CBR_xAM_X) mean(vol_CBR_xAM_Y)]);

b_X = [b(1).XEndPoints; b(2).XEndPoints];
hold on; 
line(repmat(b_X(:,1), [1 6]), [vol_CBR_point_X; vol_CBR_point_Y], 'Color', 0.7*[1 1 1], 'Marker', '.', 'MarkerSize', 10, 'LineStyle', '--');
line(repmat(b_X(:,2), [1 6]), [vol_CBR_sheet_X; vol_CBR_sheet_Y], 'Color', 0.7*[1 1 1], 'Marker', '.', 'MarkerSize', 10, 'LineStyle', '--');
line(repmat(b_X(:,3), [1 6]), [vol_CBR_xAM_X; vol_CBR_xAM_Y], 'Color', 0.7*[1 1 1], 'Marker', '.', 'MarkerSize', 10, 'LineStyle', '--');
ylabel('CBR (dB)'); legend({'along X axis', 'along Y axis'})

%%

function vol_CBR = get_volCBR(input_image, bmask, dir)
    
    % change the size of mask to that of input image
    crop_mask = zeros(size(input_image));
    s1 = size(input_image);
    s2 = size(bmask);
    crop_mask((s1(1)-s2(1))/2+1:s1(1)-(s1(1)-s2(1))/2, max(0, (s1(2)-s2(2))/2)+1:s1(2)-max((s1(2)-s2(2))/2, 0), :) = ...
        bmask(:, max(0, (s2(2)-s1(2))/2)+1:s2(2)-max((s2(2)-s1(2))/2, 0), :);

    % divide to left and right
    bmask1 = crop_mask;
    bmask2 = crop_mask;
     se = strel('sphere', 5);
    
    if dir == 2
        bmask1(1:ceil(size(crop_mask,1)/2),:,:) = 0;
        bmask2(ceil(size(crop_mask,1)/2)+1:end,:,:) = 0;
        
    elseif dir == 1
        bmask1(:,1:ceil(size(crop_mask,2)/2),:) = 0;
        bmask2(:,ceil(size(crop_mask,2)/2)+1:end,:) = 0;
    end
    eroded = imerode(crop_mask, se);
    % volshow(input_image, OverlayData = eroded); pause

    well_mean = mean(input_image(logical(eroded)))
    noise_roi = zeros(size(input_image));

    % figure;  
    if dir == 1 % Y axis
        for j = 1:size(input_image, dir) % slice by slice

            frame = squeeze(input_image(j,:,:));
            eroded_mask = squeeze(eroded(j,:,:));
            
            well_roi1 = squeeze(bmask1(j,:,:));
            well_roi2 = squeeze(bmask2(j,:,:));
            % imagesc(well_roi1); pause

            s1 = regionprops(well_roi1, 'Centroid', 'BoundingBox');
            s2 = regionprops(well_roi2, 'Centroid', 'BoundingBox');
            
            if length(s1) == 1 & length(s2) == 1 & sum(eroded_mask(:)) > 1400 % when both wells' mask are full circles
    
                [r1, c1] = find(well_roi1);
                [r2, c2] = find(well_roi2);

                noise_roi(j, max(r2)+7:min(r1)-7, max(min(c1), min(c2))+7:min(max(c1), max(c2))) = 1; % interspace between two wells
                s3 = regionprops(squeeze(noise_roi(j,:,:)), 'BoundingBox');

                if length(s3) == 1 % for visualizaiton
                    % j
                    % imagesc(frame); hold on; rectangle('Position', s1.BoundingBox); rectangle('Position', s2.BoundingBox); rectangle('Position', s3.BoundingBox); pause
                end    
            end
        end

    elseif dir == 2 % X axis

        % apply the mask (truncated due to small FOV) from 3DxAM to the others
        if s1(2)>105 
            eroded(:, [1:(s1(2)+1)/2-52 (s1(2)+1)/2+52:end], :) = 0;
        end

        for j = 1:size(input_image, dir) % slice by slice

            frame = squeeze(input_image(:,j,:));
            well_roi1 = squeeze(bmask1(:,j,:));
            well_roi2 = squeeze(bmask2(:,j,:));
            s1 = regionprops(well_roi1, 'Centroid', 'BoundingBox');
            s2 = regionprops(well_roi2, 'Centroid', 'BoundingBox');
            eroded_mask = squeeze(eroded(:,j,:));
            % imagesc(eroded_mask); pause
            % imagesc(well_roi1); pause 

            if length(s1) == 1 & length(s2) == 1 & sum(eroded_mask(:)) > 1400 % when both wells' mask are full circles
    
                [r1, c1] = find(well_roi1);
                [r2, c2] = find(well_roi2);

                noise_roi(max(r2)+7:min(r1)-7, j, max(min(c1), min(c2))+7:min(max(c1), max(c2))) = 1; % interspace between two wells

                s3 = regionprops(squeeze(noise_roi(:,j,:)), 'BoundingBox');
                if length(s3) == 1 % for visualizaiton
                    % j
                    % imagesc(frame); hold on; rectangle('Position', s1.BoundingBox); rectangle('Position', s2.BoundingBox); rectangle('Position', s3.BoundingBox); 
                    % figure; imagesc(eroded_mask); 
                    % pause;
                end

            end
        end
    end
    
    noise_in_bbox = input_image(logical(noise_roi));
    noise_mean = mean(noise_in_bbox(:))

    vol_CBR = (well_mean - noise_mean)/std(noise_in_bbox(:));

end