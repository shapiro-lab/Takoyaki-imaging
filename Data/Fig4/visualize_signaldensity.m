basedir = './TakoyakiAM_live_imaging';

folder_start = 'pAM3D_*';

folderList = dir(fullfile(basedir, folder_start));

for i = 1:numel(folderList)
    
    folderName = folderList(i).name;
    filePath = fullfile(basedir, folderName);

    data = load(filePath);
    pAM3D(:,:,:,i) = data.(erase(folderName, ".mat"));

    alpha = regexp(folderName, '([\d.]+).mat', 'tokens');
    a(i) = str2double(alpha{1}{1});

end

pAM3D(:,:,:,a) = pAM3D;

%%
load('./labels_2.mat')
load('./TakoyakiAM_live_imaging/timelog.mat')
T = Tend-Tend(1);
label1 = labels == 'Label1';
label2 = labels == 'Label2';
label3 = labels == 'Label3';
label1_c = squeeze(sum(pAM3D.*label1, [1 2 3]));
label2_c = squeeze(sum(pAM3D.*label2, [1 2 3]));
label3_c = squeeze(sum(pAM3D.*label3, [1 2 3]));
figure; plot(T, label1_c/sum(label1(:))); hold on; plot(T, label2_c/sum(label2(:))); plot(T, label3_c/sum(label3(:)))
xlim([0 T(end)]); xlabel('Time (s)'); ylabel('GV signal density (A.U.)');
legends({'Ipsi. LV', 'Contra. LV', '3rd ventricle (dorsal)'})
