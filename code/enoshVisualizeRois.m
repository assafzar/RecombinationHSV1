%%
function [] = enoshVisualizeRois(cellData,roiVisDname,ifov)
close all;
debugNuc = combineDebugRois(cellData.debugNuc);
debugCh0 = combineDebugRois(cellData.debugCh0);
debugCh1 = combineDebugRois(cellData.debugCh1);

h = figure; 
imshowpair(cellData.Inuc,debugNuc,'montage');
saveas(h,[roiVisDname filesep num2str(ifov) '_nuclei.jpg']);

h = figure; 
imshowpair(cellData.Ich0,debugCh0,'montage');
saveas(h,[roiVisDname filesep num2str(ifov) '_ch0.jpg']);

h = figure; 
imshowpair(cellData.Ich1,debugCh1,'montage');
saveas(h,[roiVisDname filesep num2str(ifov) '_ch1.jpg']);

h = figure; 
imshowpair(debugCh0,debugCh1,'montage');
saveas(h,[roiVisDname filesep num2str(ifov) '_compare.jpg']);

% 2 channels 
I2chs = uint8(zeros(size(debugCh0,1),size(debugCh0,2),3));
I2chs(:,:,1) = (debugCh0 > 0) .* 120;
I2chs(:,:,2) = (debugCh1 > 0) .* 120;
I2chs(:,:,3) = (debugNuc > 0) .* 80;
imwrite(I2chs,[roiVisDname filesep num2str(ifov) '_overlay.tif']);

% for future marking of individual cells
Iann = uint8(zeros(size(debugCh0,1),size(debugCh0,2),3));
Iann(:,:,3) = combinePair(zeros(size(debugNuc)),debugNuc);% cellData.Inuc
Iann(:,:,1) = combinePair(cellData.Ich0,debugCh0);
Iann(:,:,2) = combinePair(cellData.Ich1,debugCh1);

save([roiVisDname filesep num2str(ifov) '_annotation.mat'],'Iann');
imwrite(Iann,[roiVisDname filesep num2str(ifov) '_annotation.tif']);

%%
Iann_ch0 = uint8(zeros(size(debugCh0,1),size(debugCh0,2),3));
Iann_ch0(:,:,3) = combinePair(zeros(size(debugNuc)),debugNuc);% cellData.Inuc
Iann_ch0(:,:,1) = combinePair(cellData.Ich0,zeros(size(debugNuc)));
Iann_ch0(:,:,2) = combinePair(zeros(size(debugNuc)),debugCh0);

save([roiVisDname filesep num2str(ifov) '_annotation_ch0.mat'],'Iann_ch0');
imwrite(Iann_ch0,[roiVisDname filesep num2str(ifov) '_annotation_ch0.tif']);

%%
Iann_ch1 = uint8(zeros(size(debugCh0,1),size(debugCh0,2),3));
Iann_ch1(:,:,3) = combinePair(zeros(size(debugNuc)),debugNuc);% cellData.Inuc
Iann_ch1(:,:,2) = combinePair(cellData.Ich1,zeros(size(debugNuc)));
Iann_ch1(:,:,1) = combinePair(zeros(size(debugNuc)),debugCh1);

save([roiVisDname filesep num2str(ifov) '_annotation_ch1.mat'],'Iann_ch1');
imwrite(Iann_ch1,[roiVisDname filesep num2str(ifov) '_annotation_ch1.tif']);


% h = figure; 
% imshowpair(cellData.otsu4Ch0,cellData.otsu4Ch1,'montage');
% saveas(h,[roiVisDname filesep num2str(ifov) '_compareOtsu4.jpg']);

end

%%
function out = combinePair(I,ROI)
if sum(I(:)) == 0
    tmp = zeros(size(I));
else
    tmp = (((double(I-prctile(I(:),0.01)))./double(prctile(I(:),99.99) - prctile(I(:),0.01))));
end
out = uint8(min(255,max(0,150 * tmp)));
perim = bwperim(ROI,8);
mask = imdilate(perim,strel('square',1));
out(mask) = 255;%min(255,out(mask) + 70);
end

function out = combineDebugRois(rois)
out = zeros(size(rois{1}));
curVal = 0;
for i = 1 : length(rois)
    curRoi = rois{i};
    out(logical(curRoi)) = curVal;
    out = out + curRoi;    
    curVal = curVal + max(curRoi(:));
end
end