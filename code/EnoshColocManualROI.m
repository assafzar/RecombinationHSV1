%% EnoshColocManualROI
%   select ROI for analysis - what cells to analyze based on both channels
%   output accumulated measures for cells and for RCs
%
%
function [cellMeasures,rcMeasures] = EnoshColocManualROI(dname,always)

close all;

imgDname = [dname filesep 'imgs'];

cellMeasures = [];
rcMeasures = [];

for ifov = 0 : 100 % assuming <=100 fovs per folder
    prefix = [num2str(ifov) 'C'];
    
    if ~exist([imgDname filesep prefix '0.tif'],'file')
        continue;
    end
    
    roiDname = [dname filesep 'rois'];
    
    fprintf(sprintf('\n\nEnoshColocManualROI: %s\n',prefix));
    
    load([roiDname filesep prefix '_rois.mat']); % 'curFovRois'
    
    if curFovRois.n == 0
        warning('%s: n ROIs = 0',[roiDname filesep prefix '_rois.mat']);
        continue;
    end
    
    fprintf(sprintf('EnoshColocManualROI: n nuclei = %d\n',curFovRois.n));
    
    roiVisDname = [dname filesep 'roisVis'];
    load([roiVisDname filesep num2str(ifov) '_annotation_ch0.mat']); % 'Iann_ch0'
    roi_ch0 = EnoshManualAnnotation(Iann_ch0,[roiVisDname filesep num2str(ifov) '_annotatedROI_ch0'],always);
    
    load([roiVisDname filesep num2str(ifov) '_annotation_ch1.mat']); % 'Iann_ch1'
    roi_ch1 = EnoshManualAnnotation(Iann_ch1,[roiVisDname filesep num2str(ifov) '_annotatedROI_ch1'],always);
    
    fprintf('EnoshColocManualROI: done annotation 2 channels\n');
    
    roiIntersect = roi_ch0 & roi_ch1;
    
    %% TODO: move to after manual ROI annotation! (and exclude output)
    curMeasures = calcCellMeasures(curFovRois,roiIntersect);
    curRCMeasures = calcRCMeasures(curFovRois,roiIntersect);
    fprintf('EnoshColocManualROI: assign current measures\n');
    cellMeasures = [cellMeasures,curMeasures];
    rcMeasures = [rcMeasures,curRCMeasures];
    fprintf('EnoshColocManualROI: updated cell measures\n');
end
save([dname filesep 'stats.mat'],'cellMeasures','rcMeasures');
end

%%
function [curRCMeasures]  = calcRCMeasures(curFovRois,roiIntersect)
curRCMeasures = {};

if sum(roiIntersect(:) == 0)
    fprintf('calcCellMeasures: ROI intersection is empty\n');
end

ncells = 0;
for icell = 1 : curFovRois.n
    curCellCh0 = curFovRois.debugCh0{icell} > 0 & roiIntersect;
    curCellCh1 = curFovRois.debugCh1{icell} > 0 & roiIntersect;
    curCellNuc = curFovRois.debugNuc{icell} > 0 & roiIntersect;
    
    if sum(curCellCh0(:)) == 0 && sum(curCellCh1(:)) == 0
        continue;
    end
    
    ncells = ncells + 1;
    
    nucArea = sum(curCellNuc(:));
    intersectArea = (curCellCh0 > 0) & (curCellCh1 > 0);
    uniteArea = curCellCh0 | curCellCh1;
    
    curRCMeasures{ncells}.nucArea = sum(nucArea(:));
    curRCMeasures{ncells}.ch0TotalAreaNorm = sum(curCellCh0(:) > 0) ./ nucArea;
    curRCMeasures{ncells}.ch1TotalAreaNorm = sum(curCellCh1(:) > 0) ./ nucArea;
    
    % Measures per replication center
    curRCMeasures{ncells}.RCs.ch0 = calcRcMeasures(curCellCh0,curCellCh1);
    curRCMeasures{ncells}.RCs.ch1 = calcRcMeasures(curCellCh1,curCellCh0);
end

end

%%
% function ID  = findIntersectionCellID(curMeasures,nucArea)
% ID = 0;
% for icell = 1 : length(curMeasures)
%     if curMeasures{icell}.nucArea == nucArea
%         ID = icell;
%         return;
%     end
% end
% end

%%
function curMeasures = calcRcMeasures(curCellCh,curCellOtherCh)

[labelsRC,nRC] = bwlabel(curCellCh,8);

curMeasures.nCh = nRC;

curMeasures.areaRC = [];
curMeasures.dist = []; % Distance of each RC in one channel to RC in the other channel (0 if interacting)
curMeasures.doIntersect = [];

for i = 1 : nRC
    curRepCh = labelsRC == i;
    area = (curRepCh > 0);
    curMeasures.areaRC = [curMeasures.areaRC sum(area(:))];
    
    if sum(curCellOtherCh(:)) == 0
        curMeasures.dist = [curMeasures.dist inf];
        curMeasures.doIntersect = [curMeasures.doIntersect false];
    else
        distMap = bwdist(curRepCh);
        curMeasures.dist = [curMeasures.dist min(distMap(curCellOtherCh))];
        curMeasures.doIntersect = [curMeasures.doIntersect curMeasures.dist(end) == 0];
        assert(length(curMeasures.dist) == length(curMeasures.doIntersect));
    end
end
end

%%
function [curMeasures] = calcCellMeasures(curFovRois,roiIntersect)

curMeasures = {};


if sum(roiIntersect(:) == 0)
    fprintf('calcCellMeasures: ROI intersection is empty\n');
end

ncells = 0;
for icell = 1 : curFovRois.n
    curCellCh0 = curFovRois.debugCh0{icell} > 0  & roiIntersect;
    curCellCh1 = curFovRois.debugCh1{icell} > 0 & roiIntersect;
    curCellNuc = curFovRois.debugNuc{icell} > 0 & roiIntersect;
    
    if sum(curCellCh0(:)) == 0 && sum(curCellCh1(:)) == 0
        continue;
    end
    
    ncells = ncells + 1;
    
    nucArea = sum(curCellNuc(:));
    intersectArea = (curCellCh0 > 0) & (curCellCh1 > 0);
    uniteArea = curCellCh0 | curCellCh1;
    curMeasures{ncells}.nucArea = sum(nucArea(:));
    curMeasures{ncells}.intersectArea = sum(intersectArea(:));
    curMeasures{ncells}.uniteArea = sum(uniteArea(:));
    curMeasures{ncells}.ch0Area = sum(curCellCh0(:) > 0) ./ nucArea; % normalized
    curMeasures{ncells}.ch1Area = sum(curCellCh1(:) > 0) ./ nucArea;
    
    % Measures per replication center
    curMeasures{ncells}.coloc = calcReplicationCentersColocalization(curCellCh0,curCellCh1);
end
end

%% RC-intersection measures within each cell
function curMeasures = calcReplicationCentersColocalization(curCellCh0,curCellCh1)

[labelsCh0,nCh0] = bwlabel(curCellCh0,8);
[labelsCh1,nCh1] = bwlabel(curCellCh1,8);

curMeasures.nCh0 = nCh0;
curMeasures.nCh1 = nCh1;

curMeasures.areaCh0 = [];
curMeasures.areaCh1 = [];

curMeasures.colocMin = []; % minimum of two replication centers
curMeasures.colocMax = [];
curMeasures.colocUnite = [];

for i0 = 1 : nCh0
    curRepCh0 = labelsCh0 == i0;
    for i1 = 1 : nCh1
        curRepCh1 = labelsCh1 == i1;
        area0 = (curRepCh0 > 0);
        area1 = (curRepCh1 > 0);
        intersect = area0 & area1;
        intersectArea = sum(intersect(:));
        if intersectArea > 0
            unite = area0 | area1;
            minArea = intersectArea / min(sum(area0(:)),sum(area1(:)));
            maxArea = intersectArea / max(sum(area0(:)),sum(area1(:)));
            uniteArea = intersectArea / sum(unite(:));
            
            curMeasures.areaCh0 = [curMeasures.areaCh0 sum(area0(:))];
            curMeasures.areaCh1 = [curMeasures.areaCh1 sum(area1(:))];
            
            curMeasures.colocMin = [curMeasures.colocMin minArea];
            curMeasures.colocMax = [curMeasures.colocMax maxArea];
            curMeasures.colocUnite = [curMeasures.colocUnite uniteArea];
        end
    end
end

end
