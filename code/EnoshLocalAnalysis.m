% dname = '/project/bioinformatics/Danuser_lab/shared/assaf/OrenKobilerTAU/201702/31'
% addpath(genpath('/project/bioinformatics/Danuser_lab/shared/assaf/OrenKobilerTAU/code'));
% function cellMeasures = EnoshColocAnalysisNew(dname)
function [] = EnoshLocalAnalysis(dname)

close all;

imgDname = [dname filesep 'imgs'];
setFileNames(imgDname); % changing file names to be Matlab friendly

fovRois = {};
% cellMeasures = [];
% cellMeasures.coloc = [];
% cellMeasures.union = [];
% cellMeasures.intersect = [];
% cellMeasures.relCenterAreaCh0 = [];
% cellMeasures.relCenterAreaCh1 = [];
for ifov = 0 : 100 % assuming <=100 fovs per folder
    prefix = [num2str(ifov) 'C'];
    
    if ~exist([imgDname filesep prefix '0.tif'],'file')
        continue;
    end
    
    %% segmentation
    roiDname = [dname filesep 'rois'];
    if ~exist(roiDname,'dir')
        mkdir(roiDname)
    end
    segParams.nucW = 20; % mask for nuclei segmentation
    segParams.nucPixTH = 900;
    segParams.nucGaussianSigma = 2;
    
    segParams.chW = 3;%5;
    segParams.chPixTH = 25;%40;
    segParams.chGaussianSigma = 2; % smoothing gaussian kernal size
    
    % default = 4
    segParams.nStd = 4;%4;%3;% 4; % # stds above background intensity levels
    % default = 2
    segParams.nStdInNuc = 2;%
    % default = 3
    segParams.nMultTH = 3;%5;%4;%2;%3; % 2 # population for Otsu
    % default = 1
    segParams.nMultTHForeground = 1;%2; %1 % # of populations to consider as forground 
    
    % Dialte ROI before looking for two populations 
    segParams.postDilatePixels = 5;%10;
    
    curFovRois = enoshSegmentNew(dname,prefix,segParams);   
    
    if ~isstruct(curFovRois)
        warning(['skipping ' dname filesep num2str(ifov)]);
        continue;
    end
    
    save([roiDname filesep prefix '_rois.mat'],'curFovRois');
    
    fovRois{length(fovRois)+1} = curFovRois;
    
    %% visualize segmentation
    roiVisDname = [dname filesep 'roisVis'];
    if ~exist(roiVisDname,'dir')
        mkdir(roiVisDname)
    end
    
    enoshVisualizeRois(curFovRois,roiVisDname,ifov);
    
    %% TODO: move to after manual ROI annotation! (and exclude output)
    %     curMeasures = calcCellMeasures(curFovRois);
    %     %     for icell = 1 : length(curMeasures)
    %     %         %         cellMeasures.coloc = [cellMeasures.coloc curMeasures{icell}.coloc];
    %     %         %         cellMeasures.union = [cellMeasures.union curMeasures{icell}.uniteArea];
    %     %         %         cellMeasures.intersect = [cellMeasures.intersect curMeasures{icell}.intersectArea];
    %     %         %         cellMeasures.relCenterAreaCh0 = [cellMeasures.relCenterAreaCh0 curMeasures{icell}.ch0Area];
    %     %         %         cellMeasures.relCenterAreaCh1 = [cellMeasures.relCenterAreaCh1 curMeasures{icell}.ch1Area];
    %     cellMeasures = [cellMeasures,curMeasures];
    %     %     end
end
% save([dname filesep 'stats.mat'],'cellMeasures');
end

%%
% function curMeasures = calcCellMeasures(curFovRois)
% for icell = 1 : curFovRois.n
%     curCellCh0 = curFovRois.debugCh0{icell};
%     curCellCh1 = curFovRois.debugCh1{icell};
%     curCellNuc = curFovRois.debugNuc{icell};        
%     
%     nucArea = sum(curCellNuc(:));
%     intersectArea = (curCellCh0 > 0) & (curCellCh1 > 0);
%     uniteArea = curCellCh0 | curCellCh1;
%     curMeasures{icell}.nucArea = sum(nucArea(:));
%     curMeasures{icell}.intersectArea = sum(intersectArea(:));
%     curMeasures{icell}.uniteArea = sum(uniteArea(:));    
%     curMeasures{icell}.ch0Area = sum(curCellCh0(:) > 0) ./ nucArea;
%     curMeasures{icell}.ch1Area = sum(curCellCh1(:) > 0) ./ nucArea;
%     
%     % Measures per replication center
%     curMeasures{icell}.coloc = calcReplicationCentersColocalization(curCellCh0,curCellCh1);    
% end
% end
% 
% %%
% function curMeasures = calcReplicationCentersColocalization(curCellCh0,curCellCh1)
% 
% [labelsCh0,nCh0] = bwlabel(curCellCh0,8);
% [labelsCh1,nCh1] = bwlabel(curCellCh1,8);
% 
% curMeasures.nCh0 = nCh0;
% curMeasures.nCh1 = nCh1;
% curMeasures.colocMin = []; % minimum of two replication centers
% curMeasures.colocMax = [];
% curMeasures.colocUnite = [];
% for i0 = 1 : nCh0
%     curRepCh0 = labelsCh0 == i0;
%     for i1 = 1 : nCh1
%         curRepCh1 = labelsCh1 == i1;
%         area0 = (curRepCh0 > 0);
%         area1 = (curRepCh1 > 0);
%         intersect = area0 & area1;
%         intersectArea = sum(intersect(:)); 
%         if intersectArea > 0
%             unite = area0 | area1;
%             minArea = intersectArea / min(sum(area0(:)),sum(area1(:)));
%             maxArea = intersectArea / max(sum(area0(:)),sum(area1(:)));
%             uniteArea = intersectArea / sum(unite(:));
%             
%             curMeasures.colocMin = [curMeasures.colocMin minArea];
%             curMeasures.colocMax = [curMeasures.colocMax maxArea];
%             curMeasures.colocUnite = [curMeasures.colocUnite uniteArea];
%         end
%     end
% end
% 
% end
