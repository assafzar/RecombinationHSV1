% dname = '/project/bioinformatics/Danuser_lab/shared/assaf/OrenKobilerTAU/201702/31'
% addpath(genpath('/project/bioinformatics/Danuser_lab/shared/assaf/OrenKobilerTAU/code'));
function cellMeasures = EnoshColocAnalysis(dname)

close all;

imgDname = [dname filesep 'imgs'];
setFileNames(imgDname); % changing file names to be Matlab friendly

fovRois = {};
cellMeasures.coloc = [];
cellMeasures.union = [];
cellMeasures.intersect = [];
cellMeasures.relCenterAreaCh0 = [];
cellMeasures.relCenterAreaCh1 = [];
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
    segParams.chPixTH = 40;
    segParams.chGaussianSigma = 2;
    
    curFovRois = enoshSegment(dname,prefix,segParams);   
    
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
    
    %%
    curMeasures = calcCellMeasures(curFovRois);
    for icell = 1 : length(curMeasures)
        cellMeasures.coloc = [cellMeasures.coloc curMeasures{icell}.coloc];
        cellMeasures.union = [cellMeasures.union curMeasures{icell}.uniteArea];
        cellMeasures.intersect = [cellMeasures.intersect curMeasures{icell}.intersectArea];
        cellMeasures.relCenterAreaCh0 = [cellMeasures.relCenterAreaCh0 curMeasures{icell}.ch0Area];
        cellMeasures.relCenterAreaCh1 = [cellMeasures.relCenterAreaCh1 curMeasures{icell}.ch1Area];
    end
end
save([dname filesep 'stats.mat'],'cellMeasures');
end

%%
function curMeasures = calcCellMeasures(curFovRois)
for icell = 1 : curFovRois.n
    curCellCh0 = curFovRois.debugCh0{icell};
    curCellCh1 = curFovRois.debugCh1{icell};
    curCellNuc = curFovRois.debugNuc{icell};
    
    
    
    nucArea = sum(curCellNuc(:));
    intersectArea = (curCellCh0 > 0) & (curCellCh1 > 0);
    uniteArea = curCellCh0 | curCellCh1;
    curMeasures{icell}.nucArea = sum(nucArea(:));
    curMeasures{icell}.intersectArea = sum(intersectArea(:));
    curMeasures{icell}.uniteArea = sum(uniteArea(:));
    curMeasures{icell}.coloc = curMeasures{icell}.intersectArea ./ curMeasures{icell}.uniteArea;
    curMeasures{icell}.ch0Area = sum(curCellCh0(:) > 0) ./ sum(curCellNuc(:));
    curMeasures{icell}.ch1Area = sum(curCellCh1(:) > 0) ./ sum(curCellNuc(:));
end
end