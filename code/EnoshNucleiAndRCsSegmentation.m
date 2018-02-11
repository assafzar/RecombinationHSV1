%% EnoshNucleiAndRCsSegmentation - nuclei and RCs segmentation
function [] = EnoshNucleiAndRCsSegmentation(dname)

close all;

imgDname = [dname filesep 'imgs'];
setFileNames(imgDname); % changing file names to be Matlab friendly

fovRois = {};
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
           
    %% Enosh's parameters' setting Feb. 2018
    % default = 4
    segParams.nStd = 2.5;%4;%3;% 4; % # stds above background intensity levels
    % default = 2
    segParams.nStdInNuc = 2;%
    % default = 3
    segParams.nMultTH = 2;%5;%4;%2;%3; % 2 # population for Otsu
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
end
end
