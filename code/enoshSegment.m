function cellData = enoshSegment(dname,prefix,segParams)

%% segment nuclei in image
nucFname = [dname filesep 'imgs' filesep prefix '2.tif'];
Inuc = imread(nucFname);
InucSmooth = imgaussfilt(Inuc,segParams.nucGaussianSigma);
roiNuc = imbinarize(InucSmooth);

if length(unique(roiNuc(:))) == 1
    th50 = prctile(Inuc(:),50); % background
    level = graythresh(Inuc(Inuc > th50));
    roiNuc = imbinarize(Inuc,level);
    warning([dname filesep prefix ' nuclei seg failed going to 50% threshold']);
    if length(unique(roiNuc(:))) == 1
        th80 = prctile(Inuc(:),80); % background
        level = graythresh(Inuc(Inuc > th80));
        roiNuc = imbinarize(Inuc,level);
        warning([dname filesep prefix ' nuclei seg failed going to 80% threshold']);
        if length(unique(roiNuc(:))) == 1
            warning([dname filesep prefix ' nuclei seg failed (too weak SNR?)']);
            cellData = nan;
            return;
        end
    end
end

roiNuc = cleanROI(roiNuc,segParams.nucW,segParams.nucPixTH);
% figure;
% imshowpair(Inuc, roiNuc, 'montage');

if length(unique(roiNuc(:))) == 1
    warning([dname filesep prefix ' nuclei seg failed (too dense cells?)']);
    cellData = nan;
    return;
end

%% initiate data structure
[nucMask,nNuc] = bwlabel(roiNuc,8);

cellData.n = nNuc;
cellData.Inuc = Inuc;
cellData.Ich0 = imread([dname filesep 'imgs' filesep prefix '0.tif']);
cellData.Ich1 = imread([dname filesep 'imgs' filesep prefix '1.tif']);
cellData.otsu4Ch0 = nan;
cellData.otsu4Ch1 = nan;
cellData.nuc = cell(1,nNuc);
cellData.ch0 = cell(1,nNuc);
cellData.ch1 = cell(1,nNuc);
cellData.debugNuc = cell(1,nNuc);
cellData.debugCh0 = cell(1,nNuc);
cellData.debugCh1 = cell(1,nNuc);


%% segment replication center channels
for ich = 0:1
    imgFname = [dname filesep 'imgs' filesep prefix num2str(ich) '.tif'];
    Ich = imread(imgFname);
    IchSmooth = imgaussfilt(Ich,segParams.nucGaussianSigma);
    IchOtsu4 = otsuMultTH(IchSmooth,4);
    
    if ich == 0
        cellData.otsu4Ch0 = IchOtsu4;         
    else
        cellData.otsu4Ch1 = IchOtsu4;
    end
    
    roiChMask = IchOtsu4 == 3 | IchOtsu4 == 4; % based on populations 3&4
    
    if sum(roiChMask(:)) > 0.4*sum(roiNuc(:))
        roiChMask = IchOtsu4 == 4;
    end        
    
    roi0 = roiChMask;
    roi1 = imopen(roi0,strel('square',segParams.chW)); % exclude small excesses
    roi2 = imclose(roi1,strel('square',segParams.chW)); % unite regions with small holes in-between
    roi3 = imopen(roi2,strel('square',segParams.chW)); %
    roi4 = imfill(roi3,'holes');
        
    rois = {};
    
    for il = 1 : nNuc
        outROIch = zeros(size(roi4));
        
        curNucRoi = nucMask == il; % Nucleus mask
        cellData.debugNuc{il} = curNucRoi;
        
        roi5 = roi4;
        roi5(~curNucRoi) = false;
        
        [Lch,nLch] = bwlabel(roi5,8);
        for ilch = 1 : nLch
            curRoiCh = Lch == ilch;
            if sum(curRoiCh(:)) > segParams.chPixTH && sum(curRoiCh(:)) < 0.6*sum(curNucRoi(:))
                outROIch(curRoiCh) = ilch;
                rois{length(rois)+1} = curRoiCh;
            end
        end
        if ich == 0
            cellData.debugCh0{il} = outROIch;            
            cellData.ch0{il}.nRois = length(rois);
            cellData.ch0{il}.rois = rois;
        else
            cellData.debugCh1{il} = outROIch;            
            cellData.ch1{il}.nRois = length(rois);
            cellData.ch1{il}.rois = rois;
        end
    end    
end
end


function outROI = cleanROI(ROI0,W,pixTH)
ROI1 = imopen(ROI0,strel('square',3));
ROI2 = imclose(ROI1,strel('square',W));
ROI3 = imopen(ROI2,strel('square',W));
ROI4 = imfill(ROI3,'holes');
ROI5 = imclearborder(ROI4);

outROI = false(size(ROI5));
[L,nL] = bwlabel(ROI5,8);
for il = 1 : nL
    curRoi = L == il;
    if sum(curRoi(:)) > pixTH
        outROI = outROI | curRoi;
    end
end
end
