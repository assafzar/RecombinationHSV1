function cellData = enoshSegment(dname,prefix,segParams)

% segment nuclei in image
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

[L,nL] = bwlabel(roiNuc,8);

cellData.n = nL;
cellData.Inuc = Inuc;
cellData.Ich0 = imread([dname filesep 'imgs' filesep prefix '0.tif']);
cellData.Ich1 = imread([dname filesep 'imgs' filesep prefix '1.tif']);
cellData.nuc = cell(1,nL);
cellData.ch0 = cell(1,nL);
cellData.ch1 = cell(1,nL);
cellData.debugNuc = cell(1,nL);
cellData.debugCh0 = cell(1,nL);
cellData.debugCh1 = cell(1,nL);

for il = 1 : nL
    curNucRoi = L == il; % only in Nucleus
    cellData.debugNuc{il} = curNucRoi;
    for ich = 0:1
        imgFname = [dname filesep 'imgs' filesep prefix num2str(ich) '.tif'];
        Ich = imread(imgFname);
        IchSmooth = imgaussfilt(Ich,segParams.nucGaussianSigma);
        IchInNuc = IchSmooth(curNucRoi);
        otsuTH = graythresh(IchInNuc);
        roi0 = imbinarize(IchSmooth,otsuTH);
        roi0(~curNucRoi) = false;
        
        roi1 = imopen(roi0,strel('square',3));
        roi2 = imclose(roi1,strel('square',segParams.chW));
        roi3 = imopen(roi2,strel('square',segParams.chW));
        roi4 = imfill(roi3,'holes');
        
        outROIch = zeros(size(roi4));
        rois = {};
        
        [Lch,nLch] = bwlabel(roi4,8);
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
