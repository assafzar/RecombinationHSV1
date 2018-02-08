function [ROI] = EnoshManualAnnotation(I,roiFnamePrefix,always)

Iann = I;
I1 = Iann(:,:,1);
I2 = Iann(:,:,2);
I3 = Iann(:,:,3);

% allVals = [I488norm(:);I561norm(:)];
maxVal = 255;

roiFnameMat = [roiFnamePrefix '.mat'];
roiFnameJpg = [roiFnamePrefix '.jpg'];


if exist(roiFnameMat,'file') && ~always
    load(roiFnameMat); % ROI
    return;
end

if ~exist(roiFnameMat,'file')
    ROI = false(size(Iann,1),size(Iann,2));
else
    load(roiFnameMat); %ROI
    I1(ROI) = maxVal;
    I2(ROI) = maxVal;
    I3(ROI) = maxVal;
    Iann(:,:,1) = I1;
    Iann(:,:,2) = I2;
    Iann(:,:,3) = I3;
end

tmp = Iann;
stop = false;
while (~stop)
    roi = roipoly(tmp);
    if (sum(roi(:)) == 0)
        stop = true;
    else
        ROI = ROI | roi;
    end
    I1(ROI) = maxVal;
    I2(ROI) = maxVal;
    I3(ROI) = maxVal;
    tmp(:,:,1) = I1;
    tmp(:,:,2) = I2;
    tmp(:,:,3) = I3;
end

save(roiFnameMat,'ROI');
imwrite(tmp,roiFnameJpg);

end