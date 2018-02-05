function [] = EnoshColocMetaAnalysis(expCellMeasures,expRCMeasures,names,outDname)

clc;

nConds = length(names);
condsStr = names;

% Cell measures
nucArea = cell(1,nConds);
intersectArea = cell(1,nConds);
uniteArea = cell(1,nConds);
ch0Area = cell(1,nConds);
ch1Area = cell(1,nConds);
nRepCenPerCellCh0 = cell(1,nConds);
nRepCenPerCellCh1 = cell(1,nConds);

% Replication center measures
colocMin = cell(1,nConds);
colocMax = cell(1,nConds);
colocUnite = cell(1,nConds);

RC_area = cell(1,nConds); % minimal area of replication center
NIntersectsPerCell = cell(1,nConds); % minimal area of replication center
uniteAreaNucleusAreaRatio = cell(1,nConds);


nCells = nan(1,nConds);
nRepCentersCh0 = nan(1,nConds);
nRepCentersCh1 = nan(1,nConds);
nRepCentersPerCellCh0 = nan(1,nConds);
nRepCentersPerCellCh1 = nan(1,nConds);

for i = 1 : nConds
    
    curExpCellMeasures = expCellMeasures{i};
    curExpRCMeasures = expRCMeasures{i};
    
    nRepCentersCh0(i) = 0;
    nRepCentersCh1(i) = 0;
    
    % Cell measures
    nucArea{i} = [];
    intersectArea{i} = [];
    uniteArea{i} = [];
    ch0Area{i} = [];
    ch1Area{i} = [];
    nRepCenPerCellCh0{i} = [];
    nRepCenPerCellCh1{i} = [];
    
    % Replication center measures
    RC_area{i} = [];
    NIntersectsPerCell{i} = [];
    
    %
    uniteAreaNucleusAreaRatio{i} = [];
    
    colocMin{i} = [];
    colocMax{i} = [];
    colocUnite{i} = [];
    
    matForExcel_intersections = [];
    matForExcel_rcs = [];
    matForExcel_noIntersections = [];
    
    nCells(i) = length(curExpCellMeasures);
    
    assert(nCells(i) == length(curExpRCMeasures));
    
    for icell = 1 : nCells(i)
        % Cell measures
        nucArea{i} = [nucArea{i} curExpCellMeasures{icell}.nucArea];
        intersectArea{i} = [intersectArea{i} curExpCellMeasures{icell}.intersectArea]; % ???
        uniteArea{i} = [uniteArea{i} curExpCellMeasures{icell}.uniteArea];
        ch0Area{i} = [ch0Area{i} curExpCellMeasures{icell}.ch0Area];
        ch1Area{i} = [ch1Area{i} curExpCellMeasures{icell}.ch1Area];
        nRepCenPerCellCh0{i} = [nRepCenPerCellCh0{i} curExpCellMeasures{icell}.coloc.nCh0];
        nRepCenPerCellCh1{i} = [nRepCenPerCellCh1{i} curExpCellMeasures{icell}.coloc.nCh1];
        
        % Replication center measures
        RC_area{i} = [RC_area{i} min(curExpCellMeasures{icell}.coloc.areaCh0,curExpCellMeasures{icell}.coloc.areaCh0)];
        NIntersectsPerCell{i} = [NIntersectsPerCell{i} length(curExpCellMeasures{icell}.coloc.colocMin)];
        uniteAreaNucleusAreaRatio{i} = [uniteAreaNucleusAreaRatio{i} (curExpCellMeasures{icell}.uniteArea/curExpCellMeasures{icell}.nucArea)];
        
        
        colocMin{i} = [colocMin{i} curExpCellMeasures{icell}.coloc.colocMin];
        colocMax{i} = [colocMax{i} curExpCellMeasures{icell}.coloc.colocMax];
        colocUnite{i} = [colocUnite{i} curExpCellMeasures{icell}.coloc.colocUnite];
        
        nRepCentersCh0(i) = nRepCentersCh0(i) + curExpCellMeasures{icell}.coloc.nCh0;
        nRepCentersCh1(i) = nRepCentersCh1(i) + curExpCellMeasures{icell}.coloc.nCh1;
                
        nRCs = length(curExpCellMeasures{icell}.coloc.colocMin);
        for iRC = 1 : nRCs
            matForExcel_intersections = [matForExcel_intersections;...
                [icell,...
                iRC,...
                curExpCellMeasures{icell}.nucArea,...
                curExpCellMeasures{icell}.coloc.areaCh0(iRC),...
                curExpCellMeasures{icell}.coloc.areaCh1(iRC),...
                curExpCellMeasures{icell}.coloc.colocMin(iRC)]];
        end
        
        nRCsCh0 = curExpRCMeasures{icell}.RCs.ch0.nCh;
        nRCsCh1 = curExpRCMeasures{icell}.RCs.ch1.nCh;
        for iRC0 = 1 : nRCsCh0
            matForExcel_rcs = [matForExcel_rcs;...
                [icell,0,iRC0,curExpRCMeasures{icell}.nucArea,curExpRCMeasures{icell}.RCs.ch0.areaRC(iRC0)]];
            
            found = false;
            for iRC = 1 : nRCs
                if curExpRCMeasures{icell}.RCs.ch0.areaRC(iRC0) == curExpCellMeasures{icell}.coloc.areaCh0(iRC)
                    found = true;
                end
                if ~found
                    matForExcel_noIntersections = [matForExcel_noIntersections;...
                        [icell,0,iRC0,curExpRCMeasures{icell}.nucArea,curExpRCMeasures{icell}.RCs.ch0.areaRC(iRC0)]];
                end
            end            
        end
        for iRC1 = 1 : nRCsCh1
            matForExcel_rcs = [matForExcel_rcs;...
                [icell,1,iRC1,curExpRCMeasures{icell}.nucArea,curExpRCMeasures{icell}.RCs.ch1.areaRC(iRC1)]];
            
            found = false;
            for iRC = 1 : nRCs
                if curExpRCMeasures{icell}.RCs.ch1.areaRC(iRC1) == curExpCellMeasures{icell}.coloc.areaCh1(iRC)
                    found = true;
                end
                if ~found
                    matForExcel_noIntersections = [matForExcel_noIntersections;...
                        [icell,1,iRC1,curExpRCMeasures{icell}.nucArea,curExpRCMeasures{icell}.RCs.ch1.areaRC(iRC1)]];
                end
            end
        end
        
    end    
    
    nRepCentersPerCellCh0(i)  = nRepCentersCh0(i) / nCells(i);
    nRepCentersPerCellCh1(i)  = nRepCentersCh1(i) / nCells(i);
    
    % Write excel file
    xlswrite([outDname filesep names{i} '_intersectionStats.xlsx'],matForExcel_intersections);
    xlswrite([outDname filesep names{i} '_rcStats.xlsx'],matForExcel_rcs);
    xlswrite([outDname filesep names{i} '_noIntersectionStats.xlsx'],matForExcel_intersections);
end

%% Statistics: plot distributions
for icond = 1 : nConds
    outFnamePrefix = [outDname filesep names{icond}];
    plotEnoshDistribution(colocMin{icond},0.05:0.1:0.95,'colocalization min',[0,1],0.5,[0,0.3],0.1,[outFnamePrefix '_colocMin.jpg']);
    % Histogram - why do we need it?
    %     plotEnoshHistogram(colocMin{icond},'colocalization min',[0,1],0.5,[0,0.3],0.1,[outFnamePrefix '_colocMin.jpg']);
    %     plotEnoshHistogram(colocMin{icond},'colocalization min',[outFnamePrefix '_colocMinHist.jpg']); % numbers instead of %
    plotEnoshDistribution(colocMax{icond},0.05:0.1:0.95,'colocalization max',[0,1],0.5,[0,0.3],0.1,[outFnamePrefix '_colocMax.jpg']);
    plotEnoshDistribution(colocUnite{icond},0.05:0.1:0.95,'colocalization unite',[0,1],0.5,[0,0.3],0.1,[outFnamePrefix '_colocUnite.jpg']);
    close all;
    
    % Area of RC vs. the intersection score
    plotCorr(colocMin{icond},RC_area{icond},'colocalization','(min) rep center area (pixs)',[outFnamePrefix '_corrAreaIntersection.jpg']);
    % # intersections vs. nucleus area
    plotCorr(NIntersectsPerCell{icond},nucArea{icond},'# intersections','nucleus area',[outFnamePrefix '_nIntersectionVsNucArea.jpg']);
    % # intersections vs. (unite intersection are / nucleus area)
    plotCorr(NIntersectsPerCell{icond},uniteAreaNucleusAreaRatio{icond},'# intersections','norm unite intersect area',[outFnamePrefix '_nIntersectionVsNormIntersectArea.jpg']);
    
    %% log
    %     writeLogFile(colocMin{icond},RC_area{icond},ch0Area{icond},ch1Area{icond},[outFnamePrefix '_log.txt'])
    writeLogFile(colocMin{icond},RC_area{icond},[outFnamePrefix '_log.txt'])
end


%% Statistics: compare each pair of conditions
for i = 1 : nConds
    for j = i+1 : nConds
        % Cell measures
        pNucArea = ranksum(nucArea{i},nucArea{j});
        foldColoc = mean(nucArea{i}) / mean(nucArea{j});
        
        pIntersectArea = ranksum(intersectArea{i},intersectArea{j});
        foldIntersectArea = mean(intersectArea{i}) / mean(intersectArea{j});
        
        pUniteArea = ranksum(uniteArea{i},uniteArea{j});
        foldUniteArea = mean(uniteArea{i}) / mean(uniteArea{j});
        
        pCh0Area = ranksum(ch0Area{i},ch0Area{j});
        foldCh0Area = mean(ch0Area{i}) / mean(ch0Area{j});
        
        pCh1Area = ranksum(ch1Area{i},ch1Area{j});
        foldCh1Area = mean(ch1Area{i}) / mean(ch1Area{j});
        
        pNRepCenPerCellCh0 = ranksum(nRepCenPerCellCh0{i},nRepCenPerCellCh0{j});
        foldNRepCenPerCellCh0 = mean(nRepCenPerCellCh0{i}) / mean(nRepCenPerCellCh0{j});
        
        pNRepCenPerCellCh1 = ranksum(nRepCenPerCellCh1{i},nRepCenPerCellCh1{j});
        foldNRepCenPerCellCh1 = mean(nRepCenPerCellCh1{i}) / mean(nRepCenPerCellCh1{j});
        
        pNIntersectPerCell = ranksum(NIntersectsPerCell{i},NIntersectsPerCell{j});
        foldNIntersectPerCell = mean(NIntersectsPerCell{i}) / mean(NIntersectsPerCell{j});
        
        % Replication center measures
        
        pColocMin = ranksum(colocMin{i},colocMin{j});
        foldColocMin = mean(colocMin{i}) / mean(colocMin{j});
        
        colocMinTH_I = colocMin{i}; colocMinTH_I = colocMinTH_I(colocMinTH_I > 0.5);
        colocMinTH_J = colocMin{j}; colocMinTH_J = colocMinTH_J(colocMinTH_J > 0.5);
        pColocMinTH = ranksum(colocMinTH_I,colocMinTH_J);
        foldColocMinTH = mean(colocMinTH_I) / mean(colocMinTH_J);
        
        pColocMax = ranksum(colocMax{i},colocMax{j});
        foldColocMax = mean(colocMax{i}) / mean(colocMax{j});
        
        pColocUnite = ranksum(colocUnite{i},colocUnite{j});
        foldColocUnite = mean(colocUnite{i}) / mean(colocUnite{j});
        
        disp([names{i} ' (N = ' num2str(nCells(i)) ', n = ' num2str(nRepCentersCh0(i)) ', ' num2str(nRepCentersCh1(i)) ') vs. '...
            names{j} ' (N = ' num2str(nCells(j)) ', n = ' num2str(nRepCentersCh0(j)) ', ' num2str(nRepCentersCh1(j)) '):' ]);
        
        disp('');
        disp(['Colocalization (min): pval ' num2str(pColocMin) ', fold ' num2str(foldColocMin)]);
        disp(['Colocalization (min, intersection > 0.5), n = ' num2str(length(colocMinTH_I)) ', ' num2str(length(colocMinTH_J)) ': pval ' num2str(pColocMinTH) ', fold ' num2str(foldColocMinTH)]);
        disp(['# intersection per nucleus (' num2str(mean(NIntersectsPerCell{i})) ' vs. ' num2str(mean(NIntersectsPerCell{j})) '): pval ' num2str(pNIntersectPerCell) ', fold ' num2str(foldNIntersectPerCell)]);
        %         disp(['Colocalization (max): pval ' num2str(pColocMax) ', fold ' num2str(foldColocMax)]);
        %         disp(['Colocalization (unite): pval ' num2str(pColocUnite) ', fold ' num2str(foldColocUnite)]);
        
        disp(['Nuclei area: pval ' num2str(pNucArea) ', fold ' num2str(foldColoc)]);
        disp(['Intersection area: pval ' num2str(pIntersectArea) ', fold' num2str(foldIntersectArea)]);
        disp(['Unite area: pval ' num2str(pUniteArea) ', fold ' num2str(foldUniteArea)]);
        disp(['Ch0 area: pval ' num2str(pCh0Area) ', fold ' num2str(foldCh0Area)]);
        disp(['Ch1 area: pval ' num2str(pCh1Area) ', fold ' num2str(foldCh1Area)]);
        disp(['Ch0: # rep center per cell: pval ' num2str(pNRepCenPerCellCh0) ', fold ' num2str(foldNRepCenPerCellCh0)]);
        disp(['Ch1: # rep center per cell: pval ' num2str(pNRepCenPerCellCh1) ', fold ' num2str(foldNRepCenPerCellCh1)]);
        
        disp('*****************************');
    end
end

save([outDname filesep 'stats.mat'], 'condsStr','colocMin','nCells');



% colorStr = {'c','r','g'};
% h = figure;
% hold on;
% for i = 1 : nConds
%     plot(expCellMeasures{i}.union,expCellMeasures{i}.intersect,sprintf('o%s',colorStr{i}),'MarkerSize',8,'LineWidth',2);
% end
% xlabel('Union (#pixels)');
% ylabel('Intersection (#pixels)');
% % xlim([0,5000]);
% % ylim([0,3500]);
% xlim([0,3000]);
% ylim([0,2000]);
% legend(names);
% hold off;
% saveas(h,[outDname filesep 'summary.jpg']);

% for i = 1 : nConds
%     corr(expCellMeasures{i}.union',expCellMeasures{i}.intersect')
% end



    function [] = plotEnoshDistribution(data,bins,xlabelStr,xlimVals,xlimStep,ylimVals,ylimStep,outFname)
        [nelements, ~] = hist(data,bins);
        anglesDistribution = nelements ./ sum(nelements);
        h = figure;
        hold on;
        bar(bins,anglesDistribution,'r');
        xlabel(xlabelStr,'FontSize',22);
        ylabel('Fraction','FontSize',22);
        haxes = get(h,'CurrentAxes');
        set(haxes,'XLim',xlimVals);
        set(haxes,'XTick',xlimVals(1):xlimStep:xlimVals(2));
        set(haxes,'XTickLabel',xlimVals(1):xlimStep:xlimVals(2));
        set(haxes,'YLim',ylimVals);
        set(haxes,'YTick',ylimVals(1):ylimStep:ylimVals(2));
        set(haxes,'YTickLabel',ylimVals(1):ylimStep:ylimVals(2));
        set(haxes,'FontSize',22);
        set(h,'Color','none');
        hold off;
        saveas(h,outFname);
        % export_fig(outFname);
    end
end

function [] = plotCorr(colocMin,RC_area,xlabelStr,ylabelStr,outFname)
[rho,pval] = corr(colocMin',RC_area');
h = figure;
hold on;
title(sprintf('rho = %.2f, pval = %.5f',rho,pval));
xlabel(xlabelStr);
ylabel(ylabelStr);
scatter(colocMin,RC_area);
hold off;
saveas(h,outFname);
end

function [] = writeLogFile(colocMin,RC_area,outFname)
logOutFname = [outFname(1:end-4) '_log.txt'];
fileID = fopen(logOutFname,'w');
fprintf(fileID,'Rep. Center size (min)\n');
for i = 1 : length(RC_area)
    fprintf(fileID,sprintf('%d ',RC_area(i)));
end
fprintf(fileID,'\nColocalization (min)\n');
for i = 1 : length(colocMin)
    fprintf(fileID,sprintf('%.2f ',colocMin(i)));
end
% THIS IS WRONG BECUASE THE ch0Area, ch1Area was the normalized area in
% each channel
% fprintf(fileID,'\nCh0 areas\n');
% for i = 1 : length(colocMin)
%     fprintf(fileID,sprintf('%.2f ',ch0Area(i)));
% end
% fprintf(fileID,'\nCh1 areas\n');
% for i = 1 : length(colocMin)
%     fprintf(fileID,sprintf('%.2f ',ch1Area(i)));
% end
fclose(fileID);
end