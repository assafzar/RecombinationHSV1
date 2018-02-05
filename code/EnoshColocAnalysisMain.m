% addpath(genpath('/project/bioinformatics/Danuser_lab/shared/assaf/OrenKobilerTAU/code'));
function [] = EnoshColocAnalysisMain()

always = false;

% dnamePrefix = strrep('C:\Assaf\Research\OrenKobilerTAU\090815_test','\',filesep);
% dnameSuffix = {'31','25X35','32X35'};

dnamePrefix = strrep('C:\Assaf\Research\OrenKobilerTAU\110915','\',filesep);
dnameSuffix = {'31','25X35','32X35','26X35'};

flags.localAnalysis = false;
flags.annotateRois = false;
flags.metaAnalysis = true;

cellMeasureFname = [dnamePrefix filesep 'allMeasures.mat'];
% if ~always && exist(cellMeasureFname,'file')
%     load(cellMeasureFname);
% else


if flags.localAnalysis
    %     expCellMeasures = cell(1,length(dnameSuffix));
    for i = 1 : length(dnameSuffix)
        dname = [dnamePrefix filesep dnameSuffix{i} filesep];
        close all;
        
        EnoshLocalAnalysis(dname);
    end
    
    %     save(cellMeasureFname,'expCellMeasures');
    % end
end

if flags.annotateRois
    alwaysAnnotate = false;
    expCellMeasures = cell(1,length(dnameSuffix));
    expRCMeasures = cell(1,length(dnameSuffix));
    for i = 1 : length(dnameSuffix)
        dname = [dnamePrefix filesep dnameSuffix{i} filesep];
        
        % Also calculates the measures
        [expCellMeasures{i},expRCMeasures{i}] = EnoshColocManualROI(dname,alwaysAnnotate);
    end
    save(cellMeasureFname,'expCellMeasures','expRCMeasures');
end

if flags.metaAnalysis
    load(cellMeasureFname); % expCellMeasures,expRCMeasures
    EnoshColocMetaAnalysis(expCellMeasures,expRCMeasures,dnameSuffix,dnamePrefix);
end
end