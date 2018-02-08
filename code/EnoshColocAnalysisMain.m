%% Analysis of a specific experiment with multiple conditions
% Do not forget to add the code to the Matlab path addpath(genpath('....\RecombinationHSV1\code'));
%
%   Input:
%       dnamePrefix - experiment folder
%       dnameSuffix - different experimental conditions (to be compared)
%       Parameters - segParams in EnoshLocalAnalysis
%
%   Processing:
%       1. EnoshNucleiAndRCsSegmentation - segmentation of nuclei, RCs + visualization
%       2. EnoshColocManualROI - (manually) select ROI for analysis and
%       calculate accumulated measures per cell and per RC
%       3. EnoshColocMetaAnalysis - accumulates per experimental condition, output excel files per condition 
%
% Assaf Zaritsky 2017-18 for Enosh Tomer's (Oren Kobiler lab, TAU) project 
function [] = EnoshColocAnalysisMain()

always = false;

% dnamePrefix = strrep('C:\Assaf\Research\OrenKobilerTAU\090815_test','\',filesep);
% dnameSuffix = {'31','25X35','32X35'};

dnamePrefix = strrep('C:\Assaf\Research\OrenKobilerTAU\testData','\',filesep);
dnameSuffix = {'25X35','31'};

flags.segmentation = true;
flags.coloc = true;
flags.metaAnalysis = true;

cellMeasureFname = [dnamePrefix filesep 'allMeasures.mat'];
% if ~always && exist(cellMeasureFname,'file')
%     load(cellMeasureFname);
% else


if flags.segmentation    
    for i = 1 : length(dnameSuffix)
        dname = [dnamePrefix filesep dnameSuffix{i} filesep];
        close all;
        
        EnoshNucleiAndRCsSegmentation(dname);
    end    
end

if flags.coloc
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