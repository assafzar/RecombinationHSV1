% addpath(genpath('/project/bioinformatics/Danuser_lab/shared/assaf/OrenKobilerTAU/code'));
function [] = EnoshPoolData()

clc; close all;


mainDname = strrep('C:\Assaf\Research\OrenKobilerTAU\','\',filesep);
conditions = {'31','25X35','32X35','26X35'};
nConditions = length(conditions);
experiments = {'110915'};
nExperiments = length(experiments);
% load([outDname filesep 'stats.mat']); 

colocMinPooled = cell(1,nConditions);
nCellsPooled = zeros(1,nConditions);

for iExp = 1 : nExperiments    
    curFname = [mainDname filesep experiments{iExp} filesep 'stats.mat'];
    if exist(curFname,'file')
        load(curFname); % 'condsStr','colocMin',nCells
        for iCond = 1 : nConditions
            curCond = conditions{iCond};
            condInd = find(strcmp(condsStr,curCond));
            colocMinPooled{iCond} = [colocMinPooled{iCond} colocMin{condInd}];
            nCellsPooled(iCond) = nCellsPooled(iCond)+ nCells(condInd);
        end
    end
end


for i = 1 : nConditions
    for j = i+1 : nConditions
        try
            pColocMin = ranksum(colocMinPooled{i},colocMinPooled{j});
        catch
            fprintf('colocMinPooled\n');
            colocMinPooled
            fprintf(sprintf('i = %d\n',i));
            fprintf(sprintf('j = %d\n',j));
            disp(size(colocMinPooled{i}));
            disp(size(colocMinPooled{j}));
            fprintf('colocMinPooled{i}\n');
            colocMinPooled{i}
            fprintf('colocMinPooled{j}\n');
            colocMinPooled{j}
        end
        foldColocMin = mean(colocMinPooled{i}) / mean(colocMinPooled{j});
        
        colocMinTH_I = colocMinPooled{i}; colocMinTH_I = colocMinTH_I(colocMinTH_I > 0.5);
        colocMinTH_J = colocMinPooled{j}; colocMinTH_J = colocMinTH_J(colocMinTH_J > 0.5);
        pColocMinTH = ranksum(colocMinTH_I,colocMinTH_J);
        foldColocMinTH = mean(colocMinTH_I) / mean(colocMinTH_J);
        
        fprintf([conditions{i} ' (N = ' num2str(nCellsPooled(i)) ') vs. '...
            conditions{j} ' (N = ' num2str(nCellsPooled(j)) ')\n']);                            
        fprintf(['Colocalization (min): pval ' num2str(pColocMin) ', fold ' num2str(foldColocMin) '\n']);
        fprintf(['Colocalization (min, intersection > 0.5), n = ' num2str(length(colocMinTH_I)) ', ' num2str(length(colocMinTH_J)) ': pval ' num2str(pColocMinTH) ', fold ' num2str(foldColocMinTH) '\n']);
        fprintf('\n');
    end
end

end