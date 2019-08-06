%Load whisking and neural time series struct 
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
% U = defTouchResponse(U,.95,'off');
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)~=0);
is_tuned = object_location_quantification(U,selectedCells,'pole');

%%
tunedIdx = find(is_tuned==1);
selectedArray = U(tunedIdx);

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 5;

%BUILD parameters
glmnetOpt.buildIndices = [-25:50]; %Indices around touch
glmnetOpt.touchDirection = 'protraction';

fileName = 'glm_simplified';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] = designMatrixBlocks_simplified(selectedArray,glmnetOpt,glmModel);
end
    
%GLMdesign Matrix Build
selectedFeatures = [2:6]; 
interpOption = 'off'; %linear interpolation of missing values;
selectedFeaturesOptions = fields(glmModel{1}.io.components);
selectedFeaturesTitles = selectedFeaturesOptions(selectedFeatures);
[glmModel] = designMatrixBuilder_simplified(glmModel,glmnetOpt,selectedFeatures,interpOption);

for i = 1:length(glmModel)
    disp(['iterating for neuron ' num2str(i) '/' num2str(length(tunedIdx))])
    glmModel{i} = poissonModel(glmModel{i}.io.DmatXNormalized,glmModel{i}.io.DmatY,selectedArray{i},glmnetOpt,glmModel{i});
    glmModel{i}.meta = tunedIdx(i);
    glmModel{i}.name = fileName;
end

%%
tStruct = object_location_quantification(U,cellfun(@(x) x.meta,glmModel),'pole'); 

samplesPerBin = 100; 
rc = numSubplots(length(glmModel)); 

figure(32);clf
for i = 1:length(glmModel)
    poles = normalize_var(glmModel{i}.predicted.pole,-1,1);
    rawResponses = glmModel{i}.predicted.spikeTestRaw;
    predResponses = glmModel{i}.predicted.spikeProb;
    
    numBins = round(numel(poles)./samplesPerBin); 
    [raw_sorted,raw_sortedBy] = binslin(poles,rawResponses,'equalN',numBins);
    [pred_sorted,pred_sortedBy] = binslin(poles,predResponses,'equalN',numBins);
    
    figure(32);subplot(rc(1),rc(2),i)
    hold on; plot(cellfun(@median, raw_sortedBy),normalize_var(smooth(cellfun(@mean, raw_sorted)),0,1),'b')
    hold on; plot(cellfun(@median, pred_sortedBy),normalize_var(smooth(cellfun(@mean, pred_sorted)),0,1),'r')
    
    gof_tuning(i) = corr(smooth(cellfun(@mean, raw_sorted)),smooth(cellfun(@mean, pred_sorted)));
    pct_responsive(i) = sum(glmModel{i}.io.DmatY>0) ./ numel(glmModel{i}.io.DmatY); 
    
end

gof_de = cellfun(@(x) mean(x.gof.devExplained),glmModel); 
mod_depth = cellfun(@(x) x.mod_depth,tStruct(cellfun(@(y) y.meta,glmModel))); 

figure(8);clf
subplot(2,2,1)
scatter(mod_depth,gof_de)
title(['corr = ' num2str(corr(mod_depth',gof_de'))])
ylabel('dev explained')
xlabel('modulation depth')

subplot(2,2,2)
scatter(pct_responsive,gof_de)
title(['corr = ' num2str(corr(pct_responsive',gof_de'))])
xlabel('proportion trials touch responsive')

subplot(2,2,3)
scatter(gof_tuning,gof_de)
title(['corr = ' num2str(corr(gof_tuning',gof_de'))])
ylabel('dev explained')
xlabel('tuning goodness of fit')


    
    
