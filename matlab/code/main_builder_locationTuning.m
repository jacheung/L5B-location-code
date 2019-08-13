%Load whisking and neural time series struct 
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)~=0);
pole_tuned = object_location_quantification(U,selectedCells,'pole'); %for old see object_location_v1.0 

%% population at touch pole decoding
% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 10;
glmnetOpt.numInterpPts = 10;
glmnetOpt.pctSamplingThreshold = .80; 
glmnetOpt.interpResolution = 40; %10mm / numBins (e.g. 40 = .25mm resolution) 
glmnetOpt.samplingOption = 'poisson';
glmnetOpt.numberResamples = 50;

fileName = 'glm_location_decoder';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] =designMatrixBlocks_poleDecoder(glmModel,pole_tuned,glmnetOpt);
end

mdlResults = {};

DmatXraw = cell2mat(cellfun(@(x) x.io.DmatX,glmModel,'uniformoutput',0)); 
DmatXnorm = (DmatXraw - mean(DmatXraw)) ./ std(DmatXraw);
mdlResults.io.Y.normal = glmModel{1}.io.DmatY;

mdlResults = multinomialModel(mdlResults,DmatXnorm,mdlResults.io.Y.normal,glmnetOpt);
%%
resolution_in_X = decoderPerformance(mdlResults);
usedUnits = cellfun(@(x) x.params.cellNum,glmModel);

% sampled_mm = mean(cellfun(@(x) range(x.meta.ranges),U(usedUnits)) ./ 10000);
sampled_mm = 10;
resolution_in_mm = (sampled_mm ./ glmnetOpt.interpResolution) * resolution_in_X;

subplot(1,2,1)
title(['bin resolution = ' num2str(sampled_mm ./ glmnetOpt.interpResolution) 'mm'])
subplot(1,2,2)
title(['decoding resolution = ' num2str(resolution_in_mm) 'mm'])
