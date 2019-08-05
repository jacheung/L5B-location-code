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

fileName = 'glm_simplified';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] = designMatrixBlocks_simplified(selectedArray,glmnetOpt,glmModel);
end
    
%GLMdesign Matrix Build
selectedFeatures = [2:4]; 
interpOption = 'off'; %linear interpolation of missing values;
selectedFeaturesOptions = fields(glmModel{1}.io.components);
selectedFeaturesTitles = selectedFeaturesOptions(selectedFeatures);
[glmModel] = designMatrixBuilder_hilbert(glmModel,glmnetOpt,selectedFeatures,interpOption);


for i = 1:length(glmModel)
 responseTrials = ~glmModel{i}.io.DmatY==0;
 disp(['num response trials = ' num2str(sum(responseTrials))])
 if numel(unique(glmModel{i}.io.DmatY(responseTrials)))>1
 disp(['iterating for neuron ' num2str(i) '/' num2str(length(tunedIdx))])
 glmModel{i} = poissonModel(glmModel{i}.io.DmatXNormalized(responseTrials,:),glmModel{i}.io.DmatY(responseTrials),selectedArray{i},glmnetOpt,glmModel{i});
 glmModel{i}.meta = tunedIdx(i);
 glmModel{i}.name = fileName;
 else
     disp('skipping trial. not enough unique touch responses')
 end
 
end


glmnet(glmModel{1}.io.DmatXNormalized