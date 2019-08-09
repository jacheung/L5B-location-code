%Load whisking and neural time series struct 
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
% U = defTouchResponse(U,.95,'off');
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)~=0);
tStruct = object_location_quantification(U,selectedCells,'pole');
%%
touchCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)); 
trained = find(cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U));
trained = 1:length(U);

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 25;

for i = 1:length(U)
    mask = maskBuilder(U{i});
    spikes = squeeze(U{i}.R_ntk);
    DmatY = U{i}.meta.trialType'; 
    DmatX = nansum(mask.samplingp .* spikes)';
   
    numTrials = size(DmatY,1);
    exampleIdx = 1:numTrials;
    shuffIdx = exampleIdx(randperm(length(exampleIdx)));
    trainIdxStartRaw = shuffIdx(1:round(numTrials*.7));
    testIdxStartRaw = setdiff(shuffIdx,trainIdxStartRaw);
    
    trainDmatX = DmatX(trainIdxStartRaw,:);
    trainDmatY = DmatY(trainIdxStartRaw,:);
    testDmatX = DmatX(testIdxStartRaw,:);
    testDmatY = DmatY(testIdxStartRaw,:);
    
    cv = cvglmnet(trainDmatX,trainDmatY,'binomial',glmnetOpt,[],glmnetOpt.xfoldCV);
    
    fitLambda = cv.lambda_1se;
    iLambda = find(cv.lambda == fitLambda);
    fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
    predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
    predProb = sigmoid(predicts); %this is the sigmoid of predicts which means the probabilities
    
    predProb>0