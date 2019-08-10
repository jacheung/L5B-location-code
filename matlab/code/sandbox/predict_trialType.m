%Load whisking and neural time series struct
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions
% U = defTouchResponse(U,.95,'off');
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)~=0);
tStruct = object_location_quantification(U,selectedCells,'pole');
%%
touchCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U));
nontouchCells = find(~cellfun(@(x) isfield(x.meta,'responseWindow'),U));
trained = find(cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U));
trained = 1:length(U);

selectedUnits = nontouchCells; 


%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

glmModel = [];
for i = 1:length(selectedUnits)
    cellNum = selectedUnits(i); 
    mask = maskBuilder(U{cellNum});
    spikes = squeeze(U{cellNum}.R_ntk);
    DmatY = U{cellNum}.meta.trialType';
%     DmatX = nansum(mask.avail .* mask.touch .* spikes)'; %touches in pole avail period
    DmatX = nansum(mask.avail .* spikes)'; %pole avail
    
    [unique_values,num] = unique(DmatX);
    
    if length(unique_values)>3
        predProb = cell(1,glmnetOpt.numIterations);
        true = cell(1,glmnetOpt.numIterations);
        for p = 1:glmnetOpt.numIterations  
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
            true{p} = testDmatY;
            predProb{p} = sigmoid(predicts); %this is the sigmoid of predicts which means the probabilities
        end
        
        glmModel{i}.cellNumber = cellNum; 
        glmModel{i}.test.true = true;
        glmModel{i}.test.predProb = predProb;
    else
        disp(['skipping cell ' num2str(i) ' b/c too few unique responses'])
    end
    
end

%%
builtUnits = cellfun(@(x) isfield(x,'test'),glmModel);
[tpr,far] = cellfun(@(x) roc(cell2mat(x.test.true')',cell2mat(x.test.predProb')'),glmModel(builtUnits),'uniformoutput',0);

auc = 1-cellfun(@(x,y) trapz(x,y),tpr,far);

figure(8);clf
for g = 1:length(tpr)
    if auc(g)>.5
        hold on; plot(far{g},tpr{g},'g')
    else
        hold on; plot(far{g},tpr{g},'k')
    end   
end
title(['predictive units = ' num2str(sum(auc>.5)) '/' num2str(numel(auc))])
axis square
xlabel('false positive rate');
ylabel('true positive rate'); 
set(gca,'xtick',0:.5:1,'ytick',0:.5:1)



