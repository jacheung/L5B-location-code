%% Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% neurons and their ability to decode choice 
%using average fr from pOnset to first lick or using the spike train from pole onset to median lick time. 

glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.05;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 5;
mcc = cell(1,2); 

for i = 1:length(U)
    
    if strcmp(U{i}.meta.layer,'BVL5b')        
        disp(['iterating for cell ' num2str(i) '/' num2str(length(U))])
        
        masks = maskBuilder(U{i});
        spks = squeeze(U{i}.R_ntk);
        
        meanPreDFR = (nanmean(spks.* masks.onsettolick) * 1000)';
        spikeTrainPreDecision = (spks.* masks.onsettoMedianlick)';
        selectedTPs = sum((~isnan(spikeTrainPreDecision)),1) == U{i}.k; 
        spikeTrainPreDecision = spikeTrainPreDecision(:,selectedTPs);
        
        DmatX = {meanPreDFR,spikeTrainPreDecision}; 
        DmatY = ~(masks.nonlickTrials)';
        
        if sum(DmatY==0) < 5
            disp(['skipping cell ' num2str(i) ' because not enough non-lick trials'])
        else
            for g = 1:length(DmatX)
                for d = 1:glmnetOpt.numIterations
                    [trainIdx,testIdx] = stratifiedIndices(DmatY,.7);
                    
                    trainDmatX = DmatX{g}(cell2mat(trainIdx'),:);
                    testDmatX = DmatX{g}(cell2mat(testIdx'),:);
                    trainDmatY = DmatY(cell2mat(trainIdx'),:);
                    testDmatY = DmatY(cell2mat(testIdx'),:);
                    
                    cv = cvglmnet(trainDmatX,trainDmatY,'binomial',glmnetOpt,[],glmnetOpt.xfoldCV);
                    
                    fitLambda = cv.lambda_1se;
                    iLambda = find(cv.lambda == fitLambda);
                    fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
                    predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
                    predProb = sigmoid(predicts); %this is the sigmoid of predicts which means the probabilities
                    
                    mcc{g}{i}(d) = mccCalculator(testDmatY,predProb>.5);
                end
            end
        end
    end
    
end

meanMCC = cellfun(@(x) cellfun(@mean,x) , mcc,'uniformoutput',0);

choicePred.rowNames = {'mean FR pre decision','spike train pre decision','mouse accuracy'}; 
choicePred.mat = [cell2mat(meanMCC') ; cellfun(@(x) mean(x.meta.trialCorrect),U)];

%% using the mean fr or spike train from pole onset to average pole down to decode pole location

glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.05;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 5;
testRsqd = cell(1,2); 

for i = 1:length(U)
    disp(['iterating for cell ' num2str(i) '/' num2str(length(U))])
    
    masks = maskBuilder(U{i});
    spks = squeeze(U{i}.R_ntk);
    
    meanFR = (nanmean(spks.* masks.availtoMeanOffset) * 1000)';
    spikeTrain = (spks.* masks.availtoMeanOffset)';
    selectedTPs = sum((~isnan(spikeTrain)),1) == U{i}.k;
    spikeTrain = spikeTrain(:,selectedTPs);
    
    DmatX = {meanFR,spikeTrain};
    DmatY = U{i}.meta.motorPosition';
    
    figure(23);subplot(5,11,i)
    scatter(DmatX{1},DmatY,'.k'); 
    
    for g = 1:length(DmatX)
        for d = 1:glmnetOpt.numIterations
            [trainIdx,testIdx] = randomizedIndices(DmatY,.7);
            
            trainDmatX = DmatX{g}(trainIdx',:);
            testDmatX = DmatX{g}(testIdx',:);
            trainDmatY = DmatY(trainIdx',:);
            testDmatY = DmatY(testIdx',:);
            
            cv = cvglmnet(trainDmatX,trainDmatY,'gaussian',glmnetOpt,[],glmnetOpt.xfoldCV);
            
            fitLambda = cv.lambda_1se;
            iLambda = find(cv.lambda == fitLambda);
            fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
            predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
            
            SSR = sum((predicts-mean(DmatY)).^2);
            SSE = sum((testDmatY-predicts).^2);
            SSTO = SSR + SSE;
            testRsqd{g}{i}(d) = SSR./SSTO;
        end
    end
    
end

meanRsqd = cellfun(@(x) cellfun(@mean,x) , testRsqd,'uniformoutput',0);

tmp = cell2mat(meanRsqd');