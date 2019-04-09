function binomialModel_hilbert(DmatX,DmatY,selectedArray,glmnetOpt)

for i = 1:length(DmatX)
    
    currDmatX = DmatX{i};
    currDmatY = DmatY{i};
    
    %GLMNET 
    sangles = cell(1,glmnetOpt.numIterations);
    sHeat = cell(1,glmnetOpt.numIterations);
    sFitCoeffs = cell(1,glmnetOpt.numIterations);
    fitDevExplained = nan(1,glmnetOpt.numIterations);
    fitDevExplainedFull = nan(1,glmnetOpt.numIterations);
    mccDiff = nan(glmnetOpt.numIterations,2);
    numTrials = size(currDmatY,1)/length(buildIndices);
    
    for p = 1:glmnetOpt.numIterations
        disp(['running glmnet with elasticNet ' num2str(glmnetOpt.alpha) ' for iteration ' num2str(p) '/' num2str(glmnetOpt.numIterations) '. Cell ' num2str(i) '/' num2str(length(selectedCells)) ])
        % test and train sets
        % might want to include stratification to make sure all ranges
        % are used
        exampleIdx = 1:numTrials;
        shuffIdx = exampleIdx(randperm(length(exampleIdx)));
        trainIdxStartRaw = shuffIdx(1:round(numTrials*.7));
        testIdxStartRaw = setdiff(shuffIdx,trainIdxStartRaw);
        
        trainIdxStart = ((trainIdxStartRaw-1).*length(glmnetOpt.buildIndices))+1;
        testIdxStart = ((testIdxStartRaw-1).*length(glmnetOpt.buildIndices))+1;
        
        trainIdx = trainIdxStart'+(0:length(glmnetOpt.buildIndices)-1);
        testIdx = testIdxStart'+(0:length(glmnetOpt.buildIndices)-1);
        
        trainDmatX = currDmatX(trainIdx',:);
        trainDmatY = currDmatY(trainIdx',:);
        testDmatX = currDmatX(testIdx',:);
        testDmatY = currDmatY(testIdx',:);
        
        %Check that design matrix is properly indexed when ew split to test
        %and train set
        %             figure(19);clf
        %             numTrials = 5;
        %             startTouch = datasample(1:length(trainIdxStart)-numTrials,1);
        %             imagesc(trainDmatX(76*startTouch+1:(76*startTouch+1)+76.*numTrials,:))
        
        cv = cvglmnet(trainDmatX,trainDmatY,'binomial',glmnetOpt,[],glmnetOpt.xfoldCV);
        
        fitLambda = cv.lambda_1se;
        iLambda = find(cv.lambda == cv.lambda_1se);
        fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
        predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
        predProb = sigmoid(predicts); %this is the sigmoid of predicts which means the probabilities
        
        %This is to check that our deviance calculation is correct. Trained
        %model outputs deviance explained but we need for test model which
        %we show below.
%         trainnullLL = sum(log(binopdf(trainDmatY,1,(ones(length(trainDmatY),1)*mean(trainDmatY)))));
%         trainpredProb = sigmoid([ones(length(trainDmatX),1) trainDmatX] *fitCoeffs);
%         trainfullLL = sum(log(binopdf(trainDmatY,1, trainpredProb)));
%         trainDevianceExplained = 1-(-trainfullLL)/(-trainnullLL);
%         devExplained = cv.glmnet_fit.dev(iLambda);
        
        % goodness of fit using deviance explained in selected touch
        % response windows 
        mu = mean(testDmatY);
        reshapedDmatY = reshape(testDmatY,length(glmnetOpt.buildIndices),length(testDmatY)./length(glmnetOpt.buildIndices));
        reshapedPredProb = reshape(predProb,length(glmnetOpt.buildIndices),length(predProb)./length(glmnetOpt.buildIndices));
        selIdx = [find(glmnetOpt.buildIndices==selectedArray{i}.meta.responseWindow(1)):find(glmnetOpt.buildIndices==selectedArray{i}.meta.responseWindow(2))];
        selectedDmatY = reshapedDmatY(selIdx,:);
        selectedPredProb = reshapedPredProb(selIdx,:);
        
        fullLL = sum(log(binopdf(selectedDmatY(:),1,selectedPredProb(:))));
        nullLL = sum(log(binopdf(selectedDmatY(:),1,(ones(length(selectedPredProb(:)),1)*mu))));
        saturatedLL = sum(log(binopdf(selectedDmatY(:) ,1,selectedDmatY(:) )));
        
        fitDevExplained(p) = 1 - (saturatedLL-fullLL)/(saturatedLL-nullLL);
        
        %goodness of fit using all indices
        fullLL = sum(log(binopdf(testDmatY,1,predProb)));
        nullLL = sum(log(binopdf(testDmatY,1,(ones(length(predProb),1)*mu))));
        saturatedLL = sum(log(binopdf(testDmatY,1,testDmatY)));
        fitDevExplainedFull(p) = 1 - (saturatedLL-fullLL)/(saturatedLL-nullLL)';
        
        % goodness of fit using MCC
        mu = mean(testDmatY);
        pcBins = 0:.05:1;
        spksProbNull = (ones(length(predProb),1)*mu) > pcBins;
        spksProbFull = predProb>pcBins;
        for j = 1:length(0:.05:1)
            mccFull(j) = mccCalculator(testDmatY,double(spksProbFull(:,j)));
            mccNull(j) = mccCalculator(testDmatY,double(spksProbNull(:,j)));
        end
        [mccDiff(p,1), maxidx] = max(mccFull-mccNull);
        mccDiff(p,2) = pcBins(maxidx);
        
        %variables for recreating heat map
        rawFR = reshape(testDmatY,length(glmnetOpt.buildIndices),length(testIdxStart));
        predFR = reshape(predProb,length(glmnetOpt.buildIndices),length(testIdxStart));
        bounds = [-100:2:100];
        rawAngle = angle(touchIdx);
        [sangles{p},sidx] = sort(rawAngle(testIdxStartRaw));
        sHeat{p} = predFR(:,sidx)';
        sRaw{p} = rawFR(:,sidx)';
        sFitCoeffs{p} = fitCoeffs;
        
    end
    
    
    %constructing output structure of model
    glmModel{i}.modelParams = glmnetOpt; 
    
    coeffsIter = cell2mat(sFitCoeffs);
    glmModel{i}.basisFunctions.touch = touchShiftIdx;
    glmModel{i}.basisFunctions.features = touchShiftIdxRaw;
    
    glmModel{i}.coeffs.raw = cell2mat(sFitCoeffs);
    glmModel{i}.coeffs.touchCoeffs = mean(coeffsIter(2:2+size(touchShiftIdx,2)-1,:),2);
    glmModel{i}.coeffs.phaseCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
    glmModel{i}.coeffs.ampCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
    glmModel{i}.coeffs.midpointCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
        glmModel{i}.coeffs.curvatureCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
    
    glmModel{i}.raw.spikes = reshape(DmatY,length(buildIndices),length(DmatY)/length(buildIndices))';
    glmModel{i}.raw.angle = rawAngle;
    glmModel{i}.raw.angle(trialsToRemove) = [];
    
    
    glmModel{i}.predicted.spikeTestRaw = cell2mat(sRaw');
    glmModel{i}.predicted.spikeProb = cell2mat(sHeat');
    glmModel{i}.predicted.angles = cell2mat(sangles');
    
    glmModel{i}.gof.devExplained = fitDevExplained;
    glmModel{i}.gof.mcc = mccDiff;
    
  
    %predictedHeat
    [sorted, sortedBy ,binBounds]=binslin(cell2mat(sangles'),cell2mat(sHeat'),'equalE',numel(bounds)+1,bounds(1),bounds(end));
    [sortedRaw, sortedByRaw ,binBoundsRaw]=binslin(cell2mat(sangles'),cell2mat(sRaw'),'equalE',numel(bounds)+1,bounds(1),bounds(end));
    selBins = ~cellfun(@isempty ,sortedBy);
    selBinsRaw = ~cellfun(@isempty ,sortedByRaw);
    
    spksBinned = cell2mat(cellfun(@(x) mean(x,1),sorted,'uniformoutput',0));
    spksBinnedRaw = cell2mat(cellfun(@(x) mean(x,1),sortedRaw,'uniformoutput',0));
    
    heatPlot = spksBinned(selBins,:);
    heatPlotRaw = spksBinnedRaw(selBins,:);
    
    modelBounds = bounds(selBins);
    
    normHeatPlot =(heatPlot - min(heatPlot(:))) ./ (max(heatPlot(:))-min(heatPlot(:)));
    normRawPlot = (heatPlotRaw - min(heatPlotRaw(:))) ./ (max(heatPlotRaw(:))-min(heatPlotRaw(:)));
    
    glmModel{i}.heat.angle = modelBounds;
    glmModel{i}.heat.matrix = normHeatPlot;
    glmModel{i}.heat.matrixRaw = normRawPlot;
    
    cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs')
    save(fileName,'glmModel')
end