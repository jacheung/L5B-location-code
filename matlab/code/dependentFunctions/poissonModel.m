function mdl = poissonModel(DmatX,DmatY,selectedArray,glmnetOpt,mdl)

numTrials = size(DmatY,1);
fitDevExplained = nan(1,glmnetOpt.numIterations);

for p = 1:glmnetOpt.numIterations
    disp(['running glmnet with elasticNet ' num2str(glmnetOpt.alpha) ' for iteration ' num2str(p) '/' num2str(glmnetOpt.numIterations)])
    % test and train sets
    % might want to include stratification to make sure all ranges
    % are used
    exampleIdx = 1:numTrials;
    shuffIdx = exampleIdx(randperm(length(exampleIdx)));
    trainIdxStartRaw = shuffIdx(1:round(numTrials*.7));
    testIdxStartRaw = setdiff(shuffIdx,trainIdxStartRaw);

    trainDmatX = DmatX(trainIdxStartRaw,:);
    trainDmatY = DmatY(trainIdxStartRaw,:);
    testDmatX = DmatX(testIdxStartRaw,:);
    testDmatY = DmatY(testIdxStartRaw,:);
    
    %Check that design matrix is properly indexed when we split to test
    %and train set
    %                 figure(19);clf
    %                 numTrialsToPlot = 3;
    %                 startTouch = datasample(1:length(trainIdxStart)-numTrialsToPlot,1);
    %                 imagesc(trainDmatX(length(glmnetOpt.buildIndices)*startTouch+1:(length(glmnetOpt.buildIndices)*startTouch+1)+length(glmnetOpt.buildIndices).*numTrialsToPlot,:))
    %                 set(gca,'xtick',[],'ytick',[])
    %
    cv = cvglmnet(trainDmatX,trainDmatY,'poisson',glmnetOpt,[],glmnetOpt.xfoldCV);
%     cvglmnetPlot(cv)
    
    fitLambda = cv.lambda_1se;
    iLambda = find(cv.lambda == fitLambda);
    fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
    predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
    
    
   model = exp([ones(length(testDmatX),1),testDmatX]*fitCoeffs);
   mu = mean(testDmatY); % null poisson parameter
   nullLogLikelihood = sum(log(poisspdf(testDmatY,mu)));
   saturatedLogLikelihood = sum(log(poisspdf(testDmatY,testDmatY)));
   fullLogLikelihood = sum(log(poisspdf(testDmatY,model)));
   fitDevExplained(p) = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
   devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
    
       
%     %variables for recreating heat map
    rawFR = testDmatY;
    predFR = model;

    sHeat{p} = predFR;
    sRaw{p} = testDmatY;
    sFitCoeffs{p} = fitCoeffs;
end


%constructing output structure of model
mdl.modelParams = glmnetOpt;

mdl.coeffs.raw = cell2mat(sFitCoeffs);


mdl.predicted.spikeTestRaw = cell2mat(sRaw');
mdl.predicted.spikeProb = cell2mat(sHeat');

mdl.gof.devExplained = fitDevExplained;


