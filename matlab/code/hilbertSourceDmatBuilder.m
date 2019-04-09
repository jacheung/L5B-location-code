clearvars -except U

touchCells = touchCell(U);
selectedCells = find(touchCells==1);
%%
%BUILD parameters
buildIndices = [-25:50]; %Indices around touch

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

%basis function and convolution using gaussian distribution
glmnetOpt.bf.bfwidth =7;
glmnetOpt.bf.bfstd = 5;
glmnetOpt.bf.bfspacing = 3;
basisFunction = normalize_var(normpdf(-1*glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth,0,glmnetOpt.bf.bfstd),0,1);
glmnetOpt.bf.indicesToAdd  = [-33:glmnetOpt.bf.bfspacing:20];

%empty indices for outputs. consider building a struct.
fitDevExplained = cell(length(selectedCells),1);
mccDiff = cell(length(selectedCells),1);

fileName = 'glmModelCurvature';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
end

%% BUILDING MODEL USING GLM BINOMIAL
for i = [ocellidx']
    array = U{selectedCells(i)};
    
    %Defining features
    touchmat = zeros(array.t,array.k);
    touchIdx = [find(array.S_ctk(9,:,:)==1) ;find(array.S_ctk(12,:,:)==1)];
    spikes = squeeze(array.R_ntk);
    phase = squeeze(array.S_ctk(5,:,:));
    amplitude = squeeze(array.S_ctk(3,:,:));
    midpoint = squeeze(array.S_ctk(4,:,:));
    angle = squeeze(array.S_ctk(1,:,:));
    curvature = squeeze(array.S_ctk(17,:,:));
    touchmat(touchIdx)=1;
    
    %Get rid of touches with nan values
    nanidx = unique([find(isnan(midpoint(touchIdx))) ; find(isnan(amplitude(touchIdx))) ; find(isnan(phase(touchIdx))) ; find(isnan(curvature(touchIdx))) ] );
    touchIdx(nanidx)=[];
    
    %convoluting features with basis fucntions
    phase_conv = conv2(basisFunction,1,phase,'same');
    amplitude_conv = conv2(basisFunction,1,amplitude,'same');
    midpoint_conv = conv2(basisFunction,1,midpoint,'same');
    curvature_conv = conv2(basisFunction,1,curvature,'same');
    touchmat_conv = conv2(basisFunction,1,touchmat,'same');
    
    %INDEX BUILDER FOR LAGS
    startIdx = buildIndices'+touchIdx';
    
    %LAG AND LAGS APART
    shiftIdx = nan(length(glmnetOpt.bf.indicesToAdd ),numel(startIdx));
    touchShiftIdxRaw=nan(numel(buildIndices),length(glmnetOpt.bf.indicesToAdd ));
    for b = 1:length(glmnetOpt.bf.indicesToAdd )
        tmpIdx{b} = startIdx + glmnetOpt.bf.indicesToAdd (b);
        shiftIdx(b,:) = tmpIdx{b}(:);
        
        touchShiftIdxRaw(:,b) = (buildIndices + glmnetOpt.bf.indicesToAdd (b))==0;
        touchShiftIdxRaw(find(touchShiftIdxRaw(:,b)==1)+(-glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth),b) = basisFunction;
        %         touchShiftIdxRaw(find(touchShiftIdxRaw(:,b)==1)+(0:length(basisFunction)-1),b) = basisFunction;
    end
    
    touchShiftIdx = touchShiftIdxRaw(:,(glmnetOpt.bf.indicesToAdd <=0));
    
    DmatX = [repmat(touchShiftIdx,length(touchIdx),1) phase_conv(shiftIdx') amplitude_conv(shiftIdx') midpoint_conv(shiftIdx') curvature_conv(shiftIdx')];
%         DmatX = [repmat(touchShiftIdx,length(touchIdx),1) phase_conv(shiftIdx') amplitude_conv(shiftIdx') midpoint_conv(shiftIdx')];
    DmatY = spikes(startIdx(:));
    
    %remove any NAN trials
    [row,~] = find(isnan(DmatX));
    trialsToRemove = unique(ceil(row./length(buildIndices)));
    trialsToRemoveIdx = (length(buildIndices)*(trialsToRemove-1)) + (1:(length(buildIndices))) ;
    DmatX(trialsToRemoveIdx,:) = [] ;
    DmatY(trialsToRemoveIdx,:) = [] ;
    
    %STANDARDIZATOIN
    DmatXNorm = (DmatX-nanmean(DmatX))./nanstd(DmatX);
    
    %VIEW DESIGN MAT FOR xNumTrials;
    %         figure(9);clf
    %         numTrials = 5;
    %         startTouch = datasample(1:length(touchIdx)-numTrials,1);
    %         imagesc(DmatXNorm(length(buildIndices)*startTouch+1:(length(buildIndices)*startTouch+1)+length(buildIndices).*numTrials,:))
    %         caxis([0 1])
    %GLMNET
    
    sangles = cell(1,glmnetOpt.numIterations);
    sHeat = cell(1,glmnetOpt.numIterations);
    sFitCoeffs = cell(1,glmnetOpt.numIterations);
    fitDevExplained = nan(1,glmnetOpt.numIterations);
    fitDevExplainedFull = nan(1,glmnetOpt.numIterations);
    mccDiff = nan(glmnetOpt.numIterations,2);
    numTrials = size(DmatX,1)/length(buildIndices);
    
    for p = 1:glmnetOpt.numIterations
        disp(['running glmnet with elasticNet ' num2str(glmnetOpt.alpha) ' for iteration ' num2str(p) '/' num2str(glmnetOpt.numIterations) '. Cell ' num2str(i) '/' num2str(length(selectedCells)) ])
        % test and train sets
        % might want to include stratification to make sure all ranges
        % are used
        exampleIdx = 1:numTrials;
        shuffIdx = exampleIdx(randperm(length(exampleIdx)));
        trainIdxStartRaw = shuffIdx(1:round(numTrials*.7));
        testIdxStartRaw = setdiff(shuffIdx,trainIdxStartRaw);
        
        trainIdxStart = ((trainIdxStartRaw-1).*length(buildIndices))+1;
        testIdxStart = ((testIdxStartRaw-1).*length(buildIndices))+1;
        
        trainIdx = trainIdxStart'+(0:length(buildIndices)-1);
        testIdx = testIdxStart'+(0:length(buildIndices)-1);
        
        trainDmatX = DmatXNorm(trainIdx',:);
        trainDmatY = DmatY(trainIdx',:);
        testDmatX = DmatXNorm(testIdx',:);
        testDmatY = DmatY(testIdx',:);
        
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
        
        m = size(testDmatX,1);
        predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
        predProb = sigmoid(predicts); %this is the sigmoid of predicts which means the probabilities
        
        %This is to check that our deviance calculation is correct. Trained
        %model outputs deviance explained but we need for test model which
        %we show below.
        trainnullLL = sum(log(binopdf(trainDmatY,1,(ones(length(trainDmatY),1)*mean(trainDmatY)))));
        trainpredProb = sigmoid([ones(length(trainDmatX),1) trainDmatX] *fitCoeffs);
        trainfullLL = sum(log(binopdf(trainDmatY,1, trainpredProb)));
        trainDevianceExplained = 1-(-trainfullLL)/(-trainnullLL);
        devExplained = cv.glmnet_fit.dev(iLambda);
        
        % goodness of fit using deviance
        mu = mean(testDmatY);
        
        %goodness of fit using selectedTouch Response windows
        reshapedDmatY = reshape(testDmatY,length(buildIndices),length(testDmatY)./length(buildIndices));
        reshapedPredProb = reshape(predProb,length(buildIndices),length(predProb)./length(buildIndices));
        selIdx = [find(buildIndices==array.meta.responseWindow(1)):find(buildIndices==array.meta.responseWindow(2))];
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
        rawFR = reshape(testDmatY,length(buildIndices),length(testIdxStart));
        predFR = reshape(predProb,length(buildIndices),length(testIdxStart));
        bounds = [-100:2:100];
        rawAngle = angle(touchIdx);
        [sangles{p},sidx] = sort(rawAngle(testIdxStartRaw));
        sHeat{p} = predFR(:,sidx)';
        sRaw{p} = rawFR(:,sidx)';
        sFitCoeffs{p} = fitCoeffs;
        
    end
    
    
    %constructing output structure of model
    glmModel{i}.modelParams = glmnetOpt; 
    
    
    glmModel{i}.basisFunctions.touch = touchShiftIdx;
    glmModel{i}.basisFunctions.features = touchShiftIdxRaw;
    glmModel{i}.coeffs.touchCoeffs = mean(coeffsIter(2:2+size(touchShiftIdx,2)-1,:),2);
    glmModel{i}.coeffs.phaseCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
    glmModel{i}.coeffs.ampCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
    glmModel{i}.coeffs.midpointCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
    glmModel{i}.coeffs.curvatureCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
    
    
    coeffsIter = cell2mat(sFitCoeffs);
    glmModel{i}.coeffs.raw = cell2mat(sFitCoeffs);

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

%% QUANTIFICATION OF MODEL
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\pop.mat')
OLcells = find(~(cellfun(@isempty,glmModel))); 
for i = OLcells(2:end)
    
    array = U{selectedCells(i)};
%PLOTTING BASIS FUNCTIONS
figure(480);clf
subplot(2,1,1)
plot(buildIndices,touchShiftIdx,'k')
subplot(2,1,2)
plot(buildIndices,touchShiftIdxRaw,'k')

%PLOTTING FOR time from touch onset
%     figure(800);subplot(5,6,i)
figure(484);clf
subplot(2,1,1)
%     touchSpks = spikes(repmat(touchIdx,1,length(buildIndices))+repmat(buildIndices,size(touchIdx,1),1));
%     hold on; plot(buildIndices,smooth(mean(touchSpks*1000,1)),'k')
hold on;plot(buildIndices,smooth(mean(glmModel{i}.predicted.spikeTestRaw *1000,1),length(basisFunction)),'k','linewidth',2)
hold on;plot(buildIndices,mean(glmModel{i}.predicted.spikeProb *1000,1),'b','linewidth',2)
ylabel('firing rate (hz)')
title(['mean psth of real vs model. Rsqd = ' num2str(corr(smooth(mean(glmModel{i}.predicted.spikeTestRaw *1000,1)'),mean(glmModel{i}.predicted.spikeProb *1000,1)').^2)]);
set(gca,'xtick',[])

figure(484);subplot(2,1,2);
hold on;
plot(buildIndices,smooth(std(glmModel{i}.predicted.spikeTestRaw )*1000,length(basisFunction)),'k')
plot(buildIndices,std(glmModel{i}.predicted.spikeProb )*1000,'b')
title('variance of real(black) vs model(blue)')
%hold on;bar(startIndices,mean(cell2mat(sHeat')*1000,1),'c')
%ylabel('firing rate (hz)')
%title(['mean real fr = ' num2str(mean(touchSpks(:))*1000) 'Hz'])
xlabel('time from touch onset')

thetabounds = [-100:2:100];

%plotting of original heatmap
%toss touches out that are >10 degree difference from nearest bin.


    gaussFilt = [2];
    
    toKeep = diff(pop{selectedCells(i)}.theta.range)<10;
    if ~isempty(find(toKeep==0))
        toKeep(find(toKeep==0):end) = 0;
    end
    
    figure(5480);clf
    subplot(2,3,1)
    imagesc(imgaussfilt(pop{selectedCells(i)}.theta.spikes(toKeep,:),gaussFilt,'padding','replicate'));
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',1:sum(toKeep),'yticklabel',[pop{selectedCells(i)}.theta.range(toKeep)],...
        'xtick',(0:25:length(buildIndices)),'xticklabel',[min(buildIndices):25:max(buildIndices)],'xlim',[0 length(buildIndices)]);
    for k=1:size(pop{selectedCells(i)}.theta.counts(toKeep),1)
        text(20,k,num2str(pop{selectedCells(i)}.theta.counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([find(buildIndices==0)-1 find(buildIndices==0)-1],[length(pop{selectedCells(i)}.theta.range(toKeep)) 0],'w:')
    ylabel('angle at touch')
    title('original tuning curve')
    
    selIdx = [find(buildIndices==array.meta.responseWindow(1)):find(buildIndices==array.meta.responseWindow(2))];
    hold on; plot([selIdx(1) selIdx(1)],[0 30],'-.w')
    hold on; plot([selIdx(end) selIdx(end)],[0 30],'-.w')

    %plotting of modeled heat map
    [sorted, sortedBy ,binBounds]=binslin(glmModel{i}.predicted.angles,glmModel{i}.predicted.spikeProb,'equalE',numel(thetabounds)+1,thetabounds(1),thetabounds(end));
    
    selBins = ~cellfun(@isempty ,sortedBy);
    spksBinned = cell2mat(cellfun(@(x) mean(x,1),sorted,'uniformoutput',0));
    
    modelBounds = bounds(selBins);
    heatPlot = spksBinned(selBins,:);
    ogBounds = pop{selectedCells(i)}.theta.range(toKeep);
    toPlotBounds = find(ismember(modelBounds,ogBounds));
    normHeatPlot =(heatPlot - min(heatPlot(:))) ./ (max(heatPlot(:))-min(heatPlot(:)));
    
    
    figure(5480);
    subplot(2,3,4)
    imagesc(imgaussfilt(normHeatPlot(toPlotBounds,:),gaussFilt,'padding','replicate'))
    set(gca,'ytick',1:length(toPlotBounds),'yticklabel', modelBounds(toPlotBounds),'ydir','normal',...
        'xtick',(0:25:length(buildIndices)),'xticklabel',[min(buildIndices):25:max(buildIndices)],'xlim',[0 length(buildIndices)])
    hold on;plot([find(buildIndices==0)-1 find(buildIndices==0)-1],[length(modelBounds) 0],'w:')
    caxis([0 prctile(normHeatPlot(:),99)]) %cap limit at 99th percentile in case of outliers in predictions
    xlabel('time from touch onset (ms)')
    ylabel('angle at touch')
    title('glmNet result')
    hold on; plot([selIdx(1) selIdx(1)],[0 30],'-.w')
    hold on; plot([selIdx(end) selIdx(end)],[0 30],'-.w')
    
    
    %PLOTTING LOCATION TUNING CURVE
    subplot(2,3,[2 3 5 6])
    colors = {'k','r'};
    
    spikesGroup = {glmModel{i}.raw.spikes , glmModel{i}.predicted.spikeProb};
    angleGroup = {glmModel{i}.raw.angle, glmModel{i}.predicted.angles};
    
    for d = 1:length(spikesGroup)
        allspks = spikesGroup{d};
        angleatT = angleGroup{d};
        touchResp = U{selectedCells(i)}.meta.responseWindow;
        
        
        bl = allspks(:,1:find(buildIndices==0));
        response = allspks(:,find(buildIndices==touchResp(1)):find(buildIndices==touchResp(2)));
        
        meanbase =  mean(mean(bl,2));
        stdbase = std(mean(bl,2));
        meanresp = mean(response,2);
        zscore = (meanresp - meanbase) ./ stdbase;
        
        [sorted,~,~] = binslin(angleatT,zscore,'equalE',numel(thetabounds),thetabounds(1),thetabounds(end));
        
        samps = cellfun(@numel,sorted);
        selBins = find(samps>sum(samps)./100);
        
        zraw = nan(length(sorted),2000);
        for b = 1:length(sorted)
            if sum(b == selBins)>0
                currvals = sorted{b};
                if ~isempty(sorted{b})
                    zraw(b,1:length(currvals)) = currvals';
                end
            end
        end
        
        [p,~,stats]=anova1(zraw',[],'off');
        vals =multcompare(stats,[],'off');
        popzscore{i} = cellfun(@mean,sorted);
        otune(i)=p;
        
        x=cellfun(@str2num,stats.gnames);
        
        %CI BINS
        cibins = nan(size(x,1),1);
        
        SEMg = nanstd(zraw,[],2) ./ sqrt(sum(~isnan(zraw),2));
        for p = 1:length(x)
            rawx = zraw(x(p),:);
            SEM = SEMg(x(p));
            ts = tinv(.95,sum(~isnan(rawx),2)-1);      % T-Score
            cibins(p,:) = ts.*SEM;   %confidence intervals
        end
        
        y = stats.means;
        x=thetabounds(x)';
        
        xy = [x y'];
        
        hold on; filex_shadedErrorBar(xy(:,1),xy(:,2),cibins,['-' colors{d}]);
        set(gca,'ytick',round(min(xy(:,2))-1):round(max(xy(:,2))+1),'ylim',[round(min(xy(:,2))-1) round(max(xy(:,2))+1)],...
            'xlim',[min(xy(:,1)) max(xy(:,1))],'xtick',-20:20:60)
        
        zresp{d} = xy;
    end
    
    [~,interIdx ] = intersect(zresp{1}(:,1),zresp{2}(:,1));
    [~,interIdx2 ] = intersect(zresp{2}(:,1),zresp{1}(:,1));
    
    
    xlabel('angle at touch');ylabel('z-scored response')
    title(['raw(black) vs model(red) tuning. expVar=' num2str(corr(zresp{1}(interIdx,2),zresp{2}(interIdx2,2)).^2)])

    expVar(i) = corr(zresp{1}(interIdx,2),zresp{2}(interIdx2,2)).^2; 
end

%Plotting feature weights
% coeffsIter = cell2mat(sFitCoeffs);
% touchCoeffs = mean(coeffsIter(2:2+size(touchShiftIdx,2)-1,:),2);
% phaseCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
% ampCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
% midpointCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
% figure(9);clf
% hold on; plot(startIndices,sum(touchShiftIdx.*touchCoeffs',2))
% hold on; plot(startIndices,sum(touchShiftIdxRaw.*phaseCoeffs',2))
% hold on; plot(startIndices,sum(touchShiftIdxRaw.*ampCoeffs',2))
% hold on; plot(startIndices,sum(touchShiftIdxRaw.*midpointCoeffs',2))
% set(gca,'xlim',[-25 50],'xtick',-25:25:50)
% ylabel('kernel weight')
% xlabel('time from touch onset (ms)')
% legend('touch','phase','amp','midpoint')
%

figure(9);clf
hold on; plot(-25:50,sum(glmModel{i}.coeffs.touchCoeffs'.*glmModel{i}.basisFunctions.touch,2))
hold on; plot(-25:50,sum(glmModel{i}.coeffs.phaseCoeffs'.*glmModel{i}.basisFunctions.features,2))
hold on; plot(-25:50,sum(glmModel{i}.coeffs.ampCoeffs'.*glmModel{i}.basisFunctions.features,2))
hold on; plot(-25:50,sum(glmModel{i}.coeffs.midpointCoeffs'.*glmModel{i}.basisFunctions.features,2))
% hold on; plot(-25:50,sum(glmModel{i}.coeffs.curvatureCoeffs'.*glmModel{i}.basisFunctions.features,2))
set(gca,'xlim',[-25 50],'xtick',-25:25:50)
ylabel('kernel weight')
xlabel('time from touch onset (ms)')
legend('touch','phase','amp','midpoint','curvature')
