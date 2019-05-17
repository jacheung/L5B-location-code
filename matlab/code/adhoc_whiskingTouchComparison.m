%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
%% Top level parameters and definitions
touchWindow = [-25:50]; %window for analyses around touch

touchCells = touchCell(U,'off');
selectedCells = find(touchCells==1);

% Structure for quantifying tuning and evaluating decoding
popV = touchFeatureBinned(U,touchWindow);

% Defining touch response
U = defTouchResponse(U,.99,'off');

%% Quantifying object location tuning
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
tuneStruct = tuningQuantification(U,popV,selectedCells,fieldsList([ 1 3 4 5]),whichTouches,touchWindow,'off');

%% Quantifying whisking tuning
whisking = whisking_general(U,'off');

hilbertWhisking = whisking_hilbert(U,popV,'off');

%% comparison between whisking and touch for touchCells 
fieldsToCompare = fields(tuneStruct.R_ntk.allTouches);
selectedCells = find(touchCells==1); %look at only touch cells because those are the only ones w/ OL tuning
pThresh = 0.01; 
THilbertTuning = cell(1,length(fieldsToCompare));
WHilbertTuning = cell(1,length(fieldsToCompare));

for g = 1:4
    figure(23);clf;
    currTArray = tuneStruct.R_ntk.allTouches.(fieldsToCompare{g});
    currWArray = hilbertWhisking.R_ntk.(fieldsToCompare{g});
    stimulus = hilbertWhisking.S_ctk.(fieldsToCompare{g});
    
    chosenTArrays = currTArray(selectedCells);
    TouchR = cellfun(@(x) nanmean(x,2),chosenTArrays,'uniformoutput',0);
    TouchSEM = cellfun(@(x) nanstd(x,[],2) ./ sqrt(sum(~isnan(x),2)),chosenTArrays,'uniformoutput',0);
    chosenWArrays = currWArray(selectedCells);
    WhiskR = cellfun(@(x) nanmean(x,1),chosenWArrays,'uniformoutput',0);
    WhiskSEM = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),chosenWArrays,'uniformoutput',0);
   
    for i = 1:length(WhiskR)
        [THilbertTuning{g}(i),~,~] = anova1(chosenTArrays{i}',[],'off');
        [WHilbertTuning{g}(i),~,~] = anova1(chosenWArrays{i},[],'off');
        
        subplot(4,8,i)
        if THilbertTuning{g}(i)<pThresh
             shadedErrorBar(stimulus,TouchR{i}*1000,TouchSEM{i}*1000,'b')
        else
             shadedErrorBar(stimulus,TouchR{i}*1000,TouchSEM{i}*1000,'k')
        end
        
        if WHilbertTuning{g}(i)<pThresh
             hold on; shadedErrorBar(stimulus,WhiskR{i}*1000,WhiskSEM{i}*1000,'r')
        else
             hold on; shadedErrorBar(stimulus,WhiskR{i}*1000,WhiskSEM{i}*1000,'k')
        end 
        xBounds = stimulus(~isnan(TouchR{i}*1000));
        set(gca,'xlim',[min(xBounds) max(xBounds)])
        
        if strcmp(fieldsToCompare{g},'phase')
            set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
        end
        
        if WHilbertTuning{g}(i)<pThresh & THilbertTuning{g}(i)<pThresh
            title(num2str(corr((WhiskR{i}*1000)',TouchR{i}*1000,'rows','complete')));
        end
        
    end
    suptitle(fieldsToCompare{g})
    print(['C:\Users\jacheung\Dropbox\LocationCode\Figures\hilbertCode\Hilbert_whiskingTouch\' fieldsToCompare{g}],'-dpng')
% 
%     figure(30);
%     subplot(2,2,g)
%     allMat = [(cell2mat(WhiskR)*1000)' cell2mat(TouchR')*1000];
%     allMat = allMat(~sum(isnan(allMat),2),:); 
%     scatter(allMat(:,1),allMat(:,2),'k')
%     minMax(1) = min(allMat(:));
%     minMax(2) = prctile(allMat(:),99);
%     hold on; plot(minMax,minMax,'-.k')
%     set(gca,'xlim',minMax,'ylim',minMax)
%     axis square
%     title(fieldsToCompare{g})
%     xlabel('whisking FR');ylabel('touch FR')
    
end



Ttmp = cell2mat(THilbertTuning')<.01;
Wtmp = cell2mat(WHilbertTuning')<.01;


%% comparison
naiveVSexpert = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);

OLexpert = numel(intersect(find(tuneStruct.populationQuant.theta.matrix(1,:) == 1),find(naiveVSexpert==1))) ./ nansum(tuneStruct.populationQuant.theta.matrix(1,:) == 1);
OLnaive = numel(intersect(find(tuneStruct.populationQuant.theta.matrix(1,:) == 1),find(naiveVSexpert==0))) ./ nansum(tuneStruct.populationQuant.theta.matrix(1,:) == 1);
figure(100);clf
subplot(2,2,1);pie([OLexpert OLnaive],{'expert','naive'})
title(['OL cells = ' num2str(nansum(tuneStruct.populationQuant.theta.matrix(1,:) == 1))])

wEXCOL = numel(intersect(find(tuneStruct.populationQuant.theta.matrix(1,:) == 1),find(whisking.matrix(1,:) == 1))) ./ numel(find(tuneStruct.populationQuant.theta.matrix(1,:) == 1)) ;
wINHOL = numel(intersect(find(tuneStruct.populationQuant.theta.matrix(1,:) == 1),find(whisking.matrix(1,:) == -1))) ./ numel(find(tuneStruct.populationQuant.theta.matrix(1,:) == 1));
nsOL = numel(intersect(find(tuneStruct.populationQuant.theta.matrix(1,:) == 1),find(whisking.matrix(1,:) == 0))) ./ numel(find(tuneStruct.populationQuant.theta.matrix(1,:) == 1));

subplot(2,2,1);pie([wEXCOL wINHOL nsOL],{'whisking EXC','whisking INH','whisking ns'})
title(['OL cells = ' num2str(nansum(tuneStruct.populationQuant.theta.matrix(1,:) == 1))])

wEXCnanOL = numel(intersect(find(isnan(tuneStruct.populationQuant.theta.matrix(1,:))),find(whisking.matrix(1,:) == 1))) ./ numel(find(isnan(tuneStruct.populationQuant.theta.matrix(1,:)))) ;
wINHnanOL = numel(intersect(find(isnan(tuneStruct.populationQuant.theta.matrix(1,:))),find(whisking.matrix(1,:) == -1))) ./ numel(find(isnan(tuneStruct.populationQuant.theta.matrix(1,:))));
nsnanOL = numel(intersect(find(isnan(tuneStruct.populationQuant.theta.matrix(1,:))),find(whisking.matrix(1,:) == 0))) ./ numel(find(isnan(tuneStruct.populationQuant.theta.matrix(1,:))));

subplot(2,2,2);pie([wEXCnanOL wINHnanOL nsnanOL],{'whisking EXC','whisking INH','whisking ns'})
title(['nonOL cells = ' num2str(numel(find(isnan(tuneStruct.populationQuant.theta.matrix(1,:)))) )])


wEXCexpert = numel(intersect(find(naiveVSexpert==1),find(whisking.matrix(1,:) == 1))) ./ numel(find(whisking.matrix(1,:) == 1));
wEXCnaive = numel(intersect(find(naiveVSexpert==0),find(whisking.matrix(1,:) == 1))) ./ numel(find(whisking.matrix(1,:) == 1));
subplot(2,2,3);pie([wEXCexpert wEXCnaive],{'expert','naive'})
title(['whisking EXC cells = ' num2str(numel(find(whisking.matrix(1,:) == 1)))])


wINHexpert = numel(intersect(find(naiveVSexpert==1),find(whisking.matrix(1,:) == -1))) ./ numel(find(whisking.matrix(1,:) == -1));
wINHnaive = numel(intersect(find(naiveVSexpert==0),find(whisking.matrix(1,:) == -1))) ./ numel(find(whisking.matrix(1,:) == -1));
subplot(2,2,4);pie([wINHexpert wINHnaive],{'expert','naive'})
title(['whisking INH cells = ' num2str(numel(find(whisking.matrix(1,:) == -1)))])


nsexpert = numel(intersect(find(naiveVSexpert==1),find(whisking.matrix(1,:) == 0))) ./ sum(naiveVSexpert==1);
nsnaive = numel(intersect(find(naiveVSexpert==0),find(whisking.matrix(1,:) == 0))) ./ sum(naiveVSexpert==0);

%% neurons and their ability to decode choice 
%using average fr from pOnset to first lick 

glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.05;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 5;
mcc = cell(1,length(U)); 

for i = 1:length(U)
    if strcmp(U{i}.meta.layer,'BVL5b')
        
        disp(['iterating for cell ' num2str(i) '/' num2str(length(U))])
        
        masks = maskBuilder(U{i});
        spks = squeeze(U{i}.R_ntk);
        
        meanPreDFR = (nanmean(spks.* masks.onsettolick) * 1000)';
        spikeTrainPreDecision = (spks.* masks.onsettoMedianlick)';
        selectedTPs = sum((~isnan(spikeTrainPreDecision)),1) == U{i}.k; 
        spikeTrainPreDecision = spikeTrainPreDecision(:,selectedTPs);
        
        DmatX = spikeTrainPreDecision; 
        DmatY = ~(masks.nonlickTrials)';
        
        if sum(DmatY==0) < 5
            disp(['skipping cell ' num2str(i) ' because not enough non-lick trials'])
        else
            for d = 1:glmnetOpt.numIterations
                [trainIdx,testIdx] = stratifiedIndices(DmatY,.7);
                
                trainDmatX = DmatX(cell2mat(trainIdx'));
                testDmatX = DmatX(cell2mat(testIdx'));
                trainDmatY = DmatY(cell2mat(trainIdx'));
                testDmatY = DmatY(cell2mat(testIdx'));
                
                cv = cvglmnet(trainDmatX,trainDmatY,'binomial',glmnetOpt,[],glmnetOpt.xfoldCV);
                
                fitLambda = cv.lambda_1se;
                iLambda = find(cv.lambda == fitLambda);
                fitCoeffs = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
                predicts = cvglmnetPredict(cv,testDmatX,fitLambda); % this outputs testDmatX * fitCoeffs
                predProb = sigmoid(predicts); %this is the sigmoid of predicts which means the probabilities
                
                mcc{i}(d) = mccCalculator(testDmatY,predProb>.5);
            end
        end
    end
    
end

%or using the spike sequence from pole onset to median lick time. 



%% all neuron comparison 
fieldsToCompare = fields(tuneStruct.R_ntk.allTouches);
for g = 1:4
    currArrays = hilbertWhisking.R_ntk.(fieldsToCompare{g});
    for i = 1:length(U)
        [FullWhiskTuning{g}(i),~,~] = anova1(currArrays{i},[],'off');
    end
end

hilbertWhiskTuning = cell2mat(FullWhiskTuning')<pThresh;

OLtuning = [nan(4,sum(~(touchCells==1))) (cell2mat(THilbertTuning')<pThresh)];
[~,idx] = sort(touchCells);
comparisonMat = [touchCells(idx) ;whisking.matrix(:,idx) ; hilbertWhiskTuning(:,idx) ; OLtuning];

%add choice decoding too?