%Load whisking and neural time series struct 
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
viewWindow = [-25:50]; %window for analyses around touch

% touchCells = touchCell(U,'on');
% selectedCells = find(touchCells==1);

% Structure for quantifying tuning and evaluating decoding 
popV = touchFeatureBinned(U,viewWindow);

% Defining touch response
U = defTouchResponse(U,.95,'on');
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)==1);
% selectedCells(17) = [];
%% Plotter for feature tuning around touch window
gaussFilt = 1; %smoothing function for tuning plots
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
touchFeatureBinned_plotter(U,popV,selectedCells,fieldsList(1),whichTouches,viewWindow,gaussFilt)

%% Quantifying object location tuning
fieldsList = fields(popV{1}.allTouches);
tunedCellsIdx = tuningQuantification(U,popV,selectedCells,fieldsList(1),whichTouches(1),viewWindow,'on');

%optional raster of FRs for tuned cells. 
for d = 13
    figure(8);clf
    allSpks = squeeze(U{tunedCellsIdx{1}(d)}.R_ntk);
    [~,idx] = sort(U{tunedCellsIdx{1}(d)}.meta.motorPosition);
    allSpks = allSpks(:,idx);
    for k = 1:size(allSpks,2)
        st = find(allSpks(:,k)==1);
        if ~isempty(st)
        figure(8);hold on
        scatter(st,ones(length(st),1).*k,[],'.k')
        end
    end
end

%% Builder for decoding tuning location 
% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 10;

numInterpPts = 10;

%design matrix for feature in fields list
fieldsList = fields(popV{1}.allTouches);
whichTouches = fields(popV{1});

% selectedCells = 1:length(U)
selectedCells = find(touchCells==1);
layersID = cellfun(@(x) x.meta.layer(1),U);
diffID = strfind(layersID,'NB');
expCells = selectedCells(selectedCells>diffID);
naiveCells = selectedCells(selectedCells<=diffID);

[mdl.io.X, mdl.io.Y.normal] = designMatrixBuilder_touchFeature(U,selectedCells,whichTouches(1),viewWindow,numInterpPts);

%multinomial model for decoding location 
clear mdlnull
DmatX = (mdl.io.X.outside_pole - nanmean(mdl.io.X.outside_pole)) ./ nanstd(mdl.io.X.outside_pole); %standardization
mdlnull = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt); 

DmatX = (mdl.io.X.pole_avail - nanmean(mdl.io.X.pole_avail)) ./ nanstd(mdl.io.X.pole_avail); %standardization
mdlavail = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt);

DmatX = (mdl.io.X.touch_responseR - nanmean(mdl.io.X.touch_responseR)) ./ nanstd(mdl.io.X.touch_responseR); %standardization
mdlpertouch = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt);

% DmatX = (mdl.io.X.trialMean_touch_responseR - nanmean(mdl.io.X.trialMean_touch_responseR)) ./ nanstd(mdl.io.X.trialMean_touch_responseR); %standardization
% mdlmeantouch = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt); 

DmatX = (mdl.io.X.trialSum_touch_responseR - nanmean(mdl.io.X.trialSum_touch_responseR)) ./ nanstd(mdl.io.X.trialSum_touch_responseR); %standardization
mdlmeantouch = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt); 

DmatX = (mdl.io.X.trialMean_touch_responseR - nanmean(mdl.io.X.trialMean_touch_responseR)) ./ nanstd(mdl.io.X.trialMean_touch_responseR); %standardization
DmatX = DmatX(:,sum(~isnan(DmatX))>0);
shuffDmatX = reshape(DmatX(randperm(numel(DmatX))),size(DmatX)); 
mdlmeantouchShuff = multinomialModel(mdl,shuffDmatX,mdl.io.Y.normal,glmnetOpt); 

decoderPerformance(mdlmeantouch)

mdls = {mdlmeantouchShuff,mdlnull,mdlpertouch,mdlavail,mdlmeantouch};
for d = 1:length(mdls)
    meanError(d) = mean(abs(mdls{d}.io.trueXpredicted(:,1) - mdls{d}.io.trueXpredicted(:,2)));
    stdError(d) = std(abs(mdls{d}.io.trueXpredicted(:,1) - mdls{d}.io.trueXpredicted(:,2)));
    binAccuracy(d) = mean(mdls{d}.gof.modelAccuracy);
    errs(:,d) = abs(mdls{d}.io.trueXpredicted(:,1) - mdls{d}.io.trueXpredicted(:,2)); 
end
figure(100);clf
errorbar(1:length(mdls),meanError,stdError ./ sqrt(glmnetOpt.numIterations),'ko')
set(gca,'xlim',[0 length(mdls)+1],'xtick',1:length(mdls),'xticklabel',{'shuffled','out','perTouch','in','sumTrial'},'ytick',0:length(mdls),'ylim',[0 4])
ylabel('distance (mm) from correct location')

