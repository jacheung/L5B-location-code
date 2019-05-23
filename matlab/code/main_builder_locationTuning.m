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
%% Plotter for feature tuning around touch window
gaussFilt = 1; %smoothing function for tuning plots
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
touchFeatureBinned_plotter(U,popV,selectedCells,fieldsList(1),whichTouches,touchWindow,gaussFilt)

%% Quantifying object location tuning
fieldsList = fields(popV{1}.allTouches);
tunedCellsIdx = tuningQuantification(U,popV,selectedCells,fieldsList(1),whichTouches(1),touchWindow,'on');

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
glmnetOpt.xfoldCV = 4;
glmnetOpt.numIterations = 5;

numInterpPts = 10;

%design matrix for feature in fields list
fieldsList = fields(popV{1}.allTouches);
whichTouches = fields(popV{1});

selectedCells = find(touchCells==1);
expCells = selectedCells(selectedCells>=29)
naiveCells = selectedCells(selectedCells<29)
[mdl.io.X, mdl.io.Y.normal, mdl.io.Y.shuffled] = designMatrixBuilder_touchFeature(U,selectedCells,whichTouches(1),touchWindow,numInterpPts);

%multinomial model for decoding location 

% mdl = multinomialModel(mdl,mdl.io.X.responseTouchZ,mdl.io.Y.normal,glmnetOpt); %normalizing all model tuning to within x num interpolation points
clear mdlnull
DmatX = mdl.io.X.outside_pole;
DmatX = (mdl.io.X.outside_pole - nanmean(mdl.io.X.outside_pole)) ./ nanstd(mdl.io.X.outside_pole); %standardization
mdlnull = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt); 

DmatX = (mdl.io.X.pole_avail - nanmean(mdl.io.X.pole_avail)) ./ nanstd(mdl.io.X.pole_avail); %standardization
mdlavail = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt);

DmatX = (mdl.io.X.touch_responseR - nanmean(mdl.io.X.touch_responseR)) ./ nanstd(mdl.io.X.touch_responseR); %standardization
mdlpertouch = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt);

DmatX = (mdl.io.X.trialMean_touch_responseR - nanmean(mdl.io.X.trialMean_touch_responseR)) ./ nanstd(mdl.io.X.trialMean_touch_responseR); %standardization
mdlmeantouch = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt); 

DmatX = (mdl.io.X.trialMean_touch_responseR - nanmean(mdl.io.X.trialMean_touch_responseR)) ./ nanstd(mdl.io.X.trialMean_touch_responseR); %standardization
mdlmeantouchShuff = multinomialModel(mdl,DmatX,mdl.io.Y.shuffled,glmnetOpt); 

selectedArray = popV(selectedCells);
decoderPerformance(mdlnull,selectedArray)

mdls = {mdlnull,mdlavail,mdlpertouch,mdlmeantouch,mdlmeantouchShuff};
for d = 1:length(mdls)
    meanError(d) = mean(abs(mdls{d}.io.trueXpredicted(:,1) - mdls{d}.io.trueXpredicted(:,2)));
    stdError(d) = std(abs(mdls{d}.io.trueXpredicted(:,1) - mdls{d}.io.trueXpredicted(:,2)));
    
    errs(:,d) = abs(mdls{d}.io.trueXpredicted(:,1) - mdls{d}.io.trueXpredicted(:,2)); 
end

figure(100);clf
errorbar(1:length(mdls),meanError,stdError,'ko')
set(gca,'xlim',[0 length(mdls)+1],'xtick',1:length(mdls),'xticklabel',{'out','in','perTouch','trialTouches','shuffled'},'ytick',0:1:5)
ylabel('distance (mm) from correct location')




clear mdlnull
DmatX = sum(mdl.io.X.outside_pole,2);
DmatY = mdl.io.Y.normal;
% DmatX = (mdl.io.X.outside_pole - nanmean(mdl.io.X.outside_pole)) ./ nanstd(mdl.io.X.outside_pole); %standardization
mdlnull = multinomialModel(mdl,DmatX(:),DmatY,glmnetOpt); 
decoderPerformance(mdlnull,selectedArray)

clear mdlnull
DmatX = sum(mdl.io.X.pole_avail,2);
DmatY = mdl.io.Y.normal;
mdlavail = multinomialModel(mdl,DmatX(:),DmatY,glmnetOpt); 
decoderPerformance(mdlavail,selectedArray)

clear mdlmeantouchShuff
DmatX = sum(mdl.io.X.trialMean_touch_responseR,2);
DmatY = mdl.io.Y.normal ;
mdlmeantouchShuff = multinomialModel(mdl,DmatX(:),DmatY,glmnetOpt); 
selectedArray = popV(selectedCells);
decoderPerformance(mdlmeantouchShuff,selectedArray)

