%Load whisking and neural time series struct 
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
touchWindow = [-25:50]; %window for analyses around touch
numInterpPts = 24; %used for stretching or shrinking tuning curves to within the same bounds for decoding object location

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

numInterpPts = 20;

%design matrix for feature in fields list
fieldsList = fields(popV{1}.allTouches);
whichTouches = fields(popV{1});
[mdl.io.X, mdl.io.Y.normal, mdl.io.Y.shuffled] = designMatrixBuilder_touchFeature(U,popV,selectedCells,fieldsList{1},whichTouches(1),touchWindow,numInterpPts);

%multinomial model for decoding location 

% mdl = multinomialModel(mdl,mdl.io.X.responseTouchZ,mdl.io.Y.normal,glmnetOpt); %normalizing all model tuning to within x num interpolation points
DmatX = (mdl.io.X.outside_pole - nanmean(mdl.io.X.outside_pole)) ./ nanstd(mdl.io.X.outside_pole); %standardization
mdlnull = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt); %normalizing all model tuning to within x num interpolation points

DmatX = (mdl.io.X.pole_avail - nanmean(mdl.io.X.pole_avail)) ./ nanstd(mdl.io.X.pole_avail); %standardization
mdlavail = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt); %normalizing all model tuning to within x num interpolation points

DmatX = (mdl.io.X.touch_responseR - nanmean(mdl.io.X.touch_responseR)) ./ nanstd(mdl.io.X.touch_responseR); %standardization
mdlpertouch = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt); %normalizing all model tuning to within x num interpolation points

DmatX = (mdl.io.X.trialMean_touch_responseR - nanmean(mdl.io.X.trialMean_touch_responseR)) ./ nanstd(mdl.io.X.trialMean_touch_responseR); %standardization
mdlmeantouch = multinomialModel(mdl,DmatX,mdl.io.Y.normal,glmnetOpt); %normalizing all model tuning to within x num interpolation points

DmatX = (mdl.io.X.trialMean_touch_responseR - nanmean(mdl.io.X.trialMean_touch_responseR)) ./ nanstd(mdl.io.X.trialMean_touch_responseR); %standardization
mdlmeantouchShuff = multinomialModel(mdl,DmatX,mdl.io.Y.shuffled,glmnetOpt); %normalizing all model tuning to within x num interpolation points

selectedArray = popV(selectedCells);
decoderPerformance(mdlmeantouchShuff,selectedArray)

mdls = {mdlnull,mdlavail,mdlpertouch,mdlmeantouch,mdlmeantouchShuff};
for d = 1:length(mdls)
    meanError(d) = mean(abs(mdls{d}.io.trueXpredicted(:,1) - mdls{d}.io.trueXpredicted(:,2)));
    stdError(d) = std(abs(mdls{d}.io.trueXpredicted(:,1) - mdls{d}.io.trueXpredicted(:,2)));
    
    errs(:,d) = abs(mdls{d}.io.trueXpredicted(:,1) - mdls{d}.io.trueXpredicted(:,2)); 
end

figure(100);clf
errorbar(1:length(mdls),meanError*.25,stdError*.25,'ko')
set(gca,'xlim',[0 length(mdls)+1],'xtick',1:length(mdls),'xticklabel',{'out','in','perTouch','trialTouches','shuffled'},'ytick',0:1:3)
ylabel('distance (mm) from correct location')


