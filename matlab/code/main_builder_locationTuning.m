load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\U.mat')

load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat')
%%
% Top level parameters and definitions 
touchWindow = [-25:50]; %window for analyses around touch
numInterpPts = 24; %used for stretching or shrinking tuning curves to within the same bounds for decoding object location

touchCells = touchCell(U);
selectedCells = find(touchCells==1);

% Structure for quantifying tuning and evaluating decoding 
popV = touchFeatureBinned(U,touchWindow);

% Defining touch response
U = defTouchResponse(U,.99);
%% Plotter for feature tuning around touch window
gaussFilt = 1; %smoothing function for tuning plots
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
touchFeatureBinned_plotter(U,popV,selectedCells,fieldsList(1),whichTouches,touchWindow,gaussFilt)

%% Plotter for object location tuning
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
tunedCellsIdx = tuningQuantification(U,popV,selectedCells,fieldsList(1),whichTouches,touchWindow);

%optional raster of FRs for tuned cells. 
for d = 1
    figure(8);clf
    allSpks = squeeze(U{tunedCellsIdx(d)}.R_ntk);
    [~,idx] = sort(U{tunedCellsIdx(d)}.meta.motorPosition);
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
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

%design matrix for feature in fields list
fieldsList = fields(popV{1});
[mdl.io.X, mdl.io.Y.normal, mdl.io.Y.shuffled] = designMatrixBuilder_touchFeature(U,popV,selectedCells,fieldsList{1},touchWindow,numInterpPts);

%multinomial model for decoding location 
mdl = multinomialModel(mdl,mdl.io.X,mdl.io.Y.normal,glmnetOpt); %normalizing all model tuning to within x num interpolation points
decoderPerformance(mdl)



