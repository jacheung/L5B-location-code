load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\U.mat')

%% Top level parameters and definitions 
touchWindow = [-25:50]; %window for analyses around touch
numInterpPts = 24; %used for stretching or shrinking tuning curves to within the same bounds for decoding object location

touchCells = touchCell(U);
selectedCells = find(touchCells==1);

%% Structure for quantifying tuning and evaluating decoding 
popV = touchFeatureBinned(U,touchWindow);

%% Plotter for feature tuning around touch window
gaussFilt = 1; %smoothing function for tuning plots
fieldsList = fields(popV{1});
touchFeatureBinned_plotter(U,popV,selectedCells,fieldsList(1),touchWindow,gaussFilt)

%% Plotter for object location tuning
fieldsList = fields(popV{1});
tunedCellsIdx = tuningQuantification(U,popV,selectedCells,fieldsList(1),touchWindow);

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

%% Builder for identifying hilbert components that generate tuning
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

%% Builder for hilbert source 
selectedArray = U(selectedCells);

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 1;

%BUILD parameters
glmnetOpt.buildIndices = [-25:50]; %Indices around touch

%basis function and convolution using gaussian distribution
glmnetOpt.bf.bfwidth =7;
glmnetOpt.bf.bfstd = 5;
glmnetOpt.bf.bfspacing = 3;
basisFunction = normalize_var(normpdf(-1*glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth,0,glmnetOpt.bf.bfstd),0,1);
glmnetOpt.bf.indicesToAdd  = [-33:glmnetOpt.bf.bfspacing:20];

%GLM
fileName = 'glmModelTest';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = []; 
end

[glmModel] = designMatrixBuilder_hilbert(selectedArray,glmnetOpt,glmModel);

for i = 1:length(DmatX) 
glmModel{i} = binomialModel_hilbert(glmModel{i}.io.DmatXNormalized,glmModel{i}.io.DmatY,selectedArray{i},glmnetOpt,glmModel{i});

cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs')
save(fileName,'glmModel')
end
