load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\U.mat')

% Top level parameters and definitions 
touchWindow = [-25:50]; %window for analyses around touch
numInterpPts = 24; %used for stretching or shrinking tuning curves to within the same bounds for decoding object location

touchCells = touchCell(U);
selectedCells = find(touchCells==1);

% Structure for quantifying tuning and evaluating decoding 
popV = touchFeatureBinned(U,touchWindow);

%% Plotter for object location tuning
fieldsList = fields(popV{1});
tunedCellsIdx = tuningQuantification(U,popV,selectedCells,fieldsList(1),touchWindow);

%% Builder for identifying hilbert components that generate tuning
selectedArray = U(tunedCellsIdx);

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 5;

%BUILD parameters
glmnetOpt.buildIndices = [-25:50]; %Indices around touch

%basis function and convolution using gaussian distribution
glmnetOpt.bf.bfwidth =7;
glmnetOpt.bf.bfstd = 5;
glmnetOpt.bf.bfspacing = 3;
basisFunction = normalize_var(normpdf(-1*glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth,0,glmnetOpt.bf.bfstd),0,1);
glmnetOpt.bf.indicesToAdd  = [-33:glmnetOpt.bf.bfspacing:20];

%GLMdesign Matrix Set-up
fileName = 'glmModelFull';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = []; 
    [glmModel] = designMatrixBlocks(selectedArray,glmnetOpt,glmModel);
end

%GLMdesign Matrix Build
selectedFeatures = [1 3 4 5 ]; 

selectedFeaturesOptions = fields(glmModel{1}.io.components);
selectedFeaturesTitles = selectedFeaturesOptions([selectedFeatures])
[glmModel] = designMatrixBuilder_hilbert(glmModel,glmnetOpt,selectedFeatures);


for i = 8
    i 
glmModel{i} = binomialModel_hilbert(glmModel{i}.io.DmatXNormalized,glmModel{i}.io.DmatY,selectedArray{i},glmnetOpt,glmModel{i});

glmModel{i}.meta = tunedCellsIdx(i);
end

% cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs')
% save(fileName,'glmModel','-v7.3')
