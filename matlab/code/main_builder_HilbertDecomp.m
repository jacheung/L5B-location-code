%Load whisking and neural time series struct 
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\U.mat') %L5b excitatory cells
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
touchWindow = [-25:50]; %window for analyses around touch
numInterpPts = 24; %used for stretching or shrinking tuning curves to within the same bounds for decoding object location

touchCells = touchCell(U,'off');
selectedCells = find(touchCells==1);

% Structure for quantifying tuning and evaluating decoding 
popV = touchFeatureBinned(U,touchWindow);

% Defining touch response
U = defTouchResponse(U,.99,'off');
%% Plotter for object location tuning
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
tunedCells = tuningQuantification(U,popV,selectedCells,fieldsList(1),whichTouches,touchWindow,'off');

%% Builder for identifying hilbert components that generate tuning
selectedArray = U(find(tunedCells.mat(1,:)==1));

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
fileName = 'glmModelFullInterp';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] = designMatrixBlocks(selectedArray,glmnetOpt,glmModel);
end

%GLMdesign Matrix Build
selectedFeatures = [1 2 4 5 6 7 9]; 
interpOption = 'on'; %linear interpolation of missing values;
selectedFeaturesOptions = fields(glmModel{1}.io.components);
selectedFeaturesTitles = selectedFeaturesOptions(selectedFeatures)
[glmModel] = designMatrixBuilder_hilbert(glmModel,glmnetOpt,selectedFeatures,interpOption);

%Plot correlation matrix between features of design matrix. Look for
%obvious colinearity between features
% figure(62);clf
% imagesc(corr(glmModel{datasample(1:length(glmModel),1)}.io.DmatXNormalized))
% caxis([0 .7]) ;colorbar
% axis square; set(gca,'xtick',[],'ytick',[])

for i =13
 i
 glmModel{i} = binomialModel_hilbert(glmModel{i}.io.DmatXNormalized,glmModel{i}.io.DmatY,selectedArray{i},glmnetOpt,glmModel{i});
 glmModel{i}.meta = tunedCellsIdx{1}(i);
 glmModel{i}.name = fileName;
end

% cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs')
% save(fileName,'glmModel','-v7.3')
