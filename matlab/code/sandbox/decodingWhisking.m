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

%% Quantifying whisking tuning
whisking = whisking_general(U,'off');
hilbertWhisking = whisking_hilbert(U,popV,'off');

%% decoding whisker position

% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 10;

mdl = cell(length(U),1);
% figure;histogram(cell2mat(cellfun(@(x)
% sum(~isnan(x)),hilbertWhisking.R_ntk.phase,'uniformoutput',0)))
% %distribution of num samples in each bin
%using the median number of samples in each bin to dictate resampling size
resampNum = round(nanmedian(cell2mat(cellfun(@(x) sum(~isnan(x)),hilbertWhisking.R_ntk.phase,'uniformoutput',0))));

for i = 1:length(U)
    cellfeature = hilbertWhisking.R_ntk.phase{i};
    
    resampedSpikes = nan(resampNum,size(cellfeature,2));
    for b = 1:size(cellfeature,2)
        [~,bootsam] = bootstrp(20,@mean,cellfeature(:,b));
        binresampedSpikes = cellfeature(bootsam,b);
        binresampedSpikes = binresampedSpikes(~isnan(binresampedSpikes));
        
        resampedSpikes(:,b) = binresampedSpikes(1:resampNum);
    end
    
    boostedDmatX{i} = resampedSpikes;
    boostedDmatY{i} = repmat(1:size(cellfeature,2),resampNum,1);
end
    

numCellsToSample = 2;
selCells = datasample(1:length(U),numCellsToSample);
selCellsDmatX = boostedDmatX(selCells);
selCellsDmatY = boostedDmatY{selCells(1)};
DmatX = cell2mat(cellfun(@(x) x(:),selCellsDmatX,'uniformoutput',0));
DmatY = selCellsDmatY(:);







    





    
    
    







[glmModel] = designMatrixBuilder_whisking(mdl,glmnetOpt,hilbertWhisking);