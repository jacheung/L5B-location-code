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
glmnetOpt.numIterations = 3;

% figure;histogram(cell2mat(cellfun(@(x)
% sum(~isnan(x)),hilbertWhisking.R_ntk.phase,'uniformoutput',0)))
% %distribution of num samples in each bin
%using the median number of samples in each bin to dictate resampling size
% resampNum = round(nanmedian(cell2mat(cellfun(@(x) sum(~isnan(x)),hilbertWhisking.R_ntk.phase,'uniformoutput',0))));
resampNum = 4000; 

boostedDmatX = cell(1,length(U)); 
boostedDmatY = cell(1,length(U));

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
    
numCellsToSample = [100 150 175 250 500 750 1000];

selCells = datasample(1:length(U),max(numCellsToSample));
resampCells = boostedDmatX(selCells);
shiftTuningDmatX = cellfun(@(x) x(:,circshift(1:size(x,2),randperm(size(x,2)))),resampCells,'uniformoutput',0);
shuffledResponseDmatX = cellfun(@(x) x(randperm(size(x,1)),:),shiftTuningDmatX,'uniformoutput',0);

decodingResolutionMean = zeros(1,length(numCellsToSample)); 
decodingResolutionSEM = zeros(1,length(numCellsToSample)); 

for g = 1:length(numCellsToSample)
    selCells = datasample(1:length(shuffledResponseDmatX),numCellsToSample(g));
    selCellsDmatX = shuffledResponseDmatX(selCells);

    selCellsDmatY = boostedDmatY{1};
    
    DmatX = cell2mat(cellfun(@(x) x(:),selCellsDmatX,'uniformoutput',0));
    DmatX = (DmatX - mean(DmatX)) ./ std(DmatX); %standardization
    DmatY = selCellsDmatY(:);
    
    mdl.io.X = DmatX;
    mdl.io.Y.normal = DmatY;
    
    mdl = multinomialModel(mdl,DmatX,DmatY,glmnetOpt);
    decoderPerformance(mdl)
    
    binResolution = 360/numel(unique(DmatY));
    decodingResolutionMean(g) = mean(abs(mdl.io.trueXpredicted(:,1) - mdl.io.trueXpredicted(:,2))) * binResolution;
    decodingResolutionSEM(g) = (std(abs(mdl.io.trueXpredicted(:,1) - mdl.io.trueXpredicted(:,2))) ./ glmnetOpt.numIterations) * binResolution;
    
end

figure(23);clf
shadedErrorBar(numCellsToSample,decodingResolutionMean,decodingResolutionSEM,'k')
set(gca,'ylim',[0 160],'xlim',[0 1000])
xlabel('number of cells')
ylabel('decoding resolution (degrees of phase)')



   