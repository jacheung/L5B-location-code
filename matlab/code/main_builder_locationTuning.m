%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)~=0);
pole_tuned = object_location_quantification(U,selectedCells,'pole'); %for old see object_location_v1.0

%% population at touch pole decoding
% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 3;
glmnetOpt.numIterations = 20;
glmnetOpt.pctSamplingThreshold = .80; %what percent of total pole positions must be sampled before using
glmnetOpt.interpResolution = 40; %10mm / numBins (e.g. 40 = .25mm resolution)
glmnetOpt.samplingOption = 'poisson';
glmnetOpt.numberResamples = 50;

fileName = 'glm_location_decoder';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] =designMatrixBlocks_poleDecoder(glmModel,pole_tuned,glmnetOpt);
end

mdlResults = {};

DmatXraw = cell2mat(cellfun(@(x) x.io.DmatX,glmModel,'uniformoutput',0));
mdlResults.io.Xnorm = (DmatXraw - mean(DmatXraw)) ./ std(DmatXraw);
mdlResults.io.Y.normal = glmModel{1}.io.DmatY;

mdlResults = multinomialModel(mdlResults,mdlResults.io.Xnorm,mdlResults.io.Y.normal,glmnetOpt);
%% decoding resolution and probability of guesses
decoderPerformance(mdlResults);
usedUnits = cellfun(@(x) x.params.cellNum,glmModel);

suptitle([ 'Per touch location decoding using ' num2str(size(DmatXnorm,2)) ' tuned units'])

%% number of neurons vs the above
%mdl needs to have distance_from_true = cellfun(@(x,y) abs(x-y),mdl.io.trueY,mdl.io.predY,'uniformoutput',0);
%mdl needs to have confusion matrix  mdl.gof.confusionMatrix ./ sum(mdl.gof.confusionMatrix);
numNeurons = [1 5 10 20 30];
numIterations = 20;

mdl_mean = median(cell2mat(cellfun(@(x) x(:),mdlResults.fitCoeffs,'uniformoutput',0)),2);
reshaped_coeffs = reshape(mdl_mean,size(mdlResults.fitCoeffs{1}));
 
nframe = numNeurons;
v = VideoWriter('resolution_heatmap.avi');
open(v)
% mov(1:nframe)=struct('cdata',[],'colormap',[]);

resamp_mdl = [];
for g = 1:length(numNeurons)
    for d = 1:numIterations
        neuron_to_use = datasample(1:numel(usedUnits),numNeurons(g));
        
        coeffs_to_use = reshaped_coeffs([1 neuron_to_use+1],:);
        dmatX = [ones(length(mdlResults.io.Y.normal),1) mdlResults.io.Xnorm(:,neuron_to_use)];
        dmatY = mdlResults.io.Y.normal;
        
        predicts = dmatX * coeffs_to_use;
        probability =  1 ./ (1+exp(predicts*-1)); %convert to probability by using mean function (for binomial, it is the sigmoid f(x) 1/1+exp(-predicts))
        
        %highest probability = predict class
        [~,pred{d}] = max(probability,[],2);
        true{d} = dmatY;
    end
    
    resamp_mdl{g}.numNeurons = numNeurons(g);
    resamp_mdl{g}.io.trueY = true';
    resamp_mdl{g}.io.predY = pred';
    gofmetrics{g} = decoderPerformance(resamp_mdl{g});
    suptitle(['Number of neurons = ' num2str(numNeurons(g))])
    frame=getframe(gcf);
    writeVideo(v,frame);
end
close(v)
v.FrameRate = 1;



boneMap = flipud(jet(length(numNeurons)));
figure(80);clf
for d = 1:length(numNeurons)
    hold on;
    %     shadedErrorBar(gofmetrics{d}.resolution(:,1),gofmetrics{d}.resolution(:,2),gofmetrics{d}.resolution(:,3),'linecolor','b');
    plot(gofmetrics{d}.resolution(:,1),gofmetrics{d}.resolution(:,2),'color',boneMap(d,:));
end
set(gca,'ylim',[0 1],'xlim',[0 8],'xtick',0:2:10,'xticklabel',0:.5:5,'ytick',0:.25:1) %hard coded xticklabels for single touch prediction of pole position using population of OL tuned cells
xlabel('mm w/in prediction');ylabel('p (prediction)')
axis square
title('cdf of prediction varying # of neurons')
legend([num2str(numNeurons')])




