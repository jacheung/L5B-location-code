%Load whisking and neural time series struct 
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
% U = defTouchResponse(U,.95,'off');
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)~=0);
tStruct = object_location_quantification(U,selectedCells,'pole');

%%
touchCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)); 
trained = find(cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U));
trained = 1:length(U);

selectedUnits = intersect(touchCells,trained);
selectedArray = U(selectedUnits); 

%GLMNET parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0;
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

%BUILD parameters
glmnetOpt.touchDirection = 'pro';
glmnetOpt.touchOrder = 'all';

fileName = 'glm_trial_average';
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    glmModel = [];
    [glmModel] = designMatrixBlocks_trial_average(selectedArray,glmnetOpt,glmModel);
end

selectedFeatures = 1:7;
interpOption = 'off'; %linear interpolation of missing values;
selectedFeaturesOptions = fields(glmModel{1}.io.components);
selectedFeaturesTitles = selectedFeaturesOptions(selectedFeatures);
[glmModel] = designMatrixBuilder_trial_average(glmModel,glmnetOpt,selectedFeatures);


for i = 1:length(glmModel)
    disp(['iterating for neuron ' num2str(i) '/' num2str(length(selectedArray))])
    
    if size(glmModel{i}.io.DmatXNormalized,1)<20
        disp(['skipping neuron ' num2str(i) 'b/c too few trials'])
    else
        glmModel{i} = poissonModel(glmModel{i}.io.DmatXNormalized,glmModel{i}.io.DmatY,selectedArray{i},glmnetOpt,glmModel{i});
        glmModel{i}.meta = selectedUnits(i);
        glmModel{i}.name = fileName;
    end
end

%% heatmap of input features
for cellNum = 30

    [s_pole,idx] = sort(normalize_var(glmModel{cellNum}.raw.trimmedPole',1,-1));
    features = normalize_var(glmModel{cellNum}.io.DmatXNormalized(idx,:),-1,1);
    
    figure(658);clf
    imagesc([s_pole features]')
    
    stretch_resolution = 100 ;
    stretch_map = redbluecmap;
    new_map = nan(stretch_resolution,3);
    for j = 1:3
        new_map(:,j) = interp1(linspace(1,stretch_resolution,length(stretch_map)),stretch_map(:,j),1:stretch_resolution);
    end
    
    set(gca,'yticklabel',['pole' ; glmModel{1}.io.selectedFeatures.name])
    colorbar
    colormap(new_map)
    xlabel('trials sorted by pole position')
    
end

%% goodness of fit metrics 

builtUnits = find(cellfun(@(x) isfield(x,'gof'),glmModel));
meangof = cellfun(@(x) median(x.gof.devExplained),glmModel(builtUnits));
figure(230);subplot(2,1,1)
histogram(meangof,0:.1:1)
xlabel('gof: deviance explained');
ylabel('number of units'); 


wellfitUnits = intersect(builtUnits,find(meangof>.1)); 
bestFeatures = cell2mat(cellfun(@(x) mean(x.coeffs.raw(2:end,:),2),glmModel(wellfitUnits),'uniformoutput',0));
bfNorm = normalize_var(abs(bestFeatures),0,1);

figure(230);subplot(2,1,2)
shadedErrorBar(1:size(bestFeatures,1),nanmean(bfNorm,2),nanstd(bfNorm,[],2)./sum(~isnan(sum(bfNorm))))
set(gca,'xtick',1:length(glmModel{1}.io.selectedFeatures.name),'xticklabel',glmModel{1}.io.selectedFeatures.name); 
ylabel('normalized feature weight')
title(['n = ' num2str(sum(~isnan(sum(bfNorm)))) ' well fit units'])


figure(250);clf
[s_gof,idx] = sort(meangof);
plotUnits = builtUnits(fliplr(idx));
gof = fliplr(s_gof);

rc=numSubplots(length(plotUnits));
for g = 1:length(plotUnits)
    rec = plotUnits(g);

    subplot(rc(1),rc(2),g)

    if gof(g) > .1
        scatter(glmModel{rec}.predicted.spikeTestRaw,glmModel{rec}.predicted.spikeProb,'r.')
    else
        scatter(glmModel{rec}.predicted.spikeTestRaw,glmModel{rec}.predicted.spikeProb,'k.')
    end
    title(num2str(gof(g))); 
end
suplabel('raw fr (Hz)')
suplabel('predicted fr (Hz)','y')


