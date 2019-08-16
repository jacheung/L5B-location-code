%Load whisking and neural time series struct 
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions 
% U = defTouchResponse(U,.95,'off');

%% Quantification of all units and their response properties
masks = cellfun(@(x) maskBuilder(x),U,'uniformoutput',0);
spks_in_touch = cellfun(@(x,y) nansum(nansum(squeeze(x.R_ntk).*(isnan(y.touch)))),U,masks);
spks_in_whisking = cellfun(@(x,y) nansum(nansum(squeeze(x.R_ntk).*y.whisking .*y.touch)),U,masks);
spks_in_all = cellfun(@(x) nansum(squeeze(x.R_ntk(:))),U);

mean_fr = cellfun(@(x) mean(squeeze(x.R_ntk(:))),U)*1000; %mean firing rate of all units 
prop_touch = spks_in_touch./spks_in_all; %proportion of spikes attributed to touch 
prop_whisking_touch = (spks_in_touch+spks_in_whisking)./spks_in_all;%proportion of spikes attributed to whisking + touch 
%% Quantification of touch excited units and their response properties 
touchUnits = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U)~=0);
%BUILD parameters for GLM to quantify touch units 
glmnetOpt.buildIndices = [-25:50]; %Indices around touch
glmnetOpt.touchDirection = 'protraction';
glmModel = [];
[glmModel] = designMatrixBlocks_simplified(U(touchUnits),glmnetOpt,glmModel);
excitThresh =  num2cell(cellfun(@(x) x.meta.touchProperties.baseline(1),U(touchUnits)));

%touch properties and definitions 
resp_onset = cellfun(@(x) x.meta.touchProperties.responseWindow(1),U(touchUnits)); % time of touch onset 
resp_window_length = cellfun(@(x) range(x.meta.touchProperties.responseWindow),U(touchUnits)); % length of touch response
touch_response_spks = cellfun(@(x,y) mean(x.io.DmatY),glmModel); %average number of spikes in touch response window 
pResponse_touch = cellfun(@(x,y) mean((x.io.DmatY*1000)>y),glmModel,excitThresh); %probability of generating spiking response > baseline + 95%CI

