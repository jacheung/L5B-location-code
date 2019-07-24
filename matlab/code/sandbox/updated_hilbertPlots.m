%% pearson correlation of predicted FR and true FR

rawFR = cellfun(@(x) mean(x.io.DmatY),glmModel) *1000; 
predictedFR = cellfun(@(x) mean(x.predicted.spikeProb(:)),glmModel) * 1000; 

figure(123);clf
scatter(predictedFR,rawFR,50,'.k')
xlabel('predicted firing rate (Hz)');
ylabel('raw firing rate (Hz)')
title(['pearson correlation = ' num2str(corr(predictedFR',rawFR'))])



%%
glmModel{1}.predicted.angles

%% fitting to OL tuning
% tunedCells = cellfun(@(x) x.meta,glmModel);
% object_location_quantification(U,tunedCells,'angle');

numTouchesPerBin = 150; %number of touches to assign in each bin for quantification.
alpha_value = .05; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 5; %smoothing parameter for smooth f(x) in shadedErrorBar
min_bins = 5; %minimum number of angle bins to consider quantifying
gauss_filt = .5;

figure(22);clf
figure(24); clf

for rec = 1:length(glmModel)
    array = U{glmModel{rec}.meta};
    buildIndices = glmModel{rec}.modelParams.buildIndices;
    touchResponse = array.meta.responseWindow;
    trColumns = find(buildIndices==touchResponse(1)) : find(buildIndices==touchResponse(2));
    
    predResponseRaw = glmModel{rec}.predicted.spikeProb;
    predResponse = mean(glmModel{rec}.predicted.spikeProb(:,trColumns),2) * 1000;
    
    stimulus = glmModel{rec}.predicted.angles;
    numBins = round(numel(stimulus)./numTouchesPerBin);
    
    %% Heatmap
    
    [sorted_heat, sortedBy_heat] = binslin(stimulus,predResponseRaw,'equalN',numBins);
    
    
    if numel(sortedBy_heat)>min_bins
        figure(22);subplot(4,8,rec)
        heat_resp = cell2mat(cellfun(@(x) mean(x,1),sorted_heat,'uniformoutput',0));
        smoothed_heat_resp = imgaussfilt(heat_resp,gauss_filt,'padding','replicate');
        imagesc(smoothed_heat_resp)
        hold on; plot([find(buildIndices==0) find(buildIndices==0)],[1 length(sortedBy_heat)],'-.w')
        caxis([0 prctile(smoothed_heat_resp(:),100)])
        
        set(gca,'ydir','normal','ytick',1:length(sortedBy_heat),'yticklabel',round(cellfun(@median,sortedBy_heat)),...
            'xtick',1:25:length(buildIndices),'xticklabel',0:25:50)
        
    else
        disp(['Not plotting cell ' num2str(selectedCells(rec)) ' b/c less than 5 sampled bins'])
    end
    
    %% model OL tuning
    
    [sorted, sortedBy] = binslin(stimulus,predResponse,'equalN',numBins);
    quant_ol_p(rec) = anova1(cell2nanmat(sorted),[],'off');
    
    SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),sorted);
    tscore = cellfun(@(x) tinv(.95,numel(x)-1),sorted);
    CI = SEM.*tscore;
    
    if numel(sortedBy)>min_bins
        figure(24);subplot(3,8,rec)
        shadedErrorBar(cellfun(@median, sortedBy), smooth(cellfun(@mean,sorted),smoothing_param),smooth(CI,smoothing_param),'k')
        
        if quant_ol_p(rec) < alpha_value
            %plot smoothed response first
            smooth_response = smooth(cellfun(@mean,sorted),smoothing_param);
            if strcmp(array.meta.responseType,'excited')
                [maxResponse,idx] = max(smooth_response);
                
                %plot scatter of first sig diff from max
                sd_p = nan(length(sorted),1);
                pThresh = .05;
                for g = 1:numel(sorted)
                    [~,sd_p(g)] = ttest2(sorted{idx},sorted{g});
                end
                all_idx = find(sd_p < pThresh); %all points sig diff from max
                [~,sd_idx_tmp] = min(abs(idx - all_idx));
                sd_idx = all_idx(sd_idx_tmp);
                
                if ~isempty(sd_idx) && ~isempty(idx)
                    hold on; scatter(median(sortedBy{idx}),maxResponse,'b','filled');
                    hold on; scatter(median(sortedBy{sd_idx}),smooth_response(sd_idx),'b','filled');
                end
                
            elseif strcmp(array.meta.responseType,'inhibited')
                [minResponse,idx] = min(smooth_response);
                
                %plot scatter of first sig diff from min
                sd_p = nan(length(sorted),1);
                pThresh = .05;
                for g = 1:numel(sorted)
                    [~,sd_p(g)] = ttest2(sorted{idx},sorted{g});
                end
                all_idx = find(sd_p < pThresh); %all points sig diff from min
                [~,sd_idx_tmp] = min(abs(idx - all_idx));
                sd_idx = all_idx(sd_idx_tmp);
                
                if ~isempty(sd_idx) && ~isempty(idx)
                    hold on; scatter(median(sortedBy{idx}),minResponse,'r','filled');
                    hold on; scatter(median(sortedBy{sd_idx}),smooth_response(sd_idx),'r','filled');
                end
            end
            
            is_tuned(selectedCells(rec)) = 1;
        end
        
        
        set(gca,'xlim',[min(cellfun(@median, sortedBy)) max(cellfun(@median, sortedBy))])
        
    else
        is_tuned(selectedCells(rec)) = .5;  %not enough samples
    end
end

%% kernels
figure(128);clf

for i = 1:length(glmModel)
    subplot(3,8,i)
    BI = glmModel{i}.modelParams.buildIndices;
    coeffsToPlot = fields(glmModel{i}.coeffs);
    for u = 2:length(coeffsToPlot)
        coeffsToPlotName = coeffsToPlot{u};
        if strcmp(coeffsToPlotName,'touch') || strcmp(coeffsToPlotName,'touchDur')
            hold on; plot(BI,sum(glmModel{i}.coeffs.(coeffsToPlotName)'.*glmModel{i}.basisFunctions.touch,2))
        else
            hold on; plot(BI,sum(glmModel{i}.coeffs.(coeffsToPlotName)'.*glmModel{i}.basisFunctions.features,2))
        end
    end
    set(gca,'xtick',[-25:25:50])
end
%     ylabel('kernel weight')
%     xlabel('time from touch onset (ms)')
legend(coeffsToPlot(2:end))
suptitle('kernel weights in spike prediction')