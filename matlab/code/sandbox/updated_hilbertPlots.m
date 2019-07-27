%% pearson correlation of predicted FR and true FR
gauss_std = [1 2 4 8 16 32];
pearson_corr = zeros(length(gauss_std),length(glmModel));
for i = 1:length(gauss_std)
    rawFR = cellfun(@(x) imgaussfilt(x.predicted.spikeTestRaw(:),gauss_std(i)),glmModel,'uniformoutput',0);
    predictedFR = cellfun(@(x) x.predicted.spikeProb(:),glmModel,'uniformoutput',0);
    
    pearson_corr(i,:) = cellfun(@(x,y) corr(x,y),rawFR,predictedFR);
end
%% plotting
figure(123);clf
subplot(2,1,1)
xs = repmat((1:length(gauss_std))',1,length(glmModel)); 
scatter(xs(:),pearson_corr(:),'kx')
hold on; boxplot(pearson_corr')
set(gca,'xtick',1:length(gauss_std),'xticklabel',gauss_std)
xlabel('Gaussian sigma');
ylabel('Pearson correlation')

subplot(2,1,2)
[~,idx] = sort(pearson_corr(end,:));
imagesc(pearson_corr(:,idx))
set(gca,'ytick',1:length(gauss_std),'yticklabel',gauss_std,'xdir','reverse','ydir','normal')
colorbar
xlabel('cell number');
ylabel('Gaussian sigma')

%% deviance explained
DE = cellfun(@(x) mean(x.gof.devExplained),glmModel);
figure(5221);clf
for i = 1:length(gauss_std)
hold on;scatter(DE,pearson_corr(i,:),'filled')
end


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