
%
builtUnits = find(cellfun(@(x) isfield(x,'gof'),glmModel));

%%
for i = 1:length(builtUnits)
    trial_length = glmModel{builtUnits(i)}.modelParams.trial_length;

        true = cellfun(@(x) reshape(x,trial_length,numel(x)./trial_length), glmModel{builtUnits(i)}.predicted.spikeTestRaw,'uniformoutput',0);
        predicted = cellfun(@(x) reshape(x,trial_length,numel(x)./trial_length), glmModel{builtUnits(i)}.predicted.spikeProb,'uniformoutput',0);
        polePositions = cell2mat(glmModel{builtUnits(i)}.predicted.pole);
    
        true_psth = mean(cell2mat(true)',2);
        predicted_psth =  mean(cell2mat(predicted)',2);
    
        [sort_values,idx] = sort(polePositions);
        full_true = cell2mat(true)'; 
        full_predict = cell2mat(predicted)';
        sorted_raster_true = full_true(idx,:);
%         sorted_raster_predict = poissrnd(full_predict(idx,:));
         sorted_raster_predict = full_predict(idx,:);
        
        
        figure(40);clf
        subplot(2,1,1)
        imagesc(sorted_raster_true)
                set(gca,'ytick',[])
        ylabel('close - far')
        subplot(2,1,2);
        imagesc(sorted_raster_predict)
        
        set(gca,'ytick',[])
        ylabel('close - far')
        
    meanGOF(i) = mean(glmModel{builtUnits(i)}.gof.devExplained)
end

