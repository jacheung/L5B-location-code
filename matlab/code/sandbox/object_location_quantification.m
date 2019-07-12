preDecisionTouches = preDecisionTouchMat(U);

U = defTouchResponse(U,.95,'on');
selectedCells = find(cellfun(@(x) isfield(x.meta,'responseWindow'),U)==1);
%%
viewWindow = [-25:50];
numTouchesPerBin = 50;
alpha_value = .01;
quant_ol_p = nan(length(selectedCells),1);
smoothing_param = 5; 
close all;

for rec = 1:length(selectedCells)
    array = U{selectedCells(rec)};
    [tVar] = atTouch_sorter(array,viewWindow,preDecisionTouches{selectedCells(rec)});
    angles = tVar.allTouches.S_ctk(:,1);
    
    rw = find(viewWindow == array.meta.responseWindow(1)) : find(viewWindow == array.meta.responseWindow(2));
    response = mean(tVar.allTouches.R_ntk(:,rw),2) * 1000;
    
    numBins = round(numel(angles)./numTouchesPerBin);
    
    [sorted_heat, sortedBy_heat] = binslin(angles,tVar.allTouches.R_ntk,'equalN',numBins);
    figure(22);subplot(4,8,rec)
    imagesc(cell2mat(cellfun(@(x) mean(x),sorted_heat,'uniformoutput',0)))
    set(gca,'ydir','normal','ytick',1:length(sortedBy_heat),'yticklabel',round(cellfun(@median,sortedBy_heat)),...
        'xtick',1:25:length(viewWindow),'xticklabel',-25:25:50)
    
    
    [sorted, sortedBy] = binslin(angles,response,'equalN',numBins);
    quant_ol_p(rec) = anova1(cell2nanmat(sorted),[],'off');
    
    SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),sorted);
    tscore = cellfun(@(x) tinv(.95,numel(x)-1),sorted);
    CI = SEM.*tscore;
    
    if numel(sortedBy)>4
        figure(23);subplot(4,8,rec)
        shadedErrorBar(cellfun(@median, sortedBy), smooth(cellfun(@mean,sorted),smoothing_param),smooth(CI,smoothing_param),'k')

        if quant_ol_p(rec) < alpha_value
            %plot scatter of max response
            smooth_response = smooth(cellfun(@mean,sorted),smoothing_param); 
            [maxResponse,idx] = max(smooth_response);
            hold on; scatter(median(sortedBy{idx}),maxResponse,'g','filled');
            
            %plot scatter of first sig diff from max
            sd_p = nan(length(sorted),1);
            pThresh = .05; 
            for g = 1:numel(sorted)
            [~,sd_p(g)] = ttest2(sorted{idx},sorted{g});
            end
            all_idx = find(sd_p < pThresh); %all points sig diff from max
            [~,sd_idx_tmp] = min(abs(idx - all_idx));
            sd_idx = all_idx(sd_idx_tmp);
            
            hold on; scatter(median(sortedBy{sd_idx}),smooth_response(sd_idx),'g','filled');
           
        end
    end
    
    
    
end