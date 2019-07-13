function object_location_quantification(uberarray,selectedCells)
%this function is used to plot a heat map of location (angle at touch)
%tuning for selectedCells. Specifically, use touch cells. 
%
%inputs: 
%uberarray - packaged uber array with all recorded units
%selectedCells - indices of units in uberarray that are touch cells
%
%outputs:
%heatmap for object location tuning across time
%object location tuning in touch response window as defined from
%defTouchResponse.m function

%function parameters
viewWindow = [-25:50]; %viewing window around touch 
numTouchesPerBin = 75; %number of touches to assign in each bin for quantification. 
alpha_value = .05; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 5; %smoothing parameter for smooth f(x) in shadedErrorBar
min_bins = 5; %minimum number of angle bins to consider quantifying 

%dependent function for to id all touches and pre/post decision touches
preDecisionTouches = preDecisionTouchMat(uberarray);

quant_ol_p = nan(length(selectedCells),1);
figure(22);clf
figure(23);clf

for rec = 1:length(selectedCells)
    %stimulus and response variables definitions
    array = uberarray{selectedCells(rec)};
    [tVar] = atTouch_sorter(array,viewWindow,preDecisionTouches{selectedCells(rec)});
    angles = tVar.allTouches.S_ctk(:,1);
    rw = find(viewWindow == array.meta.responseWindow(1)) : find(viewWindow == array.meta.responseWindow(2));
    response = mean(tVar.allTouches.R_ntk(:,rw),2) * 1000;
    numBins = round(numel(angles)./numTouchesPerBin);
    
    %heatmap  
    [sorted_heat, sortedBy_heat] = binslin(angles,tVar.allTouches.R_ntk,'equalN',numBins);
    if numel(sortedBy_heat)>min_bins
        figure(22);subplot(4,8,rec)
        imagesc(cell2mat(cellfun(@(x) mean(x),sorted_heat,'uniformoutput',0)))
        hold on; plot([find(viewWindow==0) find(viewWindow==0)],[1 length(sortedBy_heat)],'-.w')
        set(gca,'ydir','normal','ytick',1:length(sortedBy_heat),'yticklabel',round(cellfun(@median,sortedBy_heat)),...
            'xtick',1:25:length(viewWindow),'xticklabel',-25:25:50)
    end
     
    %OL tuning in touch response window 
    [sorted, sortedBy] = binslin(angles,response,'equalN',numBins);
    quant_ol_p(rec) = anova1(cell2nanmat(sorted),[],'off');
    
    SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),sorted);
    tscore = cellfun(@(x) tinv(.95,numel(x)-1),sorted);
    CI = SEM.*tscore;
    
    if numel(sortedBy)>min_bins
        figure(23);subplot(4,8,rec)
        shadedErrorBar(cellfun(@median, sortedBy), smooth(cellfun(@mean,sorted),smoothing_param),smooth(CI,smoothing_param),'k')

        if quant_ol_p(rec) < alpha_value
            %plot scatter of max response
            smooth_response = smooth(cellfun(@mean,sorted),smoothing_param); 
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
                hold on; scatter(median(sortedBy{idx}),maxResponse,'g','filled');
                hold on; scatter(median(sortedBy{sd_idx}),smooth_response(sd_idx),'g','filled');
            end
            
        end
        
        set(gca,'xlim',[min(cellfun(@median, sortedBy)) max(cellfun(@median, sortedBy))])
        
    end
    
    
end