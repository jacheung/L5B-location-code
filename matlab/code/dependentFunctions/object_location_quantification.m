function [tuneStruct] = object_location_quantification(uberarray,selectedCells,hilbert_feature)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%is_tuned = vector showing whether neuron is tuned to hilbert_feature at
%touch
%tc_xy = tuning curves showing x(stim) and y(responses); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%function parameters
viewWindow = [-25:50]; %viewing window around touch 
numTouchesPerBin = 75; %number of touches to assign in each bin for quantification. 
alpha_value = .05; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 5; %smoothing parameter for smooth f(x) in shadedErrorBar
min_bins = 5; %minimum number of angle bins to consider quantifying 
gauss_filt = .5; 

%dependent function to id all touches and pre/post decision touches
preDecisionTouches = preDecisionTouchMat(uberarray);
rc = numSubplots(numel(selectedCells));

quant_ol_p = nan(length(selectedCells),1);
figure(22);clf
figure(23);clf

tuneStruct = cell(1,length(uberarray));

for i = 1:length(uberarray)
    tuneStruct{i}.mod_depth = nan; 
    tuneStruct{i}.is_tuned = nan;
end


for rec = 1:length(selectedCells)
    %stimulus and response variables definitions
    array = uberarray{selectedCells(rec)};
    [tVar] = atTouch_sorter(array,viewWindow,preDecisionTouches{selectedCells(rec)});
    
    if ~isempty(hilbert_feature)
        if strcmp(hilbert_feature,'angle')
            selected_feature = tVar.allTouches.S_ctk(:,1);
        elseif strcmp(hilbert_feature,'amplitude')
            selected_feature = tVar.allTouches.S_ctk(:,3);
        elseif strcmp(hilbert_feature,'midpoint')
            selected_feature = tVar.allTouches.S_ctk(:,4);
        elseif strcmp(hilbert_feature,'phase')
            selected_feature = tVar.allTouches.S_ctk(:,5);
        elseif strcmp(hilbert_feature,'curvature')
            selected_feature = tVar.allTouches.S_ctk(:,6);
        elseif strcmp(hilbert_feature,'pole')
            selected_feature = normalize_var(tVar.allTouches.S_ctk(:,7),-1,1);
        else
            error('select features of "angle", "amplitude", "midpoint", "phase", "curvature", or "pole"')
        end
    else
        error('select features of "angle", "amplitude", "midpoint", "phase", "curvature", or "pole"')
    end
    
    rw = find(viewWindow == array.meta.responseWindow(1)) : find(viewWindow == array.meta.responseWindow(2));
    response = mean(tVar.allTouches.R_ntk(:,rw),2) * 1000;
    numBins = round(numel(selected_feature)./numTouchesPerBin);
    
    %% Heatmap
    if strcmp(hilbert_feature,'phase')
        [sorted_heat, sortedBy_heat] = binslin(selected_feature,tVar.allTouches.R_ntk,'equalE',13,-pi,pi);
    else
        [sorted_heat, sortedBy_heat] = binslin(selected_feature,tVar.allTouches.R_ntk,'equalN',numBins);
    end
    
    
    if numel(sortedBy_heat)>min_bins
        figure(22);subplot(rc(1),rc(2),rec)
        heat_resp = cell2mat(cellfun(@(x) mean(x,1),sorted_heat,'uniformoutput',0)); 
        smoothed_heat_resp = imgaussfilt(heat_resp,gauss_filt,'padding','replicate');
        imagesc(smoothed_heat_resp)
        hold on; plot([find(viewWindow==0) find(viewWindow==0)],[1 length(sortedBy_heat)],'-.w')
        caxis([0 prctile(smoothed_heat_resp(:),99)])
        if strcmp(hilbert_feature,'phase')
            set(gca,'ydir','normal','ytick',[1 12],'yticklabel',{'-\pi','\pi'},...
                'xtick',1:25:length(viewWindow),'xticklabel',-25:25:50)
        elseif strcmp(hilbert_feature,'pole')
            [~,zero_idx] = min(abs(cellfun(@median,sortedBy_heat)));
            set(gca,'ytick',zero_idx,'yticklabel',0,...
                'xtick',1:25:length(viewWindow),'xticklabel',-25:25:50)
        else
            set(gca,'ydir','normal','ytick',1:length(sortedBy_heat),'yticklabel',round(cellfun(@median,sortedBy_heat)),...
                'xtick',1:25:length(viewWindow),'xticklabel',-25:25:50)
        end

        
    else
        disp(['Not plotting cell ' num2str(selectedCells(rec)) ' b/c less than 5 sampled bins'])
    end
     
    %% Tuning in touch response window 
    if strcmp(hilbert_feature,'phase')
        [sorted, sortedBy] = binslin(selected_feature,response,'equalE',13,-pi,pi);
    else
        [sorted, sortedBy] = binslin(selected_feature,response,'equalN',numBins);
    end
    
    quant_ol_p(rec) = anova1(cell2nanmat(sorted),[],'off');
    
    SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),sorted);
    tscore = cellfun(@(x) tinv(.95,numel(x)-1),sorted);
    CI = SEM.*tscore;
    
    if numel(sortedBy)>min_bins
        figure(23);subplot(rc(1),rc(2),rec)
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
            
            
            tuneStruct{selectedCells(rec)}.mod_depth = (max(cellfun(@mean,sorted)) - min(cellfun(@mean,sorted))) ./ mean(cellfun(@mean,sorted));
            tuneStruct{selectedCells(rec)}.is_tuned = 1; 
        end
        
        if strcmp(hilbert_feature,'phase')
            set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
        elseif strcmp(hilbert_feature,'pole')
             set(gca,'xlim',[-1 1],'xtick',-1:1:1,'xdir','reverse')
        else
            set(gca,'xlim',[min(cellfun(@median, sortedBy)) max(cellfun(@median, sortedBy))])
        end
        
         tuneStruct{selectedCells(rec)}.stim_response = [cellfun(@median, sortedBy) smooth(cellfun(@mean,sorted),smoothing_param)];
         
    else
        tuneStruct{selectedCells(rec)}.is_tuned  = .5;  %not enough samples 
    end
    
    
end