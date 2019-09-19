function [tuneStruct] = object_location_quantification(uberarray,selectedCells,hilbert_feature,displayOpt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function is used to plot a heat map of location (e.g. angle at touch)
% tuning for selectedCells. Specifically, use touch cells.
%
% inputs:
% uberarray - packaged uber array with all recorded units
% selectedCells - indices of units in uberarray that are touch cells
%
% outputs:
% whisk location tuning as defined from
% tuneStruct = struct with calculations of tuning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 4), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));


%function parameters
viewWindow = [-25:50]; %viewing window around touch
numTouchesPerBin = 75; %number of touches to assign in each bin for quantification.
alpha_value = .05; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 5; %smoothing parameter for smooth f(x) in shadedErrorBar
min_bins = 5; %minimum number of angle bins to consider quantifying
gauss_filt = .5; %2D smoothing filter for heatmap

%dependent function to id all touches and pre/post decision touches
preDecisionTouches = preDecisionTouchMat(uberarray);


%populating struct for tuning quantification
tuneStruct = cell(1,length(uberarray));
for i = 1:length(uberarray)
    tuneStruct{i}.calculations = [];
    tuneStruct{i}.is_tuned = nan;
end

if willdisplay
    rc = numSubplots(numel(selectedCells));
    figure(22);clf
    figure(23);clf
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
    
    rw = find(viewWindow == array.meta.touchProperties.responseWindow(1)) : find(viewWindow == array.meta.touchProperties.responseWindow(2));
    response = mean(tVar.allTouches.R_ntk(:,rw),2) * 1000;
    numBins = round(numel(selected_feature)./numTouchesPerBin);
    
    %% Heatmap of responses around touch. Window around touch defined as view_window
    if strcmp(hilbert_feature,'phase')
        [sorted_heat, sortedBy_heat] = binslin(selected_feature,tVar.allTouches.R_ntk,'equalE',13,-pi,pi);
    else
        [sorted_heat, sortedBy_heat] = binslin(selected_feature,tVar.allTouches.R_ntk,'equalN',numBins);
    end
    
    if willdisplay
        
        if numel(sortedBy_heat)>min_bins
            heat_resp = cell2mat(cellfun(@(x) mean(x,1),sorted_heat,'uniformoutput',0));
            smoothed_heat_resp = imgaussfilt(heat_resp,gauss_filt,'padding','replicate');
            
            
            figure(22);subplot(rc(1),rc(2),rec)
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
    end
    %% Tuning in touch response window
%     if strcmp(hilbert_feature,'phase') %commenting out for modulation
%     index calculation. Need all variables to have same amount of bins
%         [sorted, sortedBy] = binslin(selected_feature,response,'equalE',13,-pi,pi);
%     else
        [sorted, sortedBy] = binslin(selected_feature,response,'equalN',numBins);
%     end
    
    quant_ol_p = anova1(cell2nanmat(sorted),[],'off');
    
    SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),sorted);
    tscore = cellfun(@(x) tinv(.95,numel(x)-1),sorted);
    CI = SEM.*tscore;
    
    % ONLY FOR MOD IDX CALCULATIONS
    smooth_response = smooth(cellfun(@mean,sorted),smoothing_param);
    [maxResponse,~] = max(smooth_response);
    minResponse = min(smooth_response);
    tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = (maxResponse - minResponse) ./ (maxResponse + minResponse);
    
    % making sure we've sampled enough bins before plotting.
    % min_bins defined as a global param above.
    if numel(sortedBy)>min_bins
        
        if willdisplay
            figure(23);subplot(rc(1),rc(2),rec)
            shadedErrorBar(cellfun(@median, sortedBy), smooth(cellfun(@mean,sorted),smoothing_param),smooth(CI,smoothing_param),'k')
        end
        
        if quant_ol_p < alpha_value
            %plot smoothed response first
            smooth_response = smooth(cellfun(@mean,sorted),smoothing_param);
            if strcmp(array.meta.touchProperties.responseType,'excited')
                [maxResponse,idx] = max(smooth_response);
                minResponse = min(smooth_response);
                
                %plot scatter of first sig diff from max and max value
                sd_p = nan(length(sorted),1);
                pThresh = .05;
                for g = 1:numel(sorted)
                    [~,sd_p(g)] = ttest2(sorted{idx},sorted{g});
                end
                all_idx = find(sd_p < pThresh); %all points sig diff from max
                groupIdx = sort([all_idx ;idx]);
                mid_idx = find(groupIdx==idx);
                [~,sd_idx_tmp] = min(abs(idx - all_idx));
                sd_idx = all_idx(sd_idx_tmp);
                %finding tuning width before peak
                if (mid_idx-1) ~= 0
                    left_idx = groupIdx(mid_idx-1);
                    tuneStruct{selectedCells(rec)}.calculations.tune_peak = median(sortedBy{idx}); %peak modulation defined as the median value of the max bin
                    tuneStruct{selectedCells(rec)}.calculations.tune_left_width = median(sortedBy{idx}) - median(sortedBy{left_idx}); %width defined as the first bin that's sig diff from peak response
                else
                    left_idx = [];
                    tuneStruct{selectedCells(rec)}.calculations.tune_left_width = nan;
                end
                
                % finding tuning width after peak
                if (mid_idx+1) <= length(groupIdx)
                    right_idx = groupIdx(mid_idx+1);
                    tuneStruct{selectedCells(rec)}.calculations.tune_peak = median(sortedBy{idx}); %peak modulation defined as the median value of the max bin
                    tuneStruct{selectedCells(rec)}.calculations.tune_right_width = median(sortedBy{right_idx}) - median(sortedBy{idx}); %width defined as the first bin that's sig diff from peak response
                else
                    right_idx = [];
                    tuneStruct{selectedCells(rec)}.calculations.tune_right_width = nan;
                end
                

                
                if ~isempty(sd_idx) && ~isempty(idx)
                    if willdisplay
                        hold on; scatter(median(sortedBy{idx}),maxResponse,'b','filled');
                        hold on; scatter(median(sortedBy{sd_idx}),smooth_response(sd_idx),'b','filled');
                        
                    end
                    %calculations for output of tuning. Built specifically for
                    %touch excited units. Touch inhibited may need a new
                    %calculation
                    tuneStruct{selectedCells(rec)}.is_tuned = 1;
                    tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = (maxResponse - minResponse) ./ (maxResponse + minResponse);
                    tuneStruct{selectedCells(rec)}.calculations.mod_idx_abs = (maxResponse - minResponse) ;
                    tuneStruct{selectedCells(rec)}.calculations.responses_at_peak = sorted{idx};
                end
                
                
                
            elseif strcmp(array.meta.touchProperties.responseType,'inhibited')
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
                
                if ~isempty(sd_idx) && ~isempty(idx) && willdisplay
                    hold on; scatter(median(sortedBy{idx}),minResponse,'r','filled');
                    hold on; scatter(median(sortedBy{sd_idx}),smooth_response(sd_idx),'r','filled');
                end
            end
            
        end
        
        if willdisplay
            if strcmp(hilbert_feature,'phase')
                set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
            elseif strcmp(hilbert_feature,'pole')
                set(gca,'xlim',[-1 1],'xtick',-1:1:1,'xdir','reverse')
            else
                set(gca,'xlim',[min(cellfun(@median, sortedBy)) max(cellfun(@median, sortedBy))])
            end
        end
        
        tuneStruct{selectedCells(rec)}.stim_response.varNames = {'median S_ctk','mean R_ntk','std R_ntk','95CI R_ntk'};
        tuneStruct{selectedCells(rec)}.stim_response.values = [cellfun(@nanmedian, sortedBy) smooth(cellfun(@nanmean,sorted),smoothing_param) smooth(cellfun(@nanstd,sorted),smoothing_param) CI];
        tuneStruct{selectedCells(rec)}.stim_response.raw_stim = sortedBy;
        tuneStruct{selectedCells(rec)}.stim_response.raw_response = sorted;
    else
        tuneStruct{selectedCells(rec)}.is_tuned  = .5;  %not enough samples
    end
    
    
end