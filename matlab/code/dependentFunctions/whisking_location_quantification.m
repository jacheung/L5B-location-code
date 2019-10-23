function [tuneStruct] = whisking_location_quantification(uberarray,selectedCells,hilbert_feature,displayOpt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function is used to plot a heat map of whisking
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

rc = numSubplots(numel(selectedCells));

%function parameters
numWhiskSamplesPerBin = 5000; %number of whisks to assign in each bin for quantification.
alpha_value = .05; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 10; %smoothing parameter for smooth f(x) in shadedErrorBar
min_bins = 5; %minimum number of angle bins to consider quantifying

%dependent function to id all touches and pre/post decision touches
% preDecisionTouches = preDecisionTouchMat(uberarray(selectedCells));

%populating struct for tuning quantification
tuneStruct = cell(1,length(uberarray));
for i = 1:length(uberarray)
    tuneStruct{i}.is_tuned = nan;
    tuneStruct{i}.calculations = [];
end

all_masks = cellfun(@maskBuilder,uberarray(selectedCells)); %build masks
%%
for rec = 1:length(selectedCells)
    array = uberarray{selectedCells(rec)};
    curr_mask = all_masks(rec);
    whisking_mask = curr_mask.whisking .*curr_mask.touch;
    spikes = squeeze(array.R_ntk);
    
    if strcmp(hilbert_feature,'angle')
        conversion_feature = squeeze(array.S_ctk(1,:,:));
    elseif strcmp(hilbert_feature,'amplitude')
        conversion_feature = squeeze(array.S_ctk(3,:,:));
    elseif strcmp(hilbert_feature,'midpoint')
        conversion_feature = squeeze(array.S_ctk(4,:,:));
    elseif strcmp(hilbert_feature,'phase')
        conversion_feature = squeeze(array.S_ctk(5,:,:));
    elseif strcmp(hilbert_feature,'curvature')
        conversion_feature = squeeze(array.S_ctk(6,:,:));
    elseif strcmp(hilbert_feature,'pole') %conversion of angle to pole using touch positions
        viewWindow = -25:50;
        [tVar] = atTouch_sorter(array,viewWindow);
        pole_at_touch = normalize_var(tVar.allTouches.S_ctk(:,end),-1,1);
        angle_at_touch = tVar.allTouches.S_ctk(:,1);
        
        [s,sby] = binslin(pole_at_touch,angle_at_touch,'equalN',12); %bin responses based on pole positions
        cleaned = cell2mat(cellfun(@(x,y) rmoutliers([x y]),s,sby,'uniformoutput',0)); %remove outliers using median
        cleaned(sum(isnan(cleaned),2)>0,:) = [];
        clean_aat = cleaned(:,1);
        clean_pat = cleaned(:,2);
        p = polyfit(clean_aat,clean_pat,2); %use polyfit to convert angle to pole position
        
        angle = squeeze(array.S_ctk(1,:,:));
        conversion_feature = polyval(p,angle);
        
        % plotting conversion of whisker angle to pole
        %         figure(48);clf
        %         scatter(cleaned(:,2)*-1,cleaned(:,1),'.k')
        %         y = min(angle(:)):1:max(angle(:));
        %         hold on; plot(polyval(p,min(angle(:)):1:max(angle(:)))*-1,y,'r');
        %         set(gca,'xlim',[-2 4],'ytick',-20:20:80)
        %         axis square
        %         title(['cell num ' num2str(selectedCells(rec))])
        %         saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
        %         fn = ['angle2pole_' num2str(selectedCells(rec)) '.eps'];
        %         export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
        %         fix_eps_fonts([saveDir, fn])
        
        
        
    end
    
    current_feature = conversion_feature(whisking_mask==1);
    filtered_spikes =spikes(whisking_mask==1);
    
    filtered_spikes = circshift(filtered_spikes,4); %12ms lag
%     filtered_spikes =spikes(whisking_mask==1);
    
    
    numBins = round(sum(~isnan(current_feature(:)))./numWhiskSamplesPerBin);
    
    %% Tuning in touch response window
    %     if strcmp(hilbert_feature,'phase')
    %         [sorted, sortedBy] = binslin(current_feature,filtered_spikes*1000,'equalE',13,-pi,pi);
    %     else
    [sorted, sortedBy] = binslin(current_feature,filtered_spikes*1000,'equalN',numBins);
    %     end
    
    nonzero_bins = sum(cellfun(@mean, sorted)~=0);
    nanmat = cell2nanmat(sorted); 
    [quant_ol_p,~,stats] = anova1(nanmat,[],'off');
%     [shuff_ol_p,~,~] = anova1(reshape(randperm(numel(nanmat)),size(nanmat)),[],'off');
    
    SEM = cellfun(@(x) std(x) ./ sqrt(numel(x)),sorted);
    tscore = cellfun(@(x) tinv(.95,numel(x)-1),sorted);
    CI = SEM.*tscore;
    
    %BARSP FOR MOD IDX CALCULATIONS
    x = cellfun(@median,sortedBy);
    y = cellfun(@mean, sorted);
    
    if strcmp(hilbert_feature,'phase') || strcmp(hilbert_feature,'pole') 
        xq = min(cellfun(@median,sortedBy)):.1:max(cellfun(@median,sortedBy));
    else
        xq = min(cellfun(@median,sortedBy)):1:max(cellfun(@median,sortedBy));
    end
    yq = interp1(x,y,xq);
    
    barsFit = [];
    try
        barsFit = barsP(yq,[min(xq) max(xq)],numWhiskSamplesPerBin);
    catch
        disp(['skipping ' num2str(rec) ' due to ill fitting of bars'])
    end
    
    if ~isempty(barsFit) && nonzero_bins>10
%         figure(9);subplot(rc(1),rc(2),rec)
%         bar(x,y,'k')
%         if quant_ol_p<0.05
%             hold on; shadedErrorBar(xq(2:end-1),barsFit.mean(2:end-1),barsFit.confBands(2:end-1,2)-barsFit.mean(2:end-1),'g')
%         else
%             hold on; shadedErrorBar(xq(2:end-1),barsFit.mean(2:end-1),barsFit.confBands(2:end-1,2)-barsFit.mean(2:end-1),'k')
%         end
        %          hold on; plot(cellfun(@median,sortedBy),smooth(cellfun(@mean,sorted),smoothing_param),'r'); %smooth fitting
        
        
        smooth_response = barsFit.mean(2:end-1);
        smooth_stimulus = xq(2:end-1);

        [maxResponse,maxidx] = max(smooth_response);
        [minResponse,minidx] = min(smooth_response);
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = (maxResponse - minResponse) ./ (maxResponse + minResponse);
        tuneStruct{selectedCells(rec)}.calculations.tune_peak = smooth_stimulus(maxidx);
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_abs = (maxResponse - minResponse);
    else
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = 0;
        tuneStruct{selectedCells(rec)}.calculations.tune_peak = nan;
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_abs = 0;
    end
    
    
    smooth_response = smooth(cellfun(@mean,sorted),smoothing_param); %OLD
    smooth_stimulus = cellfun(@median,sortedBy);
    [maxResponse,maxidx] = max(smooth_response);
    [minResponse,minidx] = min(smooth_response);
    
    if numel(sortedBy)>min_bins
        
        if willdisplay
            figure(23);subplot(rc(1),rc(2),rec)
%             shadedErrorBar(smooth_stimulus,smooth_response,smooth(CI,smoothing_param),'k')
            shadedErrorBar(barsFit.mean(2:end-1), xq(2:end-1),barsFit.confBands(2:end-1,2)-barsFit.mean(2:end-1),'k')

        end
        
        if quant_ol_p < alpha_value
            %right and left tuning by idx
            compare_table = multcompare(stats,'Display','off');
            max_compares = compare_table(any(compare_table(:,[1 2]) == maxidx,2),:);
            sig_max_compares = max_compares(max_compares(:,end) < alpha_value,:);
            compare_idx = sig_max_compares(:,[1 2]);
            other_idx = compare_idx(compare_idx ~= maxidx);
            
            left_tuning_idx = other_idx(find(other_idx<maxidx,1,'last'));
            right_tuning_idx = other_idx(find(other_idx>maxidx,1,'first'));
            
            tuneStruct{selectedCells(rec)}.calculations.tune_peak = median(sortedBy{maxidx});
            
            %finding tuning width before peak
            if ~isempty(left_tuning_idx)
                tuneStruct{selectedCells(rec)}.calculations.tune_left_width = median(sortedBy{maxidx})-median(sortedBy{left_tuning_idx});
                if willdisplay
                    hold on; scatter(median(sortedBy{maxidx}),maxResponse,'g','filled');
                    hold on; scatter(median(sortedBy{left_tuning_idx}),smooth_response(left_tuning_idx),'r','filled');
                end
            else
                tuneStruct{selectedCells(rec)}.calculations.tune_left_width = nan;
            end
            
            %finding tuning width after peak
            if ~isempty(right_tuning_idx)
                tuneStruct{selectedCells(rec)}.calculations.tune_right_width = median(sortedBy{right_tuning_idx}) - median(sortedBy{maxidx});
                if willdisplay
                    hold on; scatter(median(sortedBy{maxidx}),maxResponse,'g','filled');
                    hold on; scatter(median(sortedBy{right_tuning_idx}),smooth_response(right_tuning_idx),'r','filled');
                end
            else
                tuneStruct{selectedCells(rec)}.calculations.tune_right_width = nan;
            end
            
            tuneStruct{selectedCells(rec)}.is_tuned = 1;
            tuneStruct{selectedCells(rec)}.calculations.responses_at_peak = sorted{maxidx};
            tuneStruct{selectedCells(rec)}.calculations.responses_at_trough = sorted{minidx};
            
        end
        
        if willdisplay
            if strcmp(hilbert_feature,'phase')
                set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'-\pi','0','\pi'})
            elseif strcmp(hilbert_feature,'pole')
                set(gca,'xtick',-5:1:5,'xdir','reverse')
            else
                set(gca,'xlim',[min(cellfun(@median, sortedBy)) max(cellfun(@median, sortedBy))])
            end
        end
        
        tuneStruct{selectedCells(rec)}.stim_response.varNames = {'median S_ctk','mean R_ntk','std R_ntk','95CI R_ntk'};
        tuneStruct{selectedCells(rec)}.stim_response.values = [cellfun(@nanmedian, sortedBy) smooth(cellfun(@nanmean,sorted),smoothing_param) smooth(cellfun(@nanstd,sorted),smoothing_param) smooth(CI,smoothing_param)];
        tuneStruct{selectedCells(rec)}.stim_response.raw_stim = sortedBy;
        tuneStruct{selectedCells(rec)}.stim_response.raw_response = sorted;
        if ~isempty(barsFit)
            tuneStruct{selectedCells(rec)}.stim_response.bars_fit= barsFit;
            tuneStruct{selectedCells(rec)}.stim_response.bars_stim = xq;
        end
    else
        tuneStruct{selectedCells(rec)}.is_tuned  = .5;  %not enough samples
    end
end
