function sample_justification(uberarray,selectedCells,hilbert_feature)

whisk_num_bins = [200, 100, 50, 20, 10];
touches_per_bin = [10 25 50 75 100 125 150 175 200];

whisk_ol_p = cell(1,numel(selectedCells)); 
touch_ol_p = cell(1,numel(selectedCells)); 

for rec = 1:length(selectedCells)
    %% WHISKING 
    array = uberarray{selectedCells(rec)};
    spikes = squeeze(array.R_ntk);
    
    %whisking mask 
    timePostTouchToTrim = 30;
    touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
    touchOnIdx = touchOnIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchOffIdx = touchOffIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchEx_mask = ones(size(squeeze(array.S_ctk(1,:,:))));
    for i = 1:length(touchOnIdx)
        touchEx_mask(touchOnIdx(i):touchOffIdx(i)+timePostTouchToTrim) = NaN; %added 30 to add time from touch offset
    end
    touchEx_mask(1:100,:) = 1; %since bleedover from end of trials before, tihs ensure we keep end
    amplitude = squeeze(array.S_ctk(3,:,:));
    whisking = nan(size(squeeze(array.S_ctk(1,:,:))));
    whisking(amplitude>5)=1;
    whisking_mask = whisking .* touchEx_mask;

    if strcmp(hilbert_feature,'angle')
        conversion_feature = squeeze(array.S_ctk(1,:,:));
    elseif strcmp(hilbert_feature,'amplitude')
        conversion_feature = squeeze(array.S_ctk(3,:,:));
    elseif strcmp(hilbert_feature,'midpoint')
        conversion_feature = squeeze(array.S_ctk(4,:,:));
    elseif strcmp(hilbert_feature,'phase')
        conversion_feature = squeeze(array.S_ctk(5,:,:));
    end
    
    current_feature = conversion_feature(whisking_mask==1);
    filtered_spikes =spikes(whisking_mask==1);
    

    for b = 1:numel(whisk_num_bins)
        [sorted, ~] = binslin(current_feature,filtered_spikes*1000,'equalN',whisk_num_bins(b));
        nanmat = cell2nanmat(sorted);
        [whisk_ol_p{rec}(b),~,~] = anova1(nanmat,[],'off');
    end
    
    %% TOUCH
    %stimulus and response variables definitions
    touch_rw = array.meta.touchProperties.responseWindow;
    [tVar] = atTouch_sorter(array,touch_rw(1):touch_rw(2));
    
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
    
    response = mean(tVar.allTouches.R_ntk,2) * 1000;
    
    for b = 1:numel(touches_per_bin)
        numBins = round(numel(selected_feature)./touches_per_bin(b));
        if numBins >= 1
            [sorted_heat, ~] = binslin(selected_feature,response,'equalN',numBins);
            nanmat = cell2nanmat(sorted_heat);
            [touch_ol_p{rec}(b),~,~] = anova1(nanmat,[],'off');
        else
            touch_ol_p{rec}(b) = nan; 
        end
    end
    
end

%% plotting
p_mat = cell2mat(whisk_ol_p')';
figure(3840);clf

subplot(1,2,1);
plot(whisk_num_bins,p_mat,'color',[.8 .8 .8])
hold on; plot(whisk_num_bins,nanmean(p_mat,2),'r')
hold on; plot(whisk_num_bins,nanmedian(p_mat,2),'r--')
set(gca,'yscale','log','ylim',[0.0001 1],'xtick',fliplr(whisk_num_bins),'xticklabel',fliplr(1./whisk_num_bins)*100)
axis square
ylabel('ANOVA p-value');
xlabel('whisk data per bin (%)')

subplot(1,2,2);
t_mat = cell2mat(touch_ol_p')';
plot(touches_per_bin,t_mat,'color',[.8 .8 .8])
hold on; plot(touches_per_bin,nanmean(t_mat,2),'r')
hold on; plot(touches_per_bin,nanmedian(t_mat,2),'r--')
set(gca,'yscale','log','ylim',[0.0001 1],'xtick',touches_per_bin)
axis square
ylabel('ANOVA p-value');
xlabel('number of touches per bin')


