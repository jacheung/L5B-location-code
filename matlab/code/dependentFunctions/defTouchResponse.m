function [U] = defTouchResponse(U,confidenceThreshold,displayOpt)

%function that builds upon location code data struct. Adds fields of
%touch properties under meta. Response window is defined by a confidence
%threshold above baseline response (BL response is from -25:0ms pre touch)

if (nargin < 3), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

if willdisplay
    rc = numSubplots(length(U));
    figure(3000);clf
end

heatTouch = cell(1,length(U)); 

for rec=1:length(U)
    array = U{rec};
    window = [-50:50];
    touchIdx= [find(array.S_ctk(9,:,:)==1) ; find(array.S_ctk(12,:,:)==1)];
    spks = squeeze(array.R_ntk);
    
    within_range = ~ (logical(sum((touchIdx + window) > numel(spks),2)) | logical(sum((touchIdx + window) < 0, 2)));
    touchIdx = touchIdx(within_range); 
       
    blIdx = window(find(window==-50):find(window==0));
    
    touchSpks = spks(touchIdx+window);
    baseline_spikes = spks(1:round(mean(array.meta.poleOnset)*1000),:);
    touchSpksShuff = spks(touchIdx+blIdx);
    
    %calculating x% confidence interval
    SEM = nanstd(touchSpksShuff(:))./ sqrt(sum(~isnan(touchSpksShuff(:,1))));
    ts = tinv(confidenceThreshold,sum(~isnan(touchSpksShuff(:,1))));
    CI = SEM.*ts;
%     SEM = nanstd(baseline_spikes(:))./ sqrt(sum(~isnan(baseline_spikes(:,1))));
%     ts = tinv(confidenceThreshold,sum(~isnan(baseline_spikes(:,1))));
%     CI = SEM.*ts;
    
    touchResponse = smooth(nanmean(touchSpks),5);
    excitThreshold = mean(touchSpksShuff(:)) + CI;
    inhibThreshold = mean(touchSpksShuff(:)) - CI;
%     excitThreshold = mean(baseline_spikes(:)) + CI; %redefining baseline as pre-pole period instead of pre-touch 
%     inhibThreshold = mean(baseline_spikes(:)) - CI;
    
    
    %Excitatory threshold defined as the period > mean+x%CI
    excitthreshIdx = touchResponse'>excitThreshold;
    excitthreshIdx(1:find(window==5))=0;
    
    %Inhibitory threshold defined as the period > mean-x%CI
    inhibthreshIdx = touchResponse'<inhibThreshold;
    inhibthreshIdx(inhibthreshIdx==1)= -1;
    inhibthreshIdx(1:find(window==5))=0;
    
    U{rec}.meta.touchProperties.baseline_varNames = {'excit base+95CI','inhib base-95CI'};
    U{rec}.meta.touchProperties.baseline = [excitThreshold inhibThreshold]*1000;
    
    heatTouch{rec} = nanmean(touchSpks)'; 
    U{rec}.meta.touchProperties.responseType = 'untuned'; %deeming neuron as touch untuned less otherwise
    
    post_touch_window = find(window==0) : find(window==max(window));
    [~,response_max_index] = max(heatTouch{rec}(post_touch_window));
    U{rec}.meta.touchProperties.peak_index = response_max_index;
    
    
    %quantifying modulation index
    max_touch = max(touchResponse(find(window==0):end));
    U{rec}.meta.touchProperties.mod_idx_relative = (max_touch - mean(baseline_spikes(:))) ./ (max_touch + mean(baseline_spikes(:)));
    
    %Defining touch responses as a period between two points that are less
    %than 5ms apart.
    if sum(excitthreshIdx)>0 && sum(excitthreshIdx)>sum(inhibthreshIdx) 
%     if sum(excitthreshIdx)>0 && min([find(excitthreshIdx) 100])<min([find(inhibthreshIdx) 100])
        tps = window(excitthreshIdx);
        if ~isempty(tps)
            startPoint = tps(1);
            endPoint = tps(find(diff(diff(tps)<10)==-1,1,'first'));
            if isempty(endPoint)
                endPoint = tps(end);
            end
            
            %used below to eliminate touch responses that are sig but way
            %too small.
            if mean(touchResponse(startPoint+find(window==0):endPoint+find(window==0))*1000) > 2 && endPoint - startPoint >= 4
                U{rec}.meta.touchProperties.responseType = 'excited';
                U{rec}.meta.touchProperties.responseWindow=[startPoint endPoint];
                
                U{rec}.meta.touchProperties.SNR(rec) = log(mean(touchResponse(startPoint+find(window==0):endPoint+find(window==0))*1000) ./ (excitThreshold*1000));
                if willdisplay
                    figure(3000);subplot(rc(1),rc(2),rec)
%                     hold on; bar(window,touchResponse*1000,'b','facealpha',.2,'edgealpha',.2);
                    hold on; bar(window,touchResponse*1000,'b');
                    hold on; scatter(window(startPoint+find(window==0):endPoint+find(window==0)),touchResponse(startPoint+find(window==0):endPoint+find(window==0))*1000,'b','filled')
                    hold on; plot(window,ones(length(window),1).* excitThreshold .* 1000,'k-.')
                    set(gca,'xtick',-25:25:50,'xlim',[-25 50])
                    title(num2str(U{rec}.meta.touchProperties.mod_idx_relative));
                end
            
                
                
            else
                if willdisplay
                    figure(3000);subplot(rc(1),rc(2),rec)
                    hold on; bar(window,touchResponse*1000,'k','facealpha',.2,'edgealpha',.2);
                    hold on; plot(window,ones(length(window),1).* excitThreshold .* 1000,'k-.')
                    set(gca,'xtick',-25:25:50)
                end
            end
        end
        
    elseif sum(excitthreshIdx)<sum(inhibthreshIdx)
        tps = window(inhibthreshIdx);
        if ~isempty(tps)
            startPoint = tps(1);
            endPoint = tps(find(diff(diff(tps)<10)==-1,1,'first'));
            if isempty(endPoint)
                endPoint = tps(end);
            end
            
            %used below to eliminate touch responses that are sig but small
            %firing rates (<2Hz) or are short (responses < 4ms)
            if mean(touchResponse(startPoint+find(window==0):endPoint+find(window==0))*1000) > 2 && endPoint - startPoint >= 4
                U{rec}.meta.touchProperties.responseType = 'inhibited';
                U{rec}.meta.touchProperties.responseWindow=[startPoint endPoint];
                
                U{rec}.meta.touchProperties.SNR (rec) = log(mean(touchResponse(startPoint+find(window==0):endPoint+find(window==0))*1000) ./ (inhibThreshold*1000));
                
                if willdisplay
                    figure(3000);subplot(rc(1),rc(2),rec)
                    hold on; bar(window,touchResponse*1000,'r','facealpha',.2,'edgealpha',.2);
                    hold on; scatter(window(startPoint+find(window==0):endPoint+find(window==0)),touchResponse(startPoint+find(window==0):endPoint+find(window==0))*1000,'r','filled')
                    hold on; plot(window,ones(length(window),1).* inhibThreshold .* 1000,'k-.')
                    set(gca,'xtick',-25:25:50)
                    title(num2str(U{rec}.meta.touchProperties.mod_idx_relative));
                end
                
            end
        end
        
        
        
    else
        if willdisplay
            figure(3000);subplot(rc(1),rc(2),rec)
            hold on; bar(window,touchResponse*1000,'k','facealpha',.2,'edgealpha',.2);
            hold on; plot(window,ones(length(window),1).* excitThreshold .* 1000,'k-.')
             set(gca,'xtick',-25:25:50)

        end

    end

    
end

%% summary plot of touch responses

if willdisplay
    untuned = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'untuned'),U));
    excited = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
    inhibited = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'inhibited'),U));
    
    responseTypes = {excited,inhibited,untuned};
    
    figure(3001);clf
    
    for g = 1:length(responseTypes)
        current = heatTouch(responseTypes{g});
        post_touch_window = find(window==0) : find(window==max(window));
        response_max_index = cellfun(@(x) find(x(post_touch_window) == max(x(post_touch_window))==1,1,'first'),current);
        [~,sorted_idx] = sort(response_max_index);
        norm_mat = normalize_var(cell2mat(current(sorted_idx)),0,1);
        
        if g == 1
            subplot(5,1,[1:3]);
        elseif g == 2
            subplot(5,1,4);
        elseif g ==3
            subplot(5,1,5);   
        end
        imagesc(norm_mat')
        hold on; plot([find(window==0) find(window==0)],[0 numel(current)+1],'-.w')
%         set(gca,'xtick',1:25:length(window),'xticklabel',min(window):25:max(window),'xlim',[1 length(window)])
        set(gca,'xtick',26:25:length(window),'xticklabel',-25:25:max(window),'xlim',[26 length(window)]) %set xlim -25:50;
        colorbar
    end
    
    
%    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
%     fn = 'touch_response.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])
end

    


