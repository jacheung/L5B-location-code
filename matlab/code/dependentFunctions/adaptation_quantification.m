function [adaptation] = adaptation_quantification(U,selectedCells,displayOpt)

if (nargin < 3), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

if willdisplay
    figure(80);clf
end

rc = numSubplots(numel(selectedCells));
%%
adaptation = cell(1,numel(U));
g_vy = cell(1,numel(selectedCells));
for rec = 1:length(selectedCells)
    array=U{selectedCells(rec)};
    array.onsetLatency=array.meta.touchProperties.responseWindow(1);
    trange=1:diff(array.meta.touchProperties.responseWindow);
    % Adaptation across touches
    % Adaptation by ITI
    
    spikes = squeeze(array.R_ntk);
    allTouchIdx = find(nansum([array.S_ctk(9,:,:);array.S_ctk(12,:,:)]));
    firstTouchIdx = find(array.S_ctk(9,:,:)==1);
    lateTouchIdx = find(array.S_ctk(12,:,:)==1);
    allTouchITI = diff([0; allTouchIdx]);
    
    tspikesIdx = repmat(allTouchIdx,1,length(trange))+repmat(trange+array.onsetLatency-1,length(allTouchIdx),1); %index for touch + trange
    
    if isnan(tspikesIdx)
        sc = 0;
        lsc = 0;
    else
        sc = sum(spikes(tspikesIdx),2);
        lspikesIdx = repmat(lateTouchIdx,1,length(trange))+repmat(trange+array.onsetLatency-1,length(lateTouchIdx),1);
        lsc = sum(spikes(lspikesIdx),2);
    end
    
    %% Touch Adaptation by touch order
    touchNumber = [];
    for i = 1:length(firstTouchIdx)-1
        touchNumber = cat(2,touchNumber,1:sum(allTouchIdx >= firstTouchIdx(i) & allTouchIdx < firstTouchIdx(i+1)));
    end
    touchNumber = cat(2,touchNumber,1:sum(allTouchIdx >= firstTouchIdx(end)));
    lateTouchNumber = touchNumber(touchNumber>1);
    adaptation{selectedCells(rec)}.lh = [];
    adaptation{selectedCells(rec)}.lci = [];
    
    [lh, lci] = poissfit(sc(touchNumber == 1));
    adaptation{selectedCells(rec)}.lh(1) = lh;
    adaptation{selectedCells(rec)}.lci(1,:) = lci;
    [adaptation{selectedCells(rec)}.sorted, adaptation{selectedCells(rec)}.sortedBy, adaptation{selectedCells(rec)}.binBounds]=binslin(lateTouchNumber, lsc, 'equalN',9);
    
    for num = 1:length(adaptation{selectedCells(rec)}.sorted)
        [lh, lci] = poissfit(adaptation{selectedCells(rec)}.sorted{num});
        adaptation{selectedCells(rec)}.lh(num+1) = lh;
        adaptation{selectedCells(rec)}.lci(num+1,:) = lci;
    end
    
    [adaptation{selectedCells(rec)}.calculations.p,~,~]= anova1(cell2nanmat(adaptation{selectedCells(rec)}.sorted),[],'off');
    
    %multcompare(stats,[],'off')
    if willdisplay
        if adaptation{selectedCells(rec)}.calculations.p<.05
            
            subplot(rc(1),rc(2),selectedCells(rec))
            shadedErrorBar(1:size(adaptation{selectedCells(rec)}.lh,2),adaptation{selectedCells(rec)}.lh,adaptation{selectedCells(rec)}.lci(:,2) - adaptation{selectedCells(rec)}.lh','b')
        else
            subplot(rc(1),rc(2),selectedCells(rec))
            shadedErrorBar(1:size(adaptation{selectedCells(rec)}.lh,2),adaptation{selectedCells(rec)}.lh,adaptation{selectedCells(rec)}.lci(:,2) - adaptation{selectedCells(rec)}.lh','k')
        end
    end
    
    maxResponse = max(adaptation{selectedCells(rec)}.lh);
    minResponse = min(adaptation{selectedCells(rec)}.lh);
    adaptation{selectedCells(rec)}.calculations.mod_idx_relative = (maxResponse - minResponse) ./ (maxResponse + minResponse);
    
    %% Touch Adaptation by ITI
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REMOVE THIS BEORE RUNNING. BUILT FOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING
    %     rec=datasample(1:length(selectedCells),1);
%     rec = 7
   
    array = U{(selectedCells(rec))};
    smooth_param = 5;
    alpha_value = 0.05;
    num_touches_per_bin = 75;
    
    spks = squeeze(array.R_ntk(:,:,:));
    response_window = array.meta.touchProperties.responseWindow(1):array.meta.touchProperties.responseWindow(2);
    touch_times = [find(array.S_ctk(9,:,:)==1) ; find(array.S_ctk(12,:,:)==1)];
    vis_touch_response = spks(touch_times + (-25:50));
    calc_touch_response = spks(touch_times + response_window);
    
    numBins = round(numel(touch_times)./num_touches_per_bin);
    
    if numBins>=2
        ITI = abs(touch_times  - [nan(1,1) ; touch_times(1:end-1)]);
        ITI(1) = max(ITI);
        ITI(ITI>4000) = 4000;
        
        [sorted_vis,sortedBy_vis] = binslin(ITI,vis_touch_response*1000,'equalN',numBins);
        [sorted_calc,~] = binslin(ITI,calc_touch_response*1000,'equalN',numBins);
        
        minx = cellfun(@min,sortedBy_vis);
        maxx = cellfun(@max,sortedBy_vis);
        medianx = cellfun(@median,sortedBy_vis);
        y = cellfun(@nanmean,sorted_vis,'uniformoutput',0) ;
        err = cellfun(@(x) nanstd(x) ./ sqrt(size(sorted_vis,1)),sorted_vis,'uniformoutput',0) ;
        
        mean_responses = cellfun(@(x) mean(sum(x,2)./numel(response_window)),sorted_calc);
        sem_responses = cellfun(@(x) std(sum(x,2)./numel(response_window)) ./ sqrt(size(x,1)),sorted_calc);
        
        
        [g_x,unique_idx] = unique(medianx);
        g_y = smooth(mean_responses(unique_idx),smooth_param);
        if numel(g_x>2)
            g_vy{rec} = interp1(g_x,g_y,50:25:4000);
        end
        
        
        if willdisplay
            figure(8);clf
            plot_elements = numBins-1;
            for g = 1:plot_elements
                subplot(3,plot_elements+1,g+1)
                shadedErrorBar(-25:50,smooth(y{g},smooth_param),smooth(err{g},smooth_param))
                title([num2str(minx(g)) '-' num2str(maxx(g))])
                set(gca,'ylim',[0 max(cell2mat(y),[],'all')],'xlim',[min(response_window) max(response_window)])
                
                if g == 1
                    subplot(3,plot_elements+1,1)
                    shadedErrorBar(-25:50,smooth(y{end},smooth_param),smooth(err{end},smooth_param))
                    title('"first touch"')
                    set(gca,'ylim',[0 max(cell2mat(y),[],'all')],'xlim',[min(response_window) max(response_window)])
                end
                
            end
            
            subplot(3,plot_elements+1,[(plot_elements+2) : (2*(plot_elements+1))])
            errorbar(medianx,smooth(mean_responses,smooth_param),smooth(sem_responses,smooth_param),'vertical','ko-')
            set(gca,'xlim',[0 max(medianx)])
            xlabel('ITI')
            ylabel('response window fr')
            suptitle(['cell num = ' num2str(selectedCells(rec))])
            
            subplot(3,plot_elements+1,[(plot_elements+plot_elements+3) : (3*(plot_elements+1))])
            f = fit(ITI,sum(calc_touch_response,2),'smoothingspline','SmoothingParam',.000001);
            hold on; plot(f,ITI,sum(calc_touch_response,2))
            set(gca,'xlim',[0 max(medianx)])
        end

    end
    
    
    
    %     firing_rate = sum(spks(touch_times + response_window),2) ./ numel(response_window) .* 1000;
    %
    %     if willdisplay
    %         figure(30);
    %         subplot(rc(1),rc(2),rec);
    %         scatter(ITI,firing_rate,'k.')
    %         set(gca,'xlim',[0 prctile(ITI,99)])
    %     end
    %
    
    
    
    %% ADAPTATION plotting raster by touch order
    %     [sort_adapt,idx_adapt] = sortrows ([allTouchIdx, touchNumber'],2);
    %     adapt_raster = zeros(length(idx_adapt),151);
    %     adapt_raster = spikes(repmat(allTouchIdx(idx_adapt),1,151)+repmat([-50:100],length(idx_adapt),1));
    %
    %
    
    %         figure(52);subplot(3,4,[1 3])
    %     plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    %     hold on
    %     plot(mod(find(adapt_raster'),151), ceil(find(adapt_raster')/151) ,'k.','markersize',4)
    %     axis([1 151 1 size(adapt_raster,1)]);
    %     colormap([1 1 1;0 0 0])
    %     box off
    %     set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[-25 0 50],'ytick',[])
    %     ax1 = gca;
    %     ax1;hold on
    %     ax1_pos = get(ax1,'Position');
    %     axes('Position',ax1_pos,...
    %         'XAxisLocation','top',...
    %         'yaxislocation','right','Color','none');hold on
    %     plot(mat2gray(sort_adapt(:,2)),1:length(sort_adapt),'r')
    %     axis tight
    %     norm_adapt_range = mat2gray([sort_adapt(1,2) sort_adapt(end,2)]);
    %     plot([norm_adapt_range(2)]*[1;1],[1 size(sort_adapt,1)],'--','color',[1 .5 .5])
    %     set(gca,'ytick',[],'xcolor','r','xtick',mat2gray([sort_adapt(1,2) sort_adapt(end,2)]),'xticklabel',round([sort_adapt(1,2) (sort_adapt(end,2))]))
    %     xlabel('Touch order in trial','color','r')
    %
    %     plot(adaptation.lh, 1:size(adaptation.lh,2),'k-')%(velocity.binBounds(1:end-1)+velocity.binBounds(2:end))/2, 'color','b')
    %     hold on
    %     box off
    %     set(gca,'ylim',[0.5 10.5],'ytick',[])
    %     plot(adaptation.lci,1:size(adaptation.lh,2), 'color',[.5 .5 .5])
    %     axes(ax1)
    %     set(gca,'color','none')
    
end


%POPULATION PLOTTING

if willdisplay
    population_heat_tuned = normalize_var(cell2mat(cellfun(@(x) x.lh,adaptation(tunedCells),'uniformoutput',0)')',0,1);
    population_heat_nontuned = normalize_var(cell2mat(cellfun(@(x) x.lh,adaptation(nontunedCells),'uniformoutput',0)')',0,1);
    [~,idx] = sort(population_heat_tuned(1,:));
    [~,idx_non] = sort(population_heat_nontuned(1,:));
    figure(81);clf
    subplot(1,2,1)
    imagesc(population_heat_tuned(:,fliplr(idx)))
    caxis([0 1])
    title('location tuned')
    subplot(1,2,2);
    imagesc(population_heat_nontuned(:,fliplr(idx_non)))
    caxis([0 1])
    title('non-location tuned')
    colormap gray
    xlabel('cell number')
    ylabel('touch order')
    colorbar
    set(gca,'ydir','reverse')
    
    
    ITI_mat = normalize_var(cell2mat(g_vy')',0,1);

    [~,maxidx] = max(ITI_mat,[],1);
    [~,idx] = sort(maxidx);
    
    figure(8324);clf
    imagesc(ITI_mat(:,idx)')
    set(gca,'xtick',1:4:41,'xticklabel',50:100:1000,'xlim',[1 40])
    xlabel('touch ITI')
    ylabel('neurons sorted by time to peak response') 
    caxis([0 1])
    colorbar
    
    
    figure(8880);clf
    hold on; plot(ITI_mat,'color',[.9 .9 .9])
    errorbar(1:size(ITI_mat,1),nanmean(ITI_mat,2),nanstd(ITI_mat,[],2) ./ sqrt(sum(~isnan(ITI_mat),2)),'k')
    hold on; errorbar(1:size(ITI_mat,1),nanmean(ITI_mat(:,idx(1:17)),2),nanstd(ITI_mat(:,idx(1:17)),[],2) ./ sqrt(sum(~isnan(ITI_mat(:,idx(1:17))),2)),'r')
        hold on; errorbar(1:size(ITI_mat,1),nanmean(ITI_mat(:,idx(18:end)),2),nanstd(ITI_mat(:,idx(18:end)),[],2) ./ sqrt(sum(~isnan(ITI_mat(:,idx(18:end))),2)),'b')
    set(gca,'xtick',1:4:size(ITI_mat,1),'xticklabel',50:100:1000,'xlim',[1 40])
    xlabel('touch ITI')
    ylabel('normalized touch response') 
    
end
