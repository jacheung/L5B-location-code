function [adaptation] = adaptation_quantification(U,selectedCells,displayOpt)

% selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));


if (nargin < 3), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

if willdisplay
    figure(80);clf
end

rc = numSubplots(numel(selectedCells));

adaptation = cell(1,numel(U));
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
    adaptation{selectedCells(rec)}.calculations.mod_idx_relative = (maxResponse - minResponse) ./ mean(adaptation{selectedCells(rec)}.lh);
    
    %% Touch Adaptation by ITI
    array = U{(selectedCells(rec))};
    
    spks = squeeze(array.R_ntk(:,:,:));
    response_window = array.meta.touchProperties.responseWindow(1):array.meta.touchProperties.responseWindow(2);
    touch_times = [find(array.S_ctk(9,:,:)==1) ; find(array.S_ctk(12,:,:)==1)];
    
    ITI = touch_times  - [nan(1,1) ; touch_times(1:end-1)];
    ITI(1) = max(ITI);
    ITI(ITI>1000) = 1000;
    
    firing_rate = sum(spks(touch_times + response_window),2) ./ numel(response_window) .* 1000;
    
    if willdisplay
        figure(30);
        subplot(rc(1),rc(2),g);
        scatter(ITI,firing_rate,'k.')
        set(gca,'xlim',[0 prctile(ITI,99)])
    end
    
    
    
    
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

population_heat = normalize_var(cell2mat(cellfun(@(x) x.lh,adaptation,'uniformoutput',0)')',0,1);
[~,idx] = sort(population_heat(1,:));
if willdisplay
    figure(81);clf
%     imagesc(population_heat)
    imagesc(population_heat(:,fliplr(idx)))
    colormap bone
    xlabel('cell number')
    ylabel('touch order')
    colorbar
end
