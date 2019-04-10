function [DmatX,DmatYnorm,DmatYshuff] = designMatrixBuilder_touchFeature(U,popV,selectedCells,fields,touchWindow,numInterpPts)

%% LOCATION TUNING MODELING

[rc] = numSubplots(length(selectedCells));
plotrow = rc(1);
plotcolumn = rc(2);

bounds = popV{1}.(fields).bounds;


for d = 1:length(selectedCells)
    
    
    currCell = selectedCells(d);
    
    counts = popV{currCell}.(fields).counts;
    ranges = popV{currCell}.theta.range ;
    thetavect = [];
    for k = 1:length(ranges)
        thetavect = [thetavect ;ones(counts(k),1).*ranges(k)];
    end
    tvectNorm = normalize_var(thetavect,-1,1);
    
    
    %ZSCORE calculation
    allspikes = cell2mat(popV{currCell}.theta.raw);
    
    tOnset = find(touchWindow==0);
    responseWindow = U{currCell}.meta.responseWindow;
    meanbase = mean(mean(allspikes(:,1:tOnset),2));
    stdbase = std(mean(allspikes(:,1:tOnset),2));
    postresponses = allspikes(:,tOnset+responseWindow(1):tOnset+responseWindow(2));
    
    
    xresp = mean(postresponses,2);
    zscore = (xresp - meanbase) ./ stdbase;
    
    [sorted,~,~] = binslin(thetavect,zscore,'equalE',numel(bounds),bounds(1),bounds(end));
    
    samps = cellfun(@numel,sorted);
    selBins = find(samps>sum(samps)./100);
    
    
    zraw = nan(length(sorted),2000);
    for b = 1:length(sorted)
        if sum(b == selBins)>0
            currvals = sorted{b};
            if ~isempty(sorted{b})
                zraw(b,1:length(currvals)) = currvals';
            end
        end
    end
    
    
    [p,~,stats]=anova1(zraw',[],'off');
    vals =multcompare(stats,[],'off');
    
    
    x = str2num(cell2mat(stats.gnames));
    x=bounds(x)';
    y = stats.means;
    
    zstd = nanstd(zraw,[],2);
    ystd = zstd(~isnan(zstd))';
    
    xy = [x y'];
    
    %     kxy{d} = xy;
    %     figure(390);subplot(plotrow,plotcolumn,d);
    %     scatter(xy(:,1),xy(:,2),[],[.7 .7 .7],'filled')
    %     set(gca,'ytick',round(min(xy(:,2))-1):round(max(xy(:,2))+1),'ylim',[round(min(xy(:,2))-1) round(max(xy(:,2))+1)],...
    %         'xlim',[min(xy(:,1)) max(xy(:,1))],'xtick',-20:20:60)
    
    %     nxy = [normalize_var(x,0,1),y'];
    %     iy(:,d) = interp1(nxy(:,1),nxy(:,2),linspace(0,1,numInterpPts));
    %
    modelY(:,d) = interp1(normalize_var(x,0,1),y,linspace(0,1,numInterpPts));
    modelYstd(:,d) = interp1(normalize_var(x,0,1),ystd,linspace(0,1,numInterpPts));
    
    
end


    resampNum = 500;
    resampX = nan(size(modelY,1)*resampNum,size(modelY,2));
    for i = 1:resampNum
        resampX(size(modelY,1)*(i-1)+1:size(modelY,1)*(i-1)+size(modelY,1),:) =  normrnd(modelY,modelYstd);
    end
    DmatX = resampX;
    
    DmatYnorm = repmat([1:size(modelY,1)]',resampNum,1);
    randshuff = randperm(length(DmatYnorm));
    DmatYshuff = DmatYnorm(randshuff);

    
%     %     for d = 1:length(DmatY)
%     for d = 1
%         rando = randperm(length(DmatX));
%         tmpDmatX=DmatX(rando,:);
%         tmpDmatY=DmatY{d}(rando,:);
%         txp = [];
%         
%         cvglm
%         
%         
%         [thetas,cost,~] = ML_oneVsAll(tmpDmatX(1:end*.7,:),tmpDmatY(1:end*.7,:),numel(unique(DmatY{d})),0);
%         [pred,opt_thresh,prob]=ML_predictOneVsAll(thetas,tmpDmatX(end*.7:end,:)...
%             ,tmpDmatY(end*.7:end,:),'Max');
%         txp = [txp ; tmpDmatY(end*.7:end) pred];
%         
%         cmat = confusionmat(txp(:,1),txp(:,2));
%         ncmat = cmat./sum(cmat);
%         figure(9+d);clf;
%         imagesc(ncmat);
%         caxis([0 .40])
%         axis square
%         colorbar
%         title(['chance = ' num2str(1/size(modely,1))])
%         set(gca,'xtick',[],'ytick',[])
%         xlabel('truth');ylabel('predicted')
%         
%         predDiff(:,d) = txp(:,2) - txp(:,1);
%     end
%     model{b}=histcounts(predDiff(:,1),-16.5:16.5,'normalization','probability');
%     stats(b,:) = [model{b}(17) std(predDiff)];
%     %     shuff{b}=histcounts(predDiff(:,2),-16.5:16.5,'normalization','probability');
%     figure(40);clf
%     plot(-16:16,model{b},'k')
%     %     hold on; plot(-16:16,shuff{b},'color',[.7 .7 .7])
%     
%     set(gca,'xtick',-16:8:16,'ytick',0:.1:.5,'ylim',[0 .35])
%     %     decodeResolution = std(predDiff)* (mean(rangeExplored)./size(modely,1))
