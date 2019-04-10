%% QUANTIFICATION OF MODEL
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\pop.mat')
OLcells = find(~(cellfun(@isempty,glmModel))); 
for i = OLcells(2:end)
    
    array = U{selectedCells(i)};
%PLOTTING BASIS FUNCTIONS
figure(480);clf
subplot(2,1,1)
plot(buildIndices,touchShiftIdx,'k')
subplot(2,1,2)
plot(buildIndices,touchShiftIdxRaw,'k')

%PLOTTING FOR time from touch onset
%     figure(800);subplot(5,6,i)
figure(484);clf
subplot(2,1,1)
%     touchSpks = spikes(repmat(touchIdx,1,length(buildIndices))+repmat(buildIndices,size(touchIdx,1),1));
%     hold on; plot(buildIndices,smooth(mean(touchSpks*1000,1)),'k')
hold on;plot(buildIndices,smooth(mean(glmModel{i}.predicted.spikeTestRaw *1000,1),length(basisFunction)),'k','linewidth',2)
hold on;plot(buildIndices,mean(glmModel{i}.predicted.spikeProb *1000,1),'b','linewidth',2)
ylabel('firing rate (hz)')
title(['mean psth of real vs model. Rsqd = ' num2str(corr(smooth(mean(glmModel{i}.predicted.spikeTestRaw *1000,1)'),mean(glmModel{i}.predicted.spikeProb *1000,1)').^2)]);
set(gca,'xtick',[])

figure(484);subplot(2,1,2);
hold on;
plot(buildIndices,smooth(std(glmModel{i}.predicted.spikeTestRaw )*1000,length(basisFunction)),'k')
plot(buildIndices,std(glmModel{i}.predicted.spikeProb )*1000,'b')
title('variance of real(black) vs model(blue)')
%hold on;bar(startIndices,mean(cell2mat(sHeat')*1000,1),'c')
%ylabel('firing rate (hz)')
%title(['mean real fr = ' num2str(mean(touchSpks(:))*1000) 'Hz'])
xlabel('time from touch onset')

thetabounds = [-100:2:100];

%plotting of original heatmap
%toss touches out that are >10 degree difference from nearest bin.


    gaussFilt = [2];
    
    toKeep = diff(pop{selectedCells(i)}.theta.range)<10;
    if ~isempty(find(toKeep==0))
        toKeep(find(toKeep==0):end) = 0;
    end
    
    figure(5480);clf
    subplot(2,3,1)
    imagesc(imgaussfilt(pop{selectedCells(i)}.theta.spikes(toKeep,:),gaussFilt,'padding','replicate'));
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',1:sum(toKeep),'yticklabel',[pop{selectedCells(i)}.theta.range(toKeep)],...
        'xtick',(0:25:length(buildIndices)),'xticklabel',[min(buildIndices):25:max(buildIndices)],'xlim',[0 length(buildIndices)]);
    for k=1:size(pop{selectedCells(i)}.theta.counts(toKeep),1)
        text(20,k,num2str(pop{selectedCells(i)}.theta.counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([find(buildIndices==0)-1 find(buildIndices==0)-1],[length(pop{selectedCells(i)}.theta.range(toKeep)) 0],'w:')
    ylabel('angle at touch')
    title('original tuning curve')
    
    selIdx = [find(buildIndices==array.meta.responseWindow(1)):find(buildIndices==array.meta.responseWindow(2))];
    hold on; plot([selIdx(1) selIdx(1)],[0 30],'-.w')
    hold on; plot([selIdx(end) selIdx(end)],[0 30],'-.w')

    %plotting of modeled heat map
    [sorted, sortedBy ,binBounds]=binslin(glmModel{i}.predicted.angles,glmModel{i}.predicted.spikeProb,'equalE',numel(thetabounds)+1,thetabounds(1),thetabounds(end));
    
    selBins = ~cellfun(@isempty ,sortedBy);
    spksBinned = cell2mat(cellfun(@(x) mean(x,1),sorted,'uniformoutput',0));
    
    modelBounds = bounds(selBins);
    heatPlot = spksBinned(selBins,:);
    ogBounds = pop{selectedCells(i)}.theta.range(toKeep);
    toPlotBounds = find(ismember(modelBounds,ogBounds));
    normHeatPlot =(heatPlot - min(heatPlot(:))) ./ (max(heatPlot(:))-min(heatPlot(:)));
    
    
    figure(5480);
    subplot(2,3,4)
    imagesc(imgaussfilt(normHeatPlot(toPlotBounds,:),gaussFilt,'padding','replicate'))
    set(gca,'ytick',1:length(toPlotBounds),'yticklabel', modelBounds(toPlotBounds),'ydir','normal',...
        'xtick',(0:25:length(buildIndices)),'xticklabel',[min(buildIndices):25:max(buildIndices)],'xlim',[0 length(buildIndices)])
    hold on;plot([find(buildIndices==0)-1 find(buildIndices==0)-1],[length(modelBounds) 0],'w:')
    caxis([0 prctile(normHeatPlot(:),99)]) %cap limit at 99th percentile in case of outliers in predictions
    xlabel('time from touch onset (ms)')
    ylabel('angle at touch')
    title('glmNet result')
    hold on; plot([selIdx(1) selIdx(1)],[0 30],'-.w')
    hold on; plot([selIdx(end) selIdx(end)],[0 30],'-.w')
    
    
    %PLOTTING LOCATION TUNING CURVE
    subplot(2,3,[2 3 5 6])
    colors = {'k','r'};
    
    spikesGroup = {glmModel{i}.raw.spikes , glmModel{i}.predicted.spikeProb};
    angleGroup = {glmModel{i}.raw.angle, glmModel{i}.predicted.angles};
    
    for d = 1:length(spikesGroup)
        allspks = spikesGroup{d};
        angleatT = angleGroup{d};
        touchResp = U{selectedCells(i)}.meta.responseWindow;
        
        
        bl = allspks(:,1:find(buildIndices==0));
        response = allspks(:,find(buildIndices==touchResp(1)):find(buildIndices==touchResp(2)));
        
        meanbase =  mean(mean(bl,2));
        stdbase = std(mean(bl,2));
        meanresp = mean(response,2);
        zscore = (meanresp - meanbase) ./ stdbase;
        
        [sorted,~,~] = binslin(angleatT,zscore,'equalE',numel(thetabounds),thetabounds(1),thetabounds(end));
        
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
        popzscore{i} = cellfun(@mean,sorted);
        otune(i)=p;
        
        x=cellfun(@str2num,stats.gnames);
        
        %CI BINS
        cibins = nan(size(x,1),1);
        
        SEMg = nanstd(zraw,[],2) ./ sqrt(sum(~isnan(zraw),2));
        for p = 1:length(x)
            rawx = zraw(x(p),:);
            SEM = SEMg(x(p));
            ts = tinv(.95,sum(~isnan(rawx),2)-1);      % T-Score
            cibins(p,:) = ts.*SEM;   %confidence intervals
        end
        
        y = stats.means;
        x=thetabounds(x)';
        
        xy = [x y'];
        
        hold on; filex_shadedErrorBar(xy(:,1),xy(:,2),cibins,['-' colors{d}]);
        set(gca,'ytick',round(min(xy(:,2))-1):round(max(xy(:,2))+1),'ylim',[round(min(xy(:,2))-1) round(max(xy(:,2))+1)],...
            'xlim',[min(xy(:,1)) max(xy(:,1))],'xtick',-20:20:60)
        
        zresp{d} = xy;
    end
    
    [~,interIdx ] = intersect(zresp{1}(:,1),zresp{2}(:,1));
    [~,interIdx2 ] = intersect(zresp{2}(:,1),zresp{1}(:,1));
    
    
    xlabel('angle at touch');ylabel('z-scored response')
    title(['raw(black) vs model(red) tuning. expVar=' num2str(corr(zresp{1}(interIdx,2),zresp{2}(interIdx2,2)).^2)])

    expVar(i) = corr(zresp{1}(interIdx,2),zresp{2}(interIdx2,2)).^2; 
end

%Plotting feature weights
% coeffsIter = cell2mat(sFitCoeffs);
% touchCoeffs = mean(coeffsIter(2:2+size(touchShiftIdx,2)-1,:),2);
% phaseCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
% ampCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
% midpointCoeffs = mean(coeffsIter(2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd ):2+size(touchShiftIdx,2)+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )+length(glmnetOpt.bf.indicesToAdd )-1,:),2);
% figure(9);clf
% hold on; plot(startIndices,sum(touchShiftIdx.*touchCoeffs',2))
% hold on; plot(startIndices,sum(touchShiftIdxRaw.*phaseCoeffs',2))
% hold on; plot(startIndices,sum(touchShiftIdxRaw.*ampCoeffs',2))
% hold on; plot(startIndices,sum(touchShiftIdxRaw.*midpointCoeffs',2))
% set(gca,'xlim',[-25 50],'xtick',-25:25:50)
% ylabel('kernel weight')
% xlabel('time from touch onset (ms)')
% legend('touch','phase','amp','midpoint')
%

figure(9);clf
hold on; plot(-25:50,sum(glmModel{i}.coeffs.touchCoeffs'.*glmModel{i}.basisFunctions.touch,2))
hold on; plot(-25:50,sum(glmModel{i}.coeffs.phaseCoeffs'.*glmModel{i}.basisFunctions.features,2))
hold on; plot(-25:50,sum(glmModel{i}.coeffs.ampCoeffs'.*glmModel{i}.basisFunctions.features,2))
hold on; plot(-25:50,sum(glmModel{i}.coeffs.midpointCoeffs'.*glmModel{i}.basisFunctions.features,2))
% hold on; plot(-25:50,sum(glmModel{i}.coeffs.curvatureCoeffs'.*glmModel{i}.basisFunctions.features,2))
set(gca,'xlim',[-25 50],'xtick',-25:25:50)
ylabel('kernel weight')
xlabel('time from touch onset (ms)')
legend('touch','phase','amp','midpoint','curvature')