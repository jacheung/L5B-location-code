%% raster sorted by go and nogo
for i = 28
smoothTime = 100;
    
spikes = squeeze(U{i}.R_ntk);
go = find(U{i}.meta.trialType==1);
nogo = find(U{i}.meta.trialType==0);

figure(80);clf
shadedErrorBar(1:U{i}.t,smooth(nanmean(spikes(:,go) ,2)*1000,smoothTime),smooth( (nanstd(spikes(:,go),[] ,2)*1000) ./ numel(go) ,smoothTime),'b')
hold on;shadedErrorBar(1:U{i}.t,smooth(nanmean(spikes(:,nogo) ,2)*1000,smoothTime),smooth( (nanstd(spikes(:,nogo),[] ,2)*1000) ./ numel(nogo),smoothTime),'r')

set(gca,'xtick',0:1000:4000,'ytick',0:10:40)
end

%% features x pole position 
masks.preDecisionTouches = preDecisionTouchMat(U);
masks.allTouches = cellfun(@(x) ones(size(x)),preDecisionTouchesMask,'uniformoutput',false);

touchDirection = 'all';
touchOrder = 'all';
featureNumber = 1;

for i = 1:length(U)

    polePositions = normalize_var(U{i}.meta.motorPosition,1,-1);
    
    [numTouches, featureMatrix] = preDecisionTouchFeatures(U{i},masks.allTouches{i},featureNumber,touchDirection,touchOrder);

    [~,cols] = find(~isnan(featureMatrix));
    allPP = polePositions(cols);
    allFeat = featureMatrix(~isnan(featureMatrix));
    [~,idx] = sort(polePositions) ;
    
    
    closePoles = find(polePositions<.1 & polePositions>-.1);
    meanFeats = nanmean(featureMatrix); 
    dbtheta(i) = nanmean(meanFeats(closePoles));
    
%     figure(i);clf
%     errorbar(polePositions,nanmean(featureMatrix)-dbtheta,nanstd(featureMatrix),'ko')
%     hold on; scatter(allPP,allFeat-dbtheta,'k.');
%     hold on; scatter(polePositions,nanmean(featureMatrix)); 
    
    pt{i} = binslin(polePositions, nanmean(featureMatrix)-dbtheta(i),'equalE',11,-1,1);
    pnumT{i} = binslin(polePositions, numTouches,'equalE',11,-1,1);
end

groupThetas = cell2mat(cellfun(@(x) cellfun(@nanmean ,x),pt,'uniformoutput',false));
groupnumT = cell2mat(cellfun(@(x) cellfun(@nanmean ,x),pnumT,'uniformoutput',false));

figure(10);clf
subplot(1,2,1)
shadedErrorBar(linspace(-.9,.9,10),nanmean(groupThetas,2),nanstd(groupThetas,[],2))
set(gca,'xtick',[-1 0 1],'ytick',[-25:25:25],'ylim',[-25 25])
title(['mean DB ' num2str(mean(dbtheta)) '+/- ' num2str(std(dbtheta))])
axis square
subplot(1,2,2)
shadedErrorBar(linspace(-.9,.9,10),nanmean(groupnumT,2),nanstd(groupnumT,[],2))
set(gca,'xtick',[-1 0 1])
axis square
%% firing rate X depth of recording
figure(480);clf

jc_silent_cell = [766 819 895 631 776 815 910 871 844];

scatter( cellfun(@(y) y.meta.depth,U),cellfun(@(y) mean(y.R_ntk(:))*1000, U),'k','filled')
hold on; scatter(jc_silent_cell,ones(length(jc_silent_cell),1)*-1,[],'c','filled');
hold on; plot([700 700],[0 30],'-.k')
hold on; plot([900 900],[0 30],'-.k')
title(['n = response (' num2str(numel(U)) ') + silent (' num2str(numel(jc_silent_cell)) ')' ])

set(gca,'ytick',0:10:30,'xtick',600:100:1000,'xlim',[600 1000])
ylabel('firing rate (hz)');
xlabel('depth (um)')



