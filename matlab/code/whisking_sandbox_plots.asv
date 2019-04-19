%% whiisking vs quiet
for i = 1:length(U)
    
    [objmask]= maskBuilder(U{i});
    spikes = squeeze(U{i}.R_ntk(:,:,:)); 
    filteredSpikes = spikes.*objmask.touch;
    whiskSpikes = filteredSpikes .* objmask.whisking;
    quietSpikes = filteredSpikes .* objmask.quiet; 
        
    [~,p(i)] = ttest2(quietSpikes(:),whiskSpikes(:));

    nwfr(i) = nanmean(quietSpikes(:))*1000;
    wfr(i) = nanmean(whiskSpikes(:))*1000;
    
end

whiskTunedcells = find(p<.01);
unsig = find(p>.01);

figure(9);subplot(1,3,1)
nwc = nan(length(U),1);
wc = nan(length(U),1);

for g = 1:length(whiskTunedcells)
    curr = whiskTunedcells(g);
    if nwfr(curr) > wfr(curr)
        hold on;scatter(nwfr(curr),wfr(curr),[],'filled','r')
        nwc(curr) = 1;
    elseif nwfr(curr)<wfr(curr)
        hold on;scatter(nwfr(curr),wfr(curr),[],'filled','b')
        wc(curr) = 1;
    end
end
scatter(nwfr(unsig),wfr(unsig),[],'filled','markerfacecolor',[.7 .7 .7])
set(gca,'xlim',[0 30],'ylim',[0 30],'xtick',0:10:30,'ytick',0:10:30)
hold on; plot([0 30],[0 30],'-.k')
xlabel('no whisk spks/s');ylabel('whisk spks/s')
axis square

vals = [nwfr' wfr'];
nwcells = vals(find(nwc==1),:);
wcells = vals(find(wc==1),:);


figure(8);clf;
subplot(2,1,1);histogram(wcells(:,1),0:2.5:30,'facealpha',1,'facecolor','b')
hold on; histogram(nwcells(:,1),0:2.5:30,'facealpha',1,'facecolor','r')
set(gca,'ylim',[0 20],'ytick',0:20:20)
ylabel('nw spks')

subplot(2,1,2);histogram(wcells(:,2),0:2.5:30,'facealpha',1,'facecolor','b')
hold on; histogram(nwcells(:,2),0:2.5:30,'facealpha',1,'facecolor','r')
set(gca,'ylim',[0 20],'ytick',0:20:20)
ylabel('w spks')


%prop tuned blocking out touch periods.
length(whiskTunedcells)./length(U)


%% WHISKING ON CELL ANALYSIS PRO/RET analysis
wccells = find(wc==1);
rVSp = zeros(length(wccells),2);

pthresh = .01;

figure(9);subplot(1,3,2)
for d = 1:length(wccells)
    curr = U{wccells(d)};
    masks = assist_touchmasks(curr);
    
    amps = squeeze(curr.S_ctk(3,:,:)).*masks.touch;
    phase = squeeze(curr.S_ctk(5,:,:)).*masks.touch;
    spks = squeeze(curr.R_ntk(1,:,:)).*masks.touch;
    
    retWidx = intersect(find(amps>5),find(phase>0));
    proWidx = intersect(find(amps>5),find(phase<0));
    
    [~,rppval(d)] = ttest2(spks(retWidx),spks(proWidx));
    
    rVSp(d,:) = [mean(spks(retWidx))*1000 mean(spks(proWidx))*1000];
    
    [~,rORp]  = max(rVSp(d,:));
    
    if rppval(d)<=pthresh
        if rORp ==2
            hold on; scatter(rVSp(d,1),rVSp(d,2),[],'filled','c')
        elseif rORp ==1
            hold on; scatter(rVSp(d,1),rVSp(d,2),[],'filled','m')
        end
    elseif rppval(d)>pthresh
        hold on; scatter(rVSp(d,1),rVSp(d,2),[],'filled','markerfacecolor',[.7 .7 .7])
    end
    
end

hold on; plot([0 30],[0 30],'-.k')
axis square
xlabel('retraction spks/s'); ylabel('protraction spks/s')
set(gca,'xtick',0:10:30,'ytick',0:10:30)

propWhiskDtuned = sum(rppval<=pthresh)./numel(rppval);
dirCells = (rppval<=pthresh);
rVSp(dirCells,:) 
%% peak (blue) of whisking vs trough (cyan) of whisking tuning 
viewWindow = [-25:50];
amplitudeThreshold = 5; 
[rc] = numSubplots(length(U));
 figure(21);clf
for i = 1:length(U)
    
    [objmask]= maskBuilder(U{i});
    P = findPeaksTroughs(U{i},amplitudeThreshold); 
    spikes = squeeze(U{i}.R_ntk(:,:,:)); 
    filteredSpikes = spikes.*objmask.touch;
    whiskSpikes = filteredSpikes .* objmask.whisking;
    quietSpikes = filteredSpikes .* objmask.quiet; 
    
    peakWindow = find(P.peakIdx==1) + viewWindow;
    troughWindow = find(P.troughIdx==1) + viewWindow;
    

    %trimming for edge cases 
    peakWindow(sum(peakWindow > U{i}.t * U{i}.k,2)>0,:) = [];
    peakWindow(sum(peakWindow <= 0,2)>0,:) = [];
    troughWindow(sum(troughWindow > U{i}.t * U{i}.k,2)>0,:) = [];
    troughWindow(sum(troughWindow <= 0)>0,:) = [];
    
    %building spike mats
    peakSpikes = filteredSpikes(peakWindow);
    troughSpikes = filteredSpikes(troughWindow);
    
    %filtering spike mats with nan values which are touch windows
    %(don't want to just nanmean because it'll give false view of decrease firing) 
    peakSpikes(sum(isnan(peakSpikes),2)>0,:) = [];
    troughSpikes(sum(isnan(troughSpikes),2)>0,:) = [];
    
%     figure(21);subplot(rc(1),rc(2),i)
%     plot(viewWindow,nanmean(peakSpikes,1)*1000,'b')
%     hold on; plot(viewWindow,nanmean(troughSpikes,1)*1000,'c')
%     
% % %     shadedErrorBar(viewWindow,nanmean(peakSpikes,1)*1000,nanstd(peakSpikes,1)*1000./sqrt(size(peakWindow,1)),'b')
% % %     hold on;shadedErrorBar(viewWindow,nanmean(troughSpikes,1)*1000,nanstd(troughSpikes,1)*1000./sqrt(size(troughWindow,1)),'c')
%     hold on; plot([0 0],[0 max(nanmean(peakSpikes,1)*1000)],'-.k')
% %     

    peaks(:,i) = mean(peakSpikes);
    troughs(:,i) = mean(troughSpikes);

end


figure(9);subplot(1,3,3)
% plot(viewWindow,normalize_var(peaks,0,1),'color',[.8 .8 .8])
hold on; shadedErrorBar(viewWindow,nanmean(peaks(:,wccells),2)*1000,nanstd(peaks,[],2)*1000,'b')
hold on; shadedErrorBar(viewWindow,nanmean(troughs(:,wccells),2)*1000,nanstd(troughs,[],2)*1000,'c')
set(gca,'xlim',[min(viewWindow) max(viewWindow)],'xtick',-25:25:50)
hold on; plot([0 0],[0 25],'-.k')
title('spks/s from whisking peak(blue) and trough(cyan)') 
axis square
xlabel('time from peak or trough') 
ylabel('spks/s') 