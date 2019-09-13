%% Raster
for i = 29
    spikes = squeeze(U{i}.R_ntk);
    motors = U{i}.meta.motorPosition;
    [~,sidx] = sort(motors);
    sidx = fliplr(sidx);
    
    figure(38);clf
    for g = 1:length(motors)
        subplot(1,2,1)
        spikeIdx = find(spikes(:,sidx(g)));
        hold on; scatter(spikeIdx,ones(numel(spikeIdx),1).*g,'k.')
    end
    set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])

    
    touchOn = [find(U{i}.S_ctk(9,:,:)==1)  ;find(U{i}.S_ctk(12,:,:)==1)];
    touchOff = [find(U{i}.S_ctk(10,:,:)==1)  ;find(U{i}.S_ctk(13,:,:)==1)];
    touch_matrix = nan(size(spikes));
    for g = 1:length(touchOn)
        touch_matrix(touchOn(g):touchOff(g)) = 1;
    end
    
    touch_matrix = touch_matrix(:,sidx)';
    subplot(1,2,2)
    pcolor(touch_matrix)
    set(gca,'ylim',[1 numel(motors)],'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])
    
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig1\';
    fn = 'example_raster.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
end
%% PSTH sorted by go and nogo
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

%% firing rate X depth of recording
figure(480);clf

jc_silent_cell = [766 819 895 631 776 815 910 871 844 902   941   840   888   748   732   940   686   944   950   933]; %Phils  from 902

scatter( cellfun(@(y) y.meta.depth,U),cellfun(@(y) mean(y.R_ntk(:))*1000, U),'k','filled')
hold on; scatter(jc_silent_cell,ones(length(jc_silent_cell),1)*-1,[],'c','filled');
hold on; plot([700 700],[0 30],'-.k')
hold on; plot([900 900],[0 30],'-.k')
title(['n = response (' num2str(numel(U)) ') + silent (' num2str(numel(jc_silent_cell)) ')' ])


set(gca,'ytick',0:10:30,'xtick',600:100:1000,'xlim',[600 1000])
ylabel('firing rate (hz)');
xlabel('depth (um)')



