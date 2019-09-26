%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions
%U = defTouchResponse(U,.95,'on');
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
pole_tuned = object_location_quantification(U,selectedCells,'pole','off'); %for old see object_location_v1.0

for i = 29
    %% raster
    motors = normalize_var(U{i}.meta.motorPosition,1,-1);
    spikes = squeeze(U{i}.R_ntk);
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
    
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
    fn = 'example_raster.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])

    %% chunked psth 
    sorted_spike_mat = spikes(:,sidx);
    chunk_one = sidx(1:round(numel(sidx)/3));
    chunk_two = sidx(round(numel(sidx)/3)+1: round(numel(sidx)/3)+1 + round(numel(sidx)/3));
    chunk_three = setdiff(sidx,[chunk_one chunk_two]);
    all_chunks = {chunk_three, chunk_two,chunk_one};
    
    figure(49);clf
    colors = [.3 .3 .3];
    for b = 1:length(all_chunks)
        hold on;plot(smooth(mean(spikes(:,all_chunks{b}),2)*1000,100),'color',colors.*b)
        set(gca,'ylim',[0 40],'xtick',0:1000:4000,'ytick',0:10:40)
        title(['motor pos = ' num2str(mean(motors(all_chunks{b})))])
    end
    title(['far:close ' num2str(cellfun(@(x) mean(motors(x)),all_chunks))])
    
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
    fn = 'example_psth.eps';
    export_fig([saveDir, fn], '-depsc ', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
        
end
%% touch psth by quartiles of far, close and near 
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
pole_tuned = object_location_quantification(U,selectedCells,'pole','on'); %for old see object_location_v1.0
tuned_units = find(cellfun(@(x) x.is_tuned==1,pole_tuned));
pdm = preDecisionTouchMat(U); 

touch_window = -25:50;
chunks = 3;
units_to_plot = [6 ,30, 9];
figure(20);clf
for g = 1:length(units_to_plot)
    selected_unit = tuned_units(units_to_plot(g)); 
    
    motors = normalize_var(U{selected_unit}.meta.motorPosition,1,-1);
    leftover = mod(length(motors),chunks);
    new_motors = datasample(motors,numel(motors)-leftover,'Replace',false);
    [s_motors,sidx] = sort(new_motors);
    all_chunks = reshape(s_motors,numel(new_motors)./chunks,[]); 

    
    [tVar] = atTouch_sorter(U{selected_unit},touch_window,pdm{selected_unit});
    touch_motors = normalize_var(tVar.allTouches.S_ctk(:,end),-1,1); 
    
    chunked_idx = cell(1,size(all_chunks,2));
    chunked_responses = cell(1,size(all_chunks,2));
    figure(20);
    colors = [.33 .33 .33];
    subplot(3,1,g)
    for b = 1:size(all_chunks,2)
        [rsmall,rbig] = bounds(all_chunks(:,b));
        chunked_idx = intersect(find(touch_motors<rbig),find(touch_motors>rsmall)); 
        chunked_responses{b} = tVar.allTouches.R_ntk(chunked_idx,:);
        
        hold on; plot(touch_window,smooth(nanmean(chunked_responses{b}).*1000,10),'color',colors.*b)

    end

    set(gca,'xtick',-25:25:50,'xlim',[-25 50])
end

    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
    fn = 'touch_psth_by_location.eps';
    export_fig([saveDir, fn], '-depsc ', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
    
%% modulation width
tuned_units = find(cellfun(@(x) x.is_tuned==1,pole_tuned));

peak_response = cellfun(@(x) x.calculations.tune_peak,pole_tuned(tuned_units));
interp_norm_y = cell(1,length(peak_response));
for g = 1:length(peak_response)
    sr = pole_tuned{tuned_units(g)}.stim_response;
    centered_x = sr.values(:,1) - peak_response(g) ;
    norm_y = normalize_var(sr.values(:,2),0,1);
    
    interp_centered_x = -2:.1:2;
    raw_x = round(centered_x,2);
    [~,idx] = unique(raw_x);
    interp_norm_y{g} = interp1(raw_x(idx),norm_y(idx),interp_centered_x);
    figure(30);
    hold on; plot(interp_centered_x,interp_norm_y{g},'color',[.8 .8 .8])
end

pop_mean = nanmean(cell2mat(interp_norm_y'));
pop_sem = nanstd(cell2mat(interp_norm_y')) ./ sqrt(sum(~isnan(cell2mat(interp_norm_y'))));
hold on; shadedErrorBar(interp_centered_x,pop_mean,pop_sem,'k')

set(gca,'xlim',[-2 2],'xdir','reverse','ytick',0:.5:1,'xtick',-2:1:2)
xlabel('distance from peak (mm)')
axis square

        figure(30);
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
    fn = 'all_modulation_width.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])