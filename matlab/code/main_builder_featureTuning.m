%% TOUCH FEATURE TUNING
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% STIMULUS CORRELATION
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
preDecisionTouches = preDecisionTouchMat(U);
corr_mat = zeros(7, 7);
for rec = 1:length(touchCells)
    %stimulus and response variables definitions
    array = U{touchCells(rec)};
    [tVar] = atTouch_sorter(array,-25:50,preDecisionTouches{touchCells(rec)});
    it_features = tVar.allTouches.S_ctk(:,[7 1 5 3 4]);
    dt_features = tVar.allTouches.dtS_ctk(:,2:3);
    all_features = [it_features dt_features];
    all_features(any(isnan(all_features),2),:) = [];
    
    corr_mat = corr_mat + corr(all_features);
end

figure(430);clf
imagesc(abs(corr_mat./numel(touchCells)))
set(gca,'xtick',1:7,'xticklabel',{'pole','angle','phase','amp','midpoint','dkappa','dtheta'})
axis square
colormap gray
caxis([0 1])

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'stimulus_correlation.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% TOUCH
% U = defTouchResponse(U,.95,'off');
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
mod_idx_touch = cellfun(@(x) x.meta.touchProperties.mod_idx_relative,U(touchCells));

%% POLE + IT FEATURES
pole_tuned = object_location_quantification(U,touchCells,'pole','on');
angle_tuned = object_location_quantification(U,touchCells,'angle','off');
phase_tuned = object_location_quantification(U,touchCells,'phase','off');
amp_tuned = object_location_quantification(U,touchCells,'amplitude','off');
mp_tuned = object_location_quantification(U,touchCells,'midpoint','off');

mod_idx_pole = cellfun(@(x) x.calculations.mod_idx_relative,pole_tuned(touchCells));
mod_idx_angle = cellfun(@(x) x.calculations.mod_idx_relative,angle_tuned(touchCells));
mod_idx_phase = cellfun(@(x) x.calculations.mod_idx_relative,phase_tuned(touchCells));
mod_idx_amp = cellfun(@(x) x.calculations.mod_idx_relative,amp_tuned(touchCells));
mod_idx_midpoint = cellfun(@(x) x.calculations.mod_idx_relative,mp_tuned(touchCells));

peak_idx_pole = cellfun(@(x) x.calculations.tune_peak,pole_tuned(touchCells)) * -1;
peak_idx_angle = cellfun(@(x) x.calculations.tune_peak,angle_tuned(touchCells));
peak_idx_phase = cellfun(@(x) x.calculations.tune_peak,phase_tuned(touchCells));
peak_idx_amp = cellfun(@(x) x.calculations.tune_peak,amp_tuned(touchCells));
peak_idx_midpoint = cellfun(@(x) x.calculations.tune_peak,mp_tuned(touchCells));

% HSV PLOTTING

% [~,sort_idx] = sort(peak_idx_pole);
% 
% peaks = {peak_idx_pole,peak_idx_angle,peak_idx_phase,peak_idx_amp,peak_idx_midpoint};
% mod_idx = {mod_idx_pole,mod_idx_angle,mod_idx_phase,mod_idx_amp,mod_idx_midpoint};
% final_image = nan(numel(peaks),numel(selectedCells),3); 
% for b = 1:numel(peaks)
%     hues = normalize_var(peaks{b}(sort_idx),.7,1);
%     hsv = [hues' mod_idx{b}(sort_idx)' ones(numel(hues),1)];
%     rgb_values = hsv2rgb(hsv);
%     final_image(b,:,1) = rgb_values(:,1);
%     final_image(b,:,2) = rgb_values(:,2);
%     final_image(b,:,3) = rgb_values(:,3);
% end
% 
% figure(580);clf
% imshow(final_image)
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = 'hsv_touch.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

%% time of peak response and tuning
peak_idx_touch = cellfun(@(x) x.meta.touchProperties.peak_index,U(touchCells));  
peak_idx_pole = cellfun(@(x) x.calculations.tune_peak,pole_tuned(touchCells)) * -1;

%% DYNAMIC FEATURES

dkappa_tuned = dynamic_touch_quantification(U,touchCells,'dkappa','off');
dtheta_tuned = dynamic_touch_quantification(U,touchCells,'dtheta','off');

mod_idx_dk = cellfun(@(x) x.calculations.mod_idx_relative,dkappa_tuned(touchCells));
mod_idx_dt = cellfun(@(x) x.calculations.mod_idx_relative,dtheta_tuned(touchCells));


%% ADAPTATION
[adaptation] = adaptation_quantification(U,touchCells,'off');
mod_idx_adaptation = cellfun(@(x) x.calculations.mod_idx_relative,adaptation(touchCells));


tunedCells = find(cellfun(@(x) x.is_tuned==1,pole_tuned));
nontunedCells = setdiff(touchCells,tunedCells);
%adaptation heat
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

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'adaptation_map.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% PROTRACTION VS RETRACTION
alpha_values = 0.05;
figure(9);clf
for rec = 1:length(touchCells)
    array = U{(touchCells(rec))};
    
    spks = squeeze(array.R_ntk(:,:,:));
    response_window = array.meta.touchProperties.responseWindow(1):array.meta.touchProperties.responseWindow(2);
    touch_times = [find(array.S_ctk(9,:,:)==1) ; find(array.S_ctk(12,:,:)==1)];
    phase = squeeze(array.S_ctk(5,:,:));
    velocity = squeeze(array.S_ctk(2,:,:));
    pt_vel = mean(velocity(touch_times - [-5:-1]),2);
    touch_phase = phase(touch_times);
    
    touch_response = mean(spks(touch_times + response_window),2);
    pro_touches = intersect(find(touch_phase<0),find(pt_vel>0));
    ret_touches = intersect(find(touch_phase>0),find(pt_vel<0));
    
    pro_resp = max(touch_response(pro_touches));
    ret_resp = max(touch_response(ret_touches));
    dir_mod{(touchCells(rec))}.mod_idx_relative = (pro_resp-ret_resp) ./  (pro_resp+ret_resp);
    
    pro_responses = touch_response(pro_touches)*1000;
    ret_responses = touch_response(ret_touches)*1000;
    [~,p] = ttest2(pro_responses,ret_responses);
    ret_error = std(ret_responses)./sqrt(numel(ret_touches));
    pro_error = std(pro_responses)./ sqrt(numel(pro_touches));
    hold on; errorbar(mean(ret_responses),mean(pro_responses),pro_error,pro_error,ret_error,ret_error,'k','CapSize',0)
    if p < alpha_values
        if mean(pro_responses)>mean(ret_responses)
            hold on; scatter(mean(ret_responses),mean(pro_responses),'filled','markerfacecolor','r')
        else
            hold on; scatter(mean(ret_responses),mean(pro_responses),'filled','markerfacecolor','b')
        end
    else
        hold on; scatter(mean(ret_responses),mean(pro_responses),'filled','markerfacecolor',[.8 .8 .8])
    end
    
end

hold on; plot([0 100],[0 100],'--k')
set(gca,'xlim',[0 100],'ylim',[0 100],'ytick',0:25:100,'xtick',0:25:100)
axis square
xlabel('retraction responses');ylabel('protraction responses')

% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = 'pro_vs_ret_unity_scatter.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])


mod_idx_directional = cellfun(@(x) x.mod_idx_relative,dir_mod(touchCells));
%% PLOTTING
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
tunedCells = find(cellfun(@(x) x.is_tuned==1,pole_tuned));
nontunedCells = setdiff(touchCells,tunedCells);


figure(8);clf
[~,idx] = sort(mod_idx_pole);
tune_map_1 = [mod_idx_touch ; mod_idx_adaptation];
tune_map_it = [mod_idx_pole ; mod_idx_angle ; mod_idx_phase ; mod_idx_amp ; mod_idx_midpoint];
tune_map_dt = [mod_idx_dk ; mod_idx_dt];
subplot(3,1,1);
imagesc(tune_map_1(:,(idx)))
set(gca,'ytick',1:2,'yticklabel',{'touch','adaptation'})
caxis([0 1])
subplot(3,1,2);
imagesc(tune_map_it(:,(idx)))
set(gca,'ytick',1:6,'yticklabel',{'pole','angle','phase','amp','midpoint'})
caxis([0 1])
subplot(3,1,3);
imagesc(tune_map_dt(:,(idx)))
set(gca,'ytick',1:2,'yticklabel',{'max dkappa','max dtheta'})
caxis([0 1])
% reverse_bone = flipud(bone(500));
colormap(gray)
colorbar
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'touch_feature_map_all.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

figure(9);clf
imagesc(mod_idx_directional((idx)))
stretch_resolution = 500 ;
stretch_map = redbluecmap;
new_map = nan(stretch_resolution,3);
for j = 1:3
    new_map(:,j) = interp1(linspace(1,stretch_resolution,length(stretch_map)),stretch_map(:,j),1:stretch_resolution);
end
colormap(new_map)
caxis([-1 1])
colorbar
set(gca,'ytick',1:5,'yticklabel',{'direction'})

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'touch_feature_map_direction.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

% Feature correlation
mod_idx_correlation = corr([tune_map_it ; tune_map_dt ; tune_map_1]');
figure(143);clf
imagesc(abs(mod_idx_correlation))
colormap(gray)
colorbar
caxis([0 1])
axis square

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'mod_depth_correlation.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% TUNED UNITS vs UNTUNED

touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
tunedCells = find(cellfun(@(x) x.is_tuned==1,pole_tuned));
nontunedCells = setdiff(touchCells,tunedCells);

[~,touchtuneidx] = intersect(touchCells,tunedCells);
[~,touchnontuneidx] = intersect(touchCells,nontunedCells);
%tuned
[~,tune_idx] = sort(mod_idx_pole(touchtuneidx));
tune_map_1 = [mod_idx_touch(touchtuneidx) ; mod_idx_adaptation(touchtuneidx)];
tune_map_it = [mod_idx_pole(touchtuneidx) ; mod_idx_angle(touchtuneidx) ; mod_idx_phase(touchtuneidx) ; mod_idx_amp(touchtuneidx) ; mod_idx_midpoint(touchtuneidx)];
tune_map_dt = [mod_idx_dk(touchtuneidx) ; mod_idx_dt(touchtuneidx)];

figure(81);clf
subplot(3,2,1);
imagesc(tune_map_1(:,(tune_idx)))
set(gca,'ytick',1:2,'yticklabel',{'touch','adaptation'})
title('tuned')
caxis([0 1])
subplot(3,2,3);
imagesc(tune_map_it(:,(tune_idx)))
set(gca,'ytick',1:6,'yticklabel',{'pole','angle','phase','amp','midpoint'})
caxis([0 1])
subplot(3,2,5);
imagesc(tune_map_dt(:,(tune_idx)))
set(gca,'ytick',1:2,'yticklabel',{'max dkappa','max dtheta'})
caxis([0 1])

%nontuned
[~,non_tune_idx] = sort(mod_idx_pole(touchnontuneidx));
tune_map_1 = [mod_idx_touch(touchnontuneidx) ; mod_idx_adaptation(touchnontuneidx)];
tune_map_it = [mod_idx_pole(touchnontuneidx) ; mod_idx_angle(touchnontuneidx) ; mod_idx_phase(touchnontuneidx) ; mod_idx_amp(touchnontuneidx) ; mod_idx_midpoint(touchnontuneidx)];
tune_map_dt = [mod_idx_dk(touchnontuneidx) ; mod_idx_dt(touchnontuneidx)];

subplot(3,2,2);
imagesc(tune_map_1(:,(non_tune_idx)))
set(gca,'ytick',1:2,'yticklabel',{'touch','adaptation'})
title('nontuned')
caxis([0 1])
subplot(3,2,4);
imagesc(tune_map_it(:,(non_tune_idx)))
set(gca,'ytick',1:6,'yticklabel',{'pole','angle','phase','amp','midpoint'})
caxis([0 1])
subplot(3,2,6);
imagesc(tune_map_dt(:,(non_tune_idx)))
set(gca,'ytick',1:2,'yticklabel',{'max dkappa','max dtheta'})
caxis([0 1])

colormap(gray)
colorbar


saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'touch_feature_map_split.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%PROTRACTION RETRACTION TUNING
figure(91);clf
tune_mod_idx_directional = mod_idx_directional(tune_idx);
nontune_mod_idx_directional = mod_idx_directional(non_tune_idx);
subplot(1,2,1);
imagesc(tune_mod_idx_directional(tune_idx))
caxis([-1 1])
subplot(1,2,2);
imagesc(nontune_mod_idx_directional(non_tune_idx))
caxis([-1 1])

stretch_resolution = 500 ;
stretch_map = redbluecmap;
new_map = nan(stretch_resolution,3);
for j = 1:3
    new_map(:,j) = interp1(linspace(1,stretch_resolution,length(stretch_map)),stretch_map(:,j),1:stretch_resolution);
end
colormap(new_map)
caxis([-1 1])
colorbar
set(gca,'ytick',1:5,'yticklabel',{'direction'})

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'touch_feature_map_direction_split.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])






%% WHISKING FEATURE TUNING
touchCells = 1:length(U)
pole_whisk = whisking_location_quantification(U,touchCells,'pole','off');
angle_whisk = whisking_location_quantification(U,touchCells,'angle','off');
midpoint_whisk = whisking_location_quantification(U,touchCells,'midpoint','off');
amp_whisk = whisking_location_quantification(U,touchCells,'amplitude','off');
phase_whisk = whisking_location_quantification(U,touchCells,'phase','off');

w_mod_idx_pole = cellfun(@(x) x.calculations.mod_idx_relative,pole_whisk(touchCells));
w_mod_idx_angle = cellfun(@(x) x.calculations.mod_idx_relative,angle_whisk(touchCells));
w_mod_idx_phase = cellfun(@(x) x.calculations.mod_idx_relative,phase_whisk(touchCells));
w_mod_idx_amp = cellfun(@(x) x.calculations.mod_idx_relative,amp_whisk(touchCells));
w_mod_idx_midpoint = cellfun(@(x) x.calculations.mod_idx_relative,midpoint_whisk(touchCells));

w_peak_idx_pole = cellfun(@(x) x.calculations.tune_peak,pole_whisk(touchCells)) * -1;
w_peak_idx_angle = cellfun(@(x) x.calculations.tune_peak,angle_whisk(touchCells));
w_peak_idx_phase = cellfun(@(x) x.calculations.tune_peak,phase_whisk(touchCells));
w_peak_idx_amp = cellfun(@(x) x.calculations.tune_peak,amp_whisk(touchCells));
w_peak_idx_midpoint = cellfun(@(x) x.calculations.tune_peak,midpoint_whisk(touchCells));


% HSV PLOTTING WHISK
[~,sort_idx] = sort(w_peak_idx_pole);

peaks = {w_peak_idx_pole,w_peak_idx_angle,w_peak_idx_phase,w_peak_idx_amp,w_peak_idx_midpoint};
mod_idx = {w_mod_idx_pole,w_mod_idx_angle,w_mod_idx_phase,w_mod_idx_amp,w_mod_idx_midpoint};
final_image = nan(numel(peaks),numel(touchCells),3); 
for b = 1:numel(peaks)
    hues = normalize_var(peaks{b}(sort_idx),.7,1);
    hsv = [hues' mod_idx{b}(sort_idx)' ones(numel(hues),1)];
    rgb_values = hsv2rgb(hsv);
    final_image(b,:,1) = rgb_values(:,1);
    final_image(b,:,2) = rgb_values(:,2);
    final_image(b,:,3) = rgb_values(:,3);
end
% 
% figure(580);clf
% imshow(final_image)
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = 'hsv_whisk.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

% HEAT GRAY WHISK 
[~,idx] = sort(w_mod_idx_pole);
tune_map = [w_mod_idx_pole ; w_mod_idx_angle ; w_mod_idx_phase ; w_mod_idx_amp ; w_mod_idx_midpoint];
figure(19);clf
imagesc(tune_map(:,(idx)))
set(gca,'ytick',1:6,'yticklabel',{'pole','angle','phase','amp','midpoint'})
caxis([0 1])
colormap(gray)
colorbar
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'gray_whisk.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

figure(144);clf
imagesc(abs(corr(tune_map')))
colormap(gray)
colorbar
caxis([0 1])
axis square
set(gca,'xtick',[],'ytick',[])
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'mod_depth_correlation_whisk.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])





