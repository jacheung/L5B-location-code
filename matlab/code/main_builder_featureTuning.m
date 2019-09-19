%% TOUCH FEATURE TUNING
%% TOUCH
% defTouchResponse(U,.95,'off');
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
mod_idx_touch = cellfun(@(x) x.meta.touchProperties.mod_idx_relative,U(selectedCells));
%% POLE + IT FEATURES
pole_tuned = object_location_quantification(U,selectedCells,'pole','off'); 
angle_tuned = object_location_quantification(U,selectedCells,'angle','off'); 
phase_tuned = object_location_quantification(U,selectedCells,'phase','off'); 
amp_tuned = object_location_quantification(U,selectedCells,'amplitude','off'); 
mp_tuned = object_location_quantification(U,selectedCells,'midpoint','off'); 

mod_idx_pole = cellfun(@(x) x.calculations.mod_idx_relative,pole_tuned(selectedCells));
mod_idx_angle = cellfun(@(x) x.calculations.mod_idx_relative,angle_tuned(selectedCells));
mod_idx_phase = cellfun(@(x) x.calculations.mod_idx_relative,phase_tuned(selectedCells));
mod_idx_amp = cellfun(@(x) x.calculations.mod_idx_relative,amp_tuned(selectedCells));
mod_idx_midpoint = cellfun(@(x) x.calculations.mod_idx_relative,mp_tuned(selectedCells));
%% DYNAMIC FEATURES

dkappa_tuned = dynamic_touch_quantification(uberarray,selectedCells,'dkappa','on');
dtheta_tuned = dynamic_touch_quantification(uberarray,selectedCells,'dtheta','on');



    
%% ADAPTATION 
[adaptation] = adaptation_quantification(U,selectedCells,'off');
mod_idx_adaptation = cellfun(@(x) x.calculations.mod_idx_relative,adaptation(selectedCells));

%adaptation heat
population_heat = normalize_var(cell2mat(cellfun(@(x) x.lh,adaptation(selectedCells),'uniformoutput',0)')',0,1);
[~,idx] = sort(population_heat(1,:));
figure(81);clf
imagesc(population_heat(:,fliplr(idx)))
colormap bone
xlabel('cell number')
ylabel('touch order')
colorbar
set(gca,'ydir','reverse')

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'adaptation_map.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% PROTRACTION VS RETRACTION

figure(9);clf
for rec = 1:length(selectedCells)
    array = U{(selectedCells(rec))};
    
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
    dir_mod{(selectedCells(rec))}.mod_idx_relative = (pro_resp-ret_resp) ./  (pro_resp+ret_resp);
    
    pro_responses = touch_response(pro_touches)*1000; 
    ret_responses = touch_response(ret_touches)*1000;
    ret_error = std(ret_responses)./sqrt(numel(ret_touches));
    pro_error = std(pro_responses)./ sqrt(numel(pro_touches));
    
    hold on; errorbar(mean(ret_responses),mean(pro_responses),pro_error,pro_error,ret_error,ret_error,'ko')
   
end

hold on; plot([0 100],[0 100],'--k')
set(gca,'xlim',[0 100],'ylim',[0 100],'ytick',0:25:100,'xtick',0:25:100)
axis square
xlabel('retraction responses');ylabel('protraction responses')

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'pro_vs_ret_unity_scatter.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


mod_idx_directional = cellfun(@(x) x.mod_idx_relative,dir_mod(selectedCells));
%% PLOTTING

figure(8);clf
[~,idx] = sort(mod_idx_touch);
tune_map = [mod_idx_touch ; mod_idx_adaptation ;mod_idx_pole ; mod_idx_angle ; mod_idx_phase ; mod_idx_amp ; mod_idx_midpoint];
imagesc(tune_map(:,fliplr(idx)))
caxis([0 1])
reverse_bone = flipud(bone(500));
set(gca,'ytick',1:7,'yticklabel',{'touch','adaptation','pole','angle','phase','amp','midpoint'})
colormap(reverse_bone)
colorbar

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'touch_feature_map.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


figure(9);clf
imagesc(mod_idx_directional(fliplr(idx)))
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


%% WHISKING FEATURE TUNING

