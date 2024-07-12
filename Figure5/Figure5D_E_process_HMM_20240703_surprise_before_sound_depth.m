
%{
----------------------------------------------------------------------------
Determine the time window for task relevant neurons
%Start of task
%Sound on
%Sound off
%Before choice
%After choice (0sec)
%After choice (1sec)
%After choice (2sec)
%p = 0.001
%Each epoch, predict the prior, sensory and choice (integration)
----------------------------------------------------------------------------
%}

function Figure5D_E_process_HMM_20240703_surprise_before_sound_depth

%%% repeat
process_HMM_20240703_surprise_before_sound_depth('repeat_OFC_20230427',1);
process_HMM_20240703_surprise_before_sound_depth('repeat_AC_20230427',1);


% % alternate
process_HMM_20240703_surprise_before_sound_depth('zigzag_OFC_20230427',1);
process_HMM_20240703_surprise_before_sound_depth('zigzag_AC_20230427',1);

return



function process_HMM_20240703_surprise_before_sound_depth(folders, kaiseki_number)

switch nargin
    case 0
        hoge
    case 1
        kaiseki_number = 1;
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end

[analysis_dir,depth_def] = eval(folders);

brainarea = folders(1:10);

all_pre_sound = [];
all_p_pre_sound = [];
all_norm_spike = [];
all_p_surprise = [];
all_sound_correct_error = [];
all_p_sound_correct_error = [];
all_neuron_choice_category = [];

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category] = ...
        HMM_ephys_20240710_surprise_before_sound_depth(temp_dir, kaiseki_number,depth_def);

    all_pre_sound = [all_pre_sound; sound_neuron];
    all_p_pre_sound = [all_p_pre_sound; p_sound_neuron];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];

    all_sound_correct_error = [all_sound_correct_error; sound_correct_error];
    all_p_sound_correct_error = [all_p_sound_correct_error; p_sound_correct_error];

    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];
end
delete(gcp('nocreate'))

disp([length(all_pre_sound),length(all_norm_spike),length(all_p_surprise)])

if contains(brainarea,'repeat')
    left_neuron = find(all_pre_sound(:,6) < 0);
    right_neuron = find(all_pre_sound(:,6) >= 0);
    sig_sound1 = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between "PREVIOUS" left and right
    sig_sound = sig_sound1;
else
    left_neuron = find(all_pre_sound(:,6) >= 0);
    right_neuron = find(all_pre_sound(:,6) < 0);
    sig_sound1 = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between "PREVIOUS" left and right
    sig_sound = sig_sound1;
end
right_sig_neuron = intersect(right_neuron, sig_sound);
left_sig_neuron = intersect(left_neuron, sig_sound);
non_sig_neuron = setdiff(1:length(all_pre_sound), sig_sound);
check_neuron = length(right_sig_neuron) + length(left_sig_neuron) + length(non_sig_neuron);
if check_neuron ~= length(all_pre_sound)
    hoge
else
   disp([length(left_sig_neuron),length(right_sig_neuron),length(non_sig_neuron)])
end
size(all_sound_correct_error)

% %Plot the tone index, previous correct or previous error
plot_index_abs_index_sig(all_pre_sound(:,6), -all_pre_sound(:,7), right_sig_neuron,left_sig_neuron,non_sig_neuron,brainarea);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, median_value,rho,pval,b_regress,p_regress] = plot_index_abs_index_sig(index1, index2, right_neuron,left_neuron,non_sig_neuron,brainarea)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_neuron = union(left_neuron, right_neuron);

figure
subplot(2,2,1)
plot([-1 1],[0 0],'k:')
hold on
plot([0 0],[-1 1],'k:')
hold on
plot([-1 1],[-1 1],'k:')
hold on

plot(index1(non_sig_neuron), index2(non_sig_neuron), 'k.')
hold on
plot(index1(left_neuron), index2(left_neuron), 'b.')
hold on
plot(index1(right_neuron), index2(right_neuron), 'r.')

[rho(1),pval(1)] = corr(index1,index2);
[rho(2),pval(2)] = corr(index1(sig_neuron),index2(sig_neuron));

%regression

[b_regress(1),~,stats] = glmfit(index1,index2,'normal','link','identity','Constant','off');
p_regress(1) = stats.p;
[b_regress(2),~,stats] = glmfit(index1(sig_neuron),index2(sig_neuron),'normal','link','identity','Constant','off');
p_regress(2) = stats.p;

subplot(2,2,2)
plot([0 1],[0 1],'k:')
hold on

plot(abs(index1(non_sig_neuron)), abs(index2(non_sig_neuron)), 'k.')
hold on
plot(abs(index1(left_neuron)), abs(index2(left_neuron)), 'b.')
hold on
plot(abs(index1(right_neuron)), abs(index2(right_neuron)), 'r.')

p(1) = signrank(abs(index1), abs(index2));
p(2) = signrank(abs(index1(sig_neuron)), abs(index2(sig_neuron)));
median_value = [median(abs(index1)), median(abs(index2))];

subplot(2,2,3)
plot([-1 1],[0 0],'k:')
hold on
plot([0 0],[-1 1],'k:')
hold on
plot([-1 1],[-1 1],'k:')
hold on
plot(index1(left_neuron), index2(left_neuron), 'b.')
hold on
plot(index1(right_neuron), index2(right_neuron), 'r.')
hold on 
plot([-1,1],[-b_regress(2),b_regress(2)],'k')

subplot(2,2,4)
plot([0 1],[0 1],'k:')
hold on
plot(abs(index1(left_neuron)), abs(index2(left_neuron)), 'b.')
hold on
plot(abs(index1(right_neuron)), abs(index2(right_neuron)), 'r.')

close 

figure
plot([-1 1],[0 0],'k:')
hold on
plot([0 0],[-1 1],'k:')
hold on
plot([-1 1],[-1 1],'k:')
hold on
plot(index1(left_neuron), index2(left_neuron), 'b.')
hold on
plot(index1(right_neuron), index2(right_neuron), 'r.')
hold on 
plot([-1,1],[-b_regress(2),b_regress(2)],'k')
text(0.5,-0.7, ['p = ';string(p_regress(2))])
box('off')
title(brainarea(8:10))
set(gca,'xlim',[-1, 1],'ylim',[-1, 1],'fontname','Arial')
set(gcf,'Position',[584,652,295,263])


return


function [pre_sound, p_pre_sound, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category] = ...
    HMM_ephys_20240710_surprise_before_sound_depth(pathname, kaiseki_number,depth_def)

switch nargin
    case 0
        pathname = pwd;
    case 3
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

% [filename1, pathname1,findex]=uigetfile('*.*','frame file');
% filename1 = [pathname1,filename1];
% load(filename1)
temp = dir('HMM_spike_count_neurons_20240420*');
if length(temp) ~= 1
    hoge
end
load(temp.name);
%neuron_index p_index

temp = dir('sig_HMM_neurons_20230310*');
if length(temp) ~= 1    
    hoge
end
load(temp.name);
%p_task: around sound
%p_task2: around choice

temp = dir('depth_spike_20230427*');
if length(temp) == 1
    load(temp.name);
    %spike_depth def_depth length_neuron
    if size(p_task,1) ~= length_neuron
        disp([length(p_task) length_neuron])
        hoge
    end
    
    if depth_def == 1
        depth_neuron = find(spike_depth <= def_depth(1));
    else
        depth_neuron = find(spike_depth > def_depth(1) & spike_depth <= def_depth(2));
    end
elseif isempty(temp)
    depth_neuron = 1:size(p_task,1); %Use all the neurons
    hoge
else
    
    hoge
end

% %Get sig neurons
new_p_thre = 10;

sig_before_sound = get_sig_neuron_time_window(p_task,new_p_thre,10:15); %600ms before sound
sig_during_sound = get_sig_neuron_time_window(p_task,new_p_thre,16:21); %During sound
sig_during_choice = get_sig_neuron_time_window(p_task2,new_p_thre,6:15); %During choice 1000ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if kaiseki_number == 1 %Before sound
    use_sig_neuron = sig_before_sound;
    use_neuron_index = neuron_index.before_sound;
    use_p_index = p_index.before_sound;
    use_norm_spike = norm_spike_group.before_sound;
    use_p_surprise = p_surprise.before_sound;

elseif kaiseki_number == 2 %during sound
    use_sig_neuron = sig_during_sound;
    use_neuron_index = neuron_index.all_sound;
    use_p_index = p_index.all_sound;
    use_norm_spike = norm_spike_group.all_sound;
    use_p_surprise = p_surprise.all_sound;
    
elseif kaiseki_number == 3 %during choice %1000ms
    use_sig_neuron = sig_during_choice;
    use_neuron_index = neuron_index.choice2;
    use_p_index = p_index.choice2;
    use_norm_spike = norm_spike_group.choice2;
    use_p_surprise = p_surprise.choice2;
    
else
    hopge
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADD depth definition
use_sig_neuron = intersect(use_sig_neuron,depth_neuron);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pre_sound(:,1) = use_neuron_index(use_sig_neuron,1); %current_choice
pre_sound(:,2) = use_neuron_index(use_sig_neuron,6); %pre sound
pre_sound(:,3) = use_neuron_index(use_sig_neuron,7); %pre choice
pre_sound(:,4) = use_neuron_index(use_sig_neuron,10); %pre SoundCorrect, current left
pre_sound(:,5) = use_neuron_index(use_sig_neuron,11); %pre SoundCorrect, current right
pre_sound(:,6) = use_neuron_index(use_sig_neuron,12); %pre SoundCorrect
pre_sound(:,7) = use_neuron_index(use_sig_neuron,13); %pre SoundError

p_pre_sound(:,1) = use_p_index(use_sig_neuron,1); %current_choice
p_pre_sound(:,2) = use_p_index(use_sig_neuron,6); %pre sound
p_pre_sound(:,3) = use_p_index(use_sig_neuron,7); %pre choice
p_pre_sound(:,4) = use_p_index(use_sig_neuron,10); %pre SoundCorrect, current left
p_pre_sound(:,5) = use_p_index(use_sig_neuron,11); %pre SoundCorrect, current right
p_pre_sound(:,6) = use_p_index(use_sig_neuron,12); %pre SoundCorrect
p_pre_sound(:,7) = use_p_index(use_sig_neuron,13); %pre SoundError

sound_correct_error(:,1) = use_neuron_index(use_sig_neuron,4); %Correct Low high
sound_correct_error(:,2) = use_neuron_index(use_sig_neuron,5); %Error Low high
p_sound_correct_error(:,1) = use_p_index(use_sig_neuron,4);
p_sound_correct_error(:,2) = use_p_index(use_sig_neuron,5);

norm_spike = use_norm_spike(use_sig_neuron,:);
p_surprise = use_p_surprise(use_sig_neuron,:);

temp = norm_spike(:,2) - norm_spike(:,1); %right-left
neuron_choice_category = zeros(length(temp),1);
temp = find(temp > 0);
neuron_choice_category(temp) = 1; %right=1, left=0

%%%
%Remove NAN
test = [pre_sound, norm_spike, sound_correct_error];
test = mean(test,2);
test = isnan(test);
test_nan = find(test == 1);
pre_sound(test_nan,:) = [];
p_pre_sound(test_nan,:) = [];
norm_spike(test_nan,:) = [];
p_surprise(test_nan,:) = [];
sound_correct_error(test_nan,:) = [];
p_sound_correct_error(test_nan,:) = [];
neuron_choice_category(test_nan,:) = [];


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_neuron = get_sig_neuron_time_window(p_task,new_p_thre,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_p_sound = min(p_task(:,window),[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return
