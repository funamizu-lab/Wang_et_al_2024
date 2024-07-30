

function FigureS8A_B_LeftPanel_HMM_20230504_scatter_choice_error

%repeat
process_HMM_20230504_surprise_at_choice_ver2_depth_error('repeat_OFC_20240427');
process_HMM_20230504_surprise_at_choice_ver2_depth_error('repeat_Hippo_20240427');
process_HMM_20230504_surprise_at_choice_ver2_depth_error('repeat_AC_20240427');
process_HMM_20230504_surprise_at_choice_ver2_depth_error('repeat_PPC_20240427');

%altern
process_HMM_20230504_surprise_at_choice_ver2_depth_error('altern_OFC_20240427');
process_HMM_20230504_surprise_at_choice_ver2_depth_error('altern_Hippo_20240427');
process_HMM_20230504_surprise_at_choice_ver2_depth_error('altern_AC_20240427');
process_HMM_20230504_surprise_at_choice_ver2_depth_error('altern_PPC_20240427');
process_HMM_20230504_surprise_at_choice_ver2_depth_error('altern_M1_20240427');
process_HMM_20230504_surprise_at_choice_ver2_depth_error('altern_STR_20240427');


function process_HMM_20230504_surprise_at_choice_ver2_depth_error(folders, kaiseki_number)

switch nargin
    case 0
        hoge
    case 1
        kaiseki_number = 3;
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end

brainarea = folders(1:10);

[analysis_dir,depth_def] = eval(folders);

all_during_choice = [];
all_p_during_choice = [];
all_norm_spike = [];
all_p_surprise = [];
all_sound_correct_error = [];
all_p_sound_correct_error = [];
all_neuron_choice_category = [];

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [choice_neuron, p_choice_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category] = ...
        HMM_ephys_20230504_surprise_at_choice_depth_error(temp_dir, kaiseki_number,depth_def);

    all_during_choice = [all_during_choice; choice_neuron];
    all_p_during_choice = [all_p_during_choice; p_choice_neuron];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];

    all_sound_correct_error = [all_sound_correct_error; sound_correct_error];
    all_p_sound_correct_error = [all_p_sound_correct_error; p_sound_correct_error];

    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];
end
delete(gcp('nocreate'))

disp([length(all_during_choice),length(all_norm_spike),length(all_p_surprise)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot about low, high neurons

right_neuron = find(all_sound_correct_error(:,1) >= 0);
left_neuron = find(all_sound_correct_error(:,1) < 0);
sig_sound = find(all_p_during_choice(:,3) < 0.01); %Sig diff. between left and right

right_sig_neuron = intersect(right_neuron, sig_sound);
left_sig_neuron = intersect(left_neuron, sig_sound);
non_sig_neuron = setdiff(1:length(all_during_choice), sig_sound);
check_neuron = length(right_sig_neuron) + length(left_sig_neuron) + length(non_sig_neuron);
if check_neuron ~= length(all_during_choice)
    hoge
else
    disp([length(left_sig_neuron),length(right_sig_neuron),length(non_sig_neuron)])
end
%PostLowCorrect_leftCorrect PostHighCorrect_leftCorrect PostLowCorrect_leftError PostHighCorrect_leftError
left_choice_matrix = [33 37 34 38];
left_surprise = [3,4];
%PostLowCorrect_rightCorrect PostHighCorrect_rightCorrect PostLowCorrect_rightError PostHighCorrect_rightError
right_choice_matrix = [35 39 36 40];
right_surprise = [5,6];

[~, ~, Lsig_PostLow_leftError, Lsig_PostHigh_leftError, ~] = ...
    test_surprise_at_choice(all_norm_spike, left_sig_neuron, left_choice_matrix, all_p_surprise(:,left_surprise));
[~, ~, Lsig_PostLow_rightError, Lsig_PostHigh_rightError, ~] = ...
    test_surprise_at_choice(all_norm_spike, left_sig_neuron, right_choice_matrix, all_p_surprise(:,right_surprise));

[~, ~, Rsig_PostLow_leftError, Rsig_PostHigh_leftError, ~] = ...
    test_surprise_at_choice(all_norm_spike, right_sig_neuron, left_choice_matrix, all_p_surprise(:,left_surprise));
[~, ~, Rsig_PostLow_rightError, Rsig_PostHigh_rightError, ~] = ...
    test_surprise_at_choice(all_norm_spike, right_sig_neuron, right_choice_matrix, all_p_surprise(:,right_surprise));

%Activity is already flipped


Repeat_prefer = [Lsig_PostLow_leftError; Rsig_PostHigh_rightError];
Repeat_nonprefer = [Lsig_PostHigh_rightError; Rsig_PostLow_leftError];
Switch_prefer = [Lsig_PostHigh_leftError; Rsig_PostLow_rightError];
Switch_nonprefer = [Lsig_PostLow_rightError; Rsig_PostHigh_leftError];

prefer_sabun = Repeat_prefer-Switch_prefer;
nonprefer_sabun = Repeat_nonprefer-Switch_nonprefer;

p_sabun(1) = signrank(prefer_sabun);
p_sabun(2) = signrank(nonprefer_sabun);

%Error trials
temp1 = length(prefer_sabun);
temp2 = length(nonprefer_sabun);

if temp1 ~= temp2
    hoge
else
    length_neuron = temp1;
end


figure
plot_surprise_choice_prefer(Lsig_PostLow_leftError, Lsig_PostHigh_leftError,[0 0 1]);
hold on
plot_surprise_choice_prefer(Rsig_PostHigh_rightError, Rsig_PostLow_rightError, [1 0 0]);
text(6,2,strcat(string(length_neuron),' neurons'))
set(gca,'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])
title({'Preferred choice';brainarea(8:10)})

figure
plot_surprise_choice_nonprefer(Lsig_PostLow_rightError, Lsig_PostHigh_rightError,[0 0 1]);
hold on
plot_surprise_choice_nonprefer(Rsig_PostHigh_leftError, Rsig_PostLow_leftError,[1 0 0]);
text(6,2,strcat(string(length_neuron),' neurons'))
set(gca,'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])
title({'Non-Preferred choice';brainarea(8:10)})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PostLow_leftCorrect, PostHigh_leftCorrect, PostLow_leftError, PostHigh_leftError, p_surprise...
    ] = test_surprise_at_choice(all_norm_spike, low_sig_neuron, use_number, p_surprise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Surprise at choice:
%PostLowCorrect_leftCorrect PostHighCorrect_leftCorrect PostLowCorrect_leftError PostHighCorrect_leftError
PostLow_leftCorrect = all_norm_spike(low_sig_neuron,use_number(1)); %PostLowCorrect_leftCorrect
PostHigh_leftCorrect = all_norm_spike(low_sig_neuron,use_number(2)); %PostHighCorrect_leftCorrect
PostLow_leftError = all_norm_spike(low_sig_neuron,use_number(3)); %PostLowCorrect_leftError
PostHigh_leftError = all_norm_spike(low_sig_neuron,use_number(4)); %PostHighCorrect_leftError

p_surprise = p_surprise(low_sig_neuron,:);

return

function plot_surprise_choice_prefer(Repeat_prefer, Switch_prefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure
% plot(predict_ave, choice_ave, 'k.')

plot(Repeat_prefer, Switch_prefer, '.','color', plot_color) %same choice, X_repeat or Y_zigzag
hold on
plot([-1, 10],[-1, 10],'k:')
set(gca,'linewidth',1,'xlim',[-1, 10],'ylim',[-1, 10], 'box','off')
set(gcf,'Position',[200 200 215 208])


function plot_surprise_choice_nonprefer(Switch_nonprefer, Repeat_nonprefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure
% plot(predict_ave, choice_ave, 'k.')

plot(Switch_nonprefer, Repeat_nonprefer, '.','color', plot_color)  %same choice, X_zigzag or Y_repeat
hold on
plot([-1, 10],[-1, 10],'k:')
set(gca,'linewidth',1,'xlim',[-1, 10],'ylim',[-1, 10], 'box','off')
set(gcf,'Position',[200 200 215 208])


function [choice_correct_error, p_choice_correct_error, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category] = ...
    HMM_ephys_20230504_surprise_at_choice_depth_error(pathname, kaiseki_number,depth_def)

switch nargin
    case 0
        pathname = pwd;
    case 3
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

temp = dir('HMM_spike_count_neurons_20230316*');
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

%Get sig neurons
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

choice_correct_error(:,1) = use_neuron_index(use_sig_neuron,4); %Correct Low high
choice_correct_error(:,2) = use_neuron_index(use_sig_neuron,5); %Error Low high
choice_correct_error(:,3) = use_neuron_index(use_sig_neuron,1); %choice_index
choice_correct_error(:,4) = use_neuron_index(use_sig_neuron,2); %sound_index
p_choice_correct_error(:,1) = use_p_index(use_sig_neuron,4);
p_choice_correct_error(:,2) = use_p_index(use_sig_neuron,5);
p_choice_correct_error(:,3) = use_p_index(use_sig_neuron,1);
p_choice_correct_error(:,4) = use_p_index(use_sig_neuron,2);

sound_correct_error(:,1) = use_neuron_index(use_sig_neuron,4); %Correct Low high
sound_correct_error(:,2) = use_neuron_index(use_sig_neuron,5); %Error Low high
p_sound_correct_error(:,1) = use_p_index(use_sig_neuron,4);
p_sound_correct_error(:,2) = use_p_index(use_sig_neuron,5);

norm_spike = use_norm_spike(use_sig_neuron,:);
p_surprise = use_p_surprise(use_sig_neuron,:);

% trial_group(1).matrix = left;
% trial_group(2).matrix = right;
temp = norm_spike(:,2) - norm_spike(:,1); %right-left
neuron_choice_category = zeros(length(temp),1);
temp = temp > 0;
neuron_choice_category(temp) = 1; %right=1, left=0

%%%
%Remove NAN
test = [choice_correct_error, norm_spike, sound_correct_error];
test = mean(test,2);
test = isnan(test);
test_nan = find(test == 1);
choice_correct_error(test_nan,:) = [];
p_choice_correct_error(test_nan,:) = [];
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



