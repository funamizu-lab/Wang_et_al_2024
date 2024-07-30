

function Figure7E_S6C_HMM_20230510_trace_at_sound_depth_prefer

% repeat
HMM_20230510_trace_at_sound_depth_prefer('repeat_OFC_20240427',2);
HMM_20230510_trace_at_sound_depth_prefer('repeat_AC_20240427',2);
HMM_20230510_trace_at_sound_depth_prefer('repeat_PPC_20240427',2);

% zigzag
HMM_20230510_trace_at_sound_depth_prefer('altern_OFC_20240427',2);
HMM_20230510_trace_at_sound_depth_prefer('altern_Hippo_20240427',2);
HMM_20230510_trace_at_sound_depth_prefer('altern_AC_20240427',2);
HMM_20230510_trace_at_sound_depth_prefer('altern_PPC_20240427',2);
HMM_20230510_trace_at_sound_depth_prefer('altern_M1_20240427',2);
HMM_20230510_trace_at_sound_depth_prefer('altern_STR_20240427',2);

return




function HMM_20230510_trace_at_sound_depth_prefer(folders, kaiseki_number)

[analysis_dir,depth_def] = eval(folders);

all_norm_spike = [];
all_p_surprise = [];
all_sound_correct_error = [];
all_p_sound_correct_error = [];
all_neuron_choice_category = [];

all_pre_sound = [];
all_p_pre_sound = [];

for j = 1:14
    all_sound_trace(j).matrix = [];
    all_choice_trace(j).matrix = [];
end

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [sound_spike, choice_spike, sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category] = ...
        HMM_ephys_20230510_trace_at_sound_depth(temp_dir, kaiseki_number,depth_def);
    
    for j = 1:14
        all_sound_trace(j).matrix = [all_sound_trace(j).matrix; sound_spike(j).matrix];
        all_choice_trace(j).matrix = [all_choice_trace(j).matrix; choice_spike(j).matrix];
    end
    
    all_pre_sound = [all_pre_sound; sound_neuron];
    all_p_pre_sound = [all_p_pre_sound; p_sound_neuron];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];

    all_sound_correct_error = [all_sound_correct_error; sound_correct_error];
    all_p_sound_correct_error = [all_p_sound_correct_error; p_sound_correct_error];

    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];
end
delete(gcp('nocreate'))


target_analysis = 1;
%target_analysis = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Based on focusing on the sound or choice, should change the analysis way
if target_analysis == 1 %Focus on choice
    %Focus on choice
    target_sig = 1;
    target_choice = [1,2]; %left, right
    target_trace = [7,8,9,10];
elseif target_analysis == 2 %Focus on sound
    %Focus on sound
    target_sig = 2;
    target_choice = [3,4]; %left, right
    target_trace = [11,12,13,14];
end

get_trace_sound_or_choice(target_sig, target_choice, target_trace, ...
    all_neuron_choice_category,all_pre_sound,all_p_pre_sound,all_sound_correct_error,all_sound_trace,all_choice_trace);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_trace_sound_or_choice(target_sig, ~, ~, ...
    ~,all_pre_sound,all_p_pre_sound,all_sound_correct_error,all_sound_trace,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right_neuron = find(all_sound_correct_error(:,1) >= 0);
left_neuron = find(all_sound_correct_error(:,1) < 0);
sig_sound = find(all_p_pre_sound(:,target_sig) < 0.01); %Sig diff. between left and right

right_sig_neuron = intersect(right_neuron, sig_sound);
left_sig_neuron = intersect(left_neuron, sig_sound);
non_sig_neuron = setdiff(1:length(all_pre_sound), sig_sound);
check_neuron = length(right_sig_neuron) + length(left_sig_neuron) + length(non_sig_neuron);
if check_neuron ~= length(all_pre_sound)
    hoge
else
    disp([length(left_sig_neuron),length(right_sig_neuron),length(non_sig_neuron)])
end

length_neuron = length(left_sig_neuron)+length(right_sig_neuron);

make_prefer_trace_2(left_sig_neuron, right_sig_neuron, all_sound_trace,length_neuron);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c_repeat_prefer,c_switch_prefer,c_switch_nonprefer,c_repeat_nonprefer, ...
          s_repeat_prefer,s_switch_prefer,s_switch_nonprefer,s_repeat_nonprefer] ...
           = get_trace_prefer(use_neuron, all_choice_trace, prefer_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_repeat_prefer = all_choice_trace(prefer_order(1,1)).matrix(use_neuron,:);
c_switch_prefer = all_choice_trace(prefer_order(1,2)).matrix(use_neuron,:);
c_switch_nonprefer = all_choice_trace(prefer_order(1,3)).matrix(use_neuron,:);
c_repeat_nonprefer = all_choice_trace(prefer_order(1,4)).matrix(use_neuron,:);

s_repeat_prefer = all_choice_trace(prefer_order(2,1)).matrix(use_neuron,:);
s_switch_prefer = all_choice_trace(prefer_order(2,2)).matrix(use_neuron,:);
s_switch_nonprefer = all_choice_trace(prefer_order(2,3)).matrix(use_neuron,:);
s_repeat_nonprefer = all_choice_trace(prefer_order(2,4)).matrix(use_neuron,:);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_prefer_trace_2(left_sig_neuron, right_sig_neuron, all_sound_trace,length_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

left_order =  [7 8 9 10; 11 12 13 14]; %repeat_prefer, switch_prefer, switch_nonprefer, repeat_nonprefer
right_order = [10 9 8 7; 14 13 12 11]; %repeat_prefer, switch_prefer, switch_nonprefer, repeat_nonprefer


[c_repeat_prefer_l,c_switch_prefer_l,c_switch_nonprefer_l,c_repeat_nonprefer_l, ...
 ~,~,~,~] ...
  = get_trace_prefer(left_sig_neuron, all_sound_trace, left_order);
[c_repeat_prefer_r,c_switch_prefer_r,c_switch_nonprefer_r,c_repeat_nonprefer_r, ...
 ~,~,~,~] ...
  = get_trace_prefer(right_sig_neuron, all_sound_trace, right_order);

c_repeat_prefer = [c_repeat_prefer_l; c_repeat_prefer_r];
c_switch_prefer = [c_switch_prefer_l; c_switch_prefer_r];
c_switch_nonprefer = [c_switch_nonprefer_l; c_switch_nonprefer_r];
c_repeat_nonprefer = [c_repeat_nonprefer_l; c_repeat_nonprefer_r];


figure
plot_mean_se_moto(c_repeat_prefer,[1 0 0],0) %
hold on
plot_mean_se_moto(c_switch_prefer,[255 102 0]./255,0) %
hold on
plot_mean_se_moto(c_switch_nonprefer,[68 168 40]./255,0) %
hold on
plot_mean_se_moto(c_repeat_nonprefer,[0 0 1],0) %
hold on
plot([2000 2000], ylim, ':')
yl = ylim;
text(4000,yl(2)*0.9, strcat(string(length_neuron),' neurons'))
xticklabels({'-2',[],'0',[],'2',[],'4'})
xlabel('Time from sound onset(s)')
ylabel('Zscored activity')
set(gca,'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])

return


function [sound_spike, choice_spike, sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category] = ...
    HMM_ephys_20230510_trace_at_sound_depth(pathname, kaiseki_number,depth_def)
    
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


%The activity depended on the current choice or sound:
%How are they modulated?
%And do they have surprised activity?

%First, focus on the choice or sound index
%Neuron index for low or high sounds
sound_neuron(:,1) = use_neuron_index(use_sig_neuron,1); %Left right: choice
sound_neuron(:,2) = use_neuron_index(use_sig_neuron,2); %Low high: sound
p_sound_neuron(:,1) = use_p_index(use_sig_neuron,1); %Left right: choice
p_sound_neuron(:,2) = use_p_index(use_sig_neuron,2); %Low high

sound_correct_error(:,1) = use_neuron_index(use_sig_neuron,4); %Correct Low high
sound_correct_error(:,2) = use_neuron_index(use_sig_neuron,5); %Error Low high
p_sound_correct_error(:,1) = use_p_index(use_sig_neuron,4);
p_sound_correct_error(:,2) = use_p_index(use_sig_neuron,5);

norm_spike = use_norm_spike(use_sig_neuron,:);
p_surprise = use_p_surprise(use_sig_neuron,:);

temp = norm_spike(:,2) - norm_spike(:,1); %right-left
neuron_choice_category = zeros(length(temp),1);
temp = temp > 0;
neuron_choice_category(temp) = 1; %right=1, left=0

%%%
%Remove NAN
test = [sound_neuron, norm_spike, sound_correct_error];
test = mean(test,2);
test = isnan(test);
test_nan = find(test == 1);
sound_neuron(test_nan,:) = [];
p_sound_neuron(test_nan,:) = [];
norm_spike(test_nan,:) = [];
p_surprise(test_nan,:) = [];
sound_correct_error(test_nan,:) = [];
p_sound_correct_error(test_nan,:) = [];
neuron_choice_category(test_nan,:) = [];

use_sig_neuron(test_nan) = [];

sound_spike = get_spike_trace('HMM_ephys1_20230315_make_ave*', use_sig_neuron);
choice_spike = get_spike_trace('HMM_choice_ephys1_20230315_make_ave*', use_sig_neuron);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sound_spike = get_spike_trace(filename, use_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = dir(filename);
if length(temp) ~= 1
    hoge
end
load(temp.name);
%norm_spike, max_time, mean_spike, std_spike
sound_spike(1).matrix = norm_spike.left(use_sig_neuron,:);
sound_spike(2).matrix = norm_spike.right(use_sig_neuron,:);
sound_spike(3).matrix = norm_spike.left_correct(use_sig_neuron,:);
sound_spike(4).matrix = norm_spike.right_correct(use_sig_neuron,:);
sound_spike(5).matrix = norm_spike.left_error(use_sig_neuron,:);
sound_spike(6).matrix = norm_spike.right_error(use_sig_neuron,:);

sound_spike(7).matrix =  norm_spike.PostLowCorrect_left(use_sig_neuron,:);
sound_spike(8).matrix =  norm_spike.PostHighCorrect_left(use_sig_neuron,:);
sound_spike(9).matrix =  norm_spike.PostLowCorrect_right(use_sig_neuron,:);
sound_spike(10).matrix = norm_spike.PostHighCorrect_right(use_sig_neuron,:);
sound_spike(11).matrix = norm_spike.PostLowCorrect_low(use_sig_neuron,:);
sound_spike(12).matrix = norm_spike.PostHighCorrect_low(use_sig_neuron,:);
sound_spike(13).matrix = norm_spike.PostLowCorrect_high(use_sig_neuron,:);
sound_spike(14).matrix = norm_spike.PostHighCorrect_high(use_sig_neuron,:);


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_neuron = get_sig_neuron_time_window(p_task,new_p_thre,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_p_sound = min(p_task(:,window),[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return

