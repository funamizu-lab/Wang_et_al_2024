
function Figure6B_S5C_20240710_trace_before_sound_depth


% % repeat
process_20240710_trace_before_sound_depth('repeat_OFC_20240427');
process_20240710_trace_before_sound_depth('repeat_AC_20240427');

% % zigzag
process_20240710_trace_before_sound_depth('altern_OFC_20240427');
process_20240710_trace_before_sound_depth('altern_AC_20240427');
process_20240710_trace_before_sound_depth('altern_PPC_20240427');
process_20240710_trace_before_sound_depth('altern_M1_20240427');
process_20240710_trace_before_sound_depth('altern_STR_20240427');

return

function process_20240710_trace_before_sound_depth(folders, kaiseki_number)

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
brainarea = folders(1:10);

[analysis_dir,depth_def] = eval(folders);


all_pre_sound = [];
all_p_pre_sound = [];
all_norm_spike = [];
all_p_surprise = [];
all_choice_correct_error = [];
all_p_choice_correct_error = [];
all_neuron_choice_category = [];

all_tone_correct_error = [];
all_p_tone_correct_error = [];

all_pre_tone_correct_error = [];
all_p_pre_tone_correct_error = [];

all_pre_choice_correct_error = [];
all_p_pre_choice_correct_error = [];

for j = 1:14
    all_sound_trace(j).matrix = [];
    all_choice_trace(j).matrix = [];
end

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [sound_spike, choice_spike, pre_sound, p_pre_sound, norm_spike, p_surprise, ...
    choice_correct_error, p_choice_correct_error, neuron_choice_category,...
    tone_correct_error, p_tone_correct_error, pre_tone_correct_error, p_pre_tone_correct_error,...
    pre_choice_correct_error, p_pre_choice_correct_error,...
     ~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
        HMM_ephys_20240710_trace_before_sound_depth(temp_dir, kaiseki_number,depth_def);

    for j = 1:14
        all_sound_trace(j).matrix = [all_sound_trace(j).matrix; sound_spike(j).matrix];
        all_choice_trace(j).matrix = [all_choice_trace(j).matrix; choice_spike(j).matrix];
    end

%     p_pre_sound(:,1) = use_p_index(use_sig_neuron,1); %current choice index
%     p_pre_sound(:,2) = use_p_index(use_sig_neuron,6); %pre tone index
%     p_pre_sound(:,3) = use_p_index(use_sig_neuron,7); %pre choice index
%     p_pre_sound(:,4) = use_p_index(use_sig_neuron,10); %pre SoundCorrect, current left
%     p_pre_sound(:,5) = use_p_index(use_sig_neuron,11); %pre SoundCorrect, current right
%     p_pre_sound(:,6) = use_p_index(use_sig_neuron,11); % current tone index

    all_pre_sound = [all_pre_sound; pre_sound];
    all_p_pre_sound = [all_p_pre_sound; p_pre_sound];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];

    all_choice_correct_error = [all_choice_correct_error; choice_correct_error];
    all_p_choice_correct_error = [all_p_choice_correct_error; p_choice_correct_error];

    all_tone_correct_error = [all_tone_correct_error; tone_correct_error]; % current tone correct index
    all_p_tone_correct_error = [all_p_tone_correct_error; p_tone_correct_error]; % current tone error index
       
    all_pre_tone_correct_error = [all_pre_tone_correct_error; pre_tone_correct_error]; % pre tone correct index
    all_p_pre_tone_correct_error = [all_p_pre_tone_correct_error; p_pre_tone_correct_error]; % pre tone error index
        
    all_pre_choice_correct_error = [all_pre_choice_correct_error; pre_choice_correct_error]; % pre choice correct index
    all_p_pre_choice_correct_error = [all_p_pre_choice_correct_error; p_pre_choice_correct_error]; % pre choice error index

    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];

end
delete(gcp('nocreate'))

disp([length(all_pre_sound),length(all_norm_spike),length(all_p_surprise)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot about low, high neurons


if contains (brainarea,'repeat')
    right_neuron = find(all_pre_choice_correct_error(:,1) >= 0);
    left_neuron = find(all_pre_choice_correct_error(:,1) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between left and right
else
    left_neuron = find(all_pre_choice_correct_error(:,1) >= 0);
    right_neuron = find(all_pre_choice_correct_error(:,1) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between left and right
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

length_neuron = length(left_sig_neuron)+length(right_sig_neuron);

make_prefer_trace(left_sig_neuron, right_sig_neuron, all_sound_trace,brainarea,length_neuron);

return


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_neuron_trace_at_choice(use_neuron, all_sound_trace, all_choice_trace)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1)
plot(mean(all_sound_trace(1).matrix(use_neuron,:)),'b') %left
hold on
plot(mean(all_sound_trace(2).matrix(use_neuron,:)),'r') %right
subplot(2,2,2)
plot(mean(all_choice_trace(1).matrix(use_neuron,:)),'b') %left
hold on
plot(mean(all_choice_trace(2).matrix(use_neuron,:)),'r') %right

subplot(2,2,3)
plot(mean(all_sound_trace(3).matrix(use_neuron,:)),'b') %left_correct
hold on
plot(mean(all_sound_trace(4).matrix(use_neuron,:)),'r') %right_correct
hold on
plot(mean(all_sound_trace(5).matrix(use_neuron,:)),'c') %left_error
hold on
plot(mean(all_sound_trace(6).matrix(use_neuron,:)),'m') %right_error
subplot(2,2,4)
plot(mean(all_choice_trace(3).matrix(use_neuron,:)),'b') %left_correct
hold on
plot(mean(all_choice_trace(4).matrix(use_neuron,:)),'r') %right_correct
hold on
plot(mean(all_choice_trace(5).matrix(use_neuron,:)),'c') %left_error
hold on
plot(mean(all_choice_trace(6).matrix(use_neuron,:)),'m') %right_error

return


function [sound_spike, choice_spike, pre_sound, p_pre_sound, norm_spike, p_surprise, ...
    choice_correct_error, p_choice_correct_error, neuron_choice_category, ...
    tone_correct_error, p_tone_correct_error, pre_tone_correct_error, p_pre_tone_correct_error,...
    pre_choice_correct_error, p_pre_choice_correct_error,...
    SameDif_preTone_correct_error,p_SameDif_preTone_correct_error,...
    preCorrectTone_index_correct_correct_error,p_preCorrectTone_index_correct_correct_error,...
    preCorrectChoice_index_correct_correct_error,p_preCorrectChoice_index_correct_correct_error,...
    pre_same_pLcL_pLcH_index,p_pre_same_pLcL_pLcH_index,...
    pre_same_pHcL_pHcH_index,p_pre_same_pHcL_pHcH_index,...
    cur_same_pLcL_pHcL_index,p_cur_same_pLcL_pHcL_index,...
    cur_same_pLcH_pHcH_index,p_cur_same_pLcH_pHcH_index] = ...
    HMM_ephys_20240710_trace_before_sound_depth(pathname, kaiseki_number,depth_def)

switch nargin
    case 0
        pathname = pwd;
    case 3
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

temp = dir('HMM_spike_count_neurons_20240120*');
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
%neuron_index.before_sound
%neuron_index.init_sound
%neuron_index.late_sound
%neuron_index.all_sound
%neuron_index.before_choice
%neuron_index.choice1
%neuron_index.choice2
%neuron_index.choice3
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

pre_sound(:,6) = use_neuron_index(use_sig_neuron,2); %current tone index


p_pre_sound(:,1) = use_p_index(use_sig_neuron,1); %current choice index
p_pre_sound(:,2) = use_p_index(use_sig_neuron,6); %pre tone index
p_pre_sound(:,3) = use_p_index(use_sig_neuron,7); %pre choice index
p_pre_sound(:,4) = use_p_index(use_sig_neuron,10); %pre SoundCorrect, current left
p_pre_sound(:,5) = use_p_index(use_sig_neuron,11); %pre SoundCorrect, current right

p_pre_sound(:,6) = use_p_index(use_sig_neuron,2); %current tone index

choice_correct_error(:,1) = use_neuron_index(use_sig_neuron,4); %Correct choice index
choice_correct_error(:,2) = use_neuron_index(use_sig_neuron,5); %Error choice index
p_choice_correct_error(:,1) = use_p_index(use_sig_neuron,4);
p_choice_correct_error(:,2) = use_p_index(use_sig_neuron,5);

tone_correct_error(:,1) = use_neuron_index(use_sig_neuron,12); %Correct choice index
tone_correct_error(:,2) = use_neuron_index(use_sig_neuron,13); %Error choice index
p_tone_correct_error(:,1) = use_p_index(use_sig_neuron,12);
p_tone_correct_error(:,2) = use_p_index(use_sig_neuron,13);

pre_tone_correct_error(:,1) = use_neuron_index(use_sig_neuron,14); %Correct choice index
pre_tone_correct_error(:,2) = use_neuron_index(use_sig_neuron,15); %Error choice index
p_pre_tone_correct_error(:,1) = use_p_index(use_sig_neuron,14);
p_pre_tone_correct_error(:,2) = use_p_index(use_sig_neuron,15);

pre_choice_correct_error(:,1) = use_neuron_index(use_sig_neuron,16); %Correct choice index
pre_choice_correct_error(:,2) = use_neuron_index(use_sig_neuron,17); %Error choice index
p_pre_choice_correct_error(:,1) = use_p_index(use_sig_neuron,16);
p_pre_choice_correct_error(:,2) = use_p_index(use_sig_neuron,17);

SameDif_preTone_correct_error(:,1) = use_neuron_index(use_sig_neuron,18);
SameDif_preTone_correct_error(:,2) = use_neuron_index(use_sig_neuron,19);
p_SameDif_preTone_correct_error(:,1) = use_p_index(use_sig_neuron,18);
p_SameDif_preTone_correct_error(:,2) = use_p_index(use_sig_neuron,19);


preCorrectTone_index_correct_correct_error(:,1) = use_neuron_index(use_sig_neuron,20);
preCorrectTone_index_correct_correct_error(:,2) = use_neuron_index(use_sig_neuron,21);
p_preCorrectTone_index_correct_correct_error(:,1) = use_p_index(use_sig_neuron,20);
p_preCorrectTone_index_correct_correct_error(:,2) = use_p_index(use_sig_neuron,21);

preCorrectChoice_index_correct_correct_error(:,1) = use_neuron_index(use_sig_neuron,22);
preCorrectChoice_index_correct_correct_error(:,2) = use_neuron_index(use_sig_neuron,23);
p_preCorrectChoice_index_correct_correct_error(:,1) = use_p_index(use_sig_neuron,22);
p_preCorrectChoice_index_correct_correct_error(:,2) = use_p_index(use_sig_neuron,23);


pre_same_pLcL_pLcH_index(:,1) = use_neuron_index(use_sig_neuron,24);
pre_same_pHcL_pHcH_index(:,1) = use_neuron_index(use_sig_neuron,25);
cur_same_pLcL_pHcL_index(:,1) = use_neuron_index(use_sig_neuron,26);
cur_same_pLcH_pHcH_index(:,1) = use_neuron_index(use_sig_neuron,27);

p_pre_same_pLcL_pLcH_index(:,1) = use_p_index(use_sig_neuron,24);
p_pre_same_pHcL_pHcH_index(:,1) = use_p_index(use_sig_neuron,25);
p_cur_same_pLcL_pHcL_index(:,1) = use_p_index(use_sig_neuron,26);
p_cur_same_pLcH_pHcH_index(:,1) = use_p_index(use_sig_neuron,27);

norm_spike = use_norm_spike(use_sig_neuron,:);
p_surprise = use_p_surprise(use_sig_neuron,:);

temp = norm_spike(:,2) - norm_spike(:,1); %right-left
neuron_choice_category = zeros(length(temp),1);
temp = find(temp > 0);
neuron_choice_category(temp) = 1; %right=1, left=0

%%%
%Remove NAN
test = [pre_sound, norm_spike, choice_correct_error,tone_correct_error,...
    pre_tone_correct_error,pre_choice_correct_error,SameDif_preTone_correct_error, p_SameDif_preTone_correct_error,...
    preCorrectTone_index_correct_correct_error,p_preCorrectTone_index_correct_correct_error,...
    preCorrectChoice_index_correct_correct_error,p_preCorrectChoice_index_correct_correct_error,...
    pre_same_pLcL_pLcH_index, p_pre_same_pLcL_pLcH_index,...
    pre_same_pHcL_pHcH_index,p_pre_same_pHcL_pHcH_index,...
    cur_same_pLcL_pHcL_index, p_cur_same_pLcL_pHcL_index,...
    cur_same_pLcH_pHcH_index,p_cur_same_pLcH_pHcH_index ];
test = mean(test,2);
test = isnan(test);
test_nan = find(test == 1);
pre_sound(test_nan,:) = [];
p_pre_sound(test_nan,:) = [];
norm_spike(test_nan,:) = [];
p_surprise(test_nan,:) = [];
choice_correct_error(test_nan,:) = [];
p_choice_correct_error(test_nan,:) = [];
neuron_choice_category(test_nan,:) = [];

tone_correct_error(test_nan,:) = []; %Correct choice index
p_tone_correct_error(test_nan,:) = [];

pre_tone_correct_error(test_nan,:) = []; %Correct choice index
p_pre_tone_correct_error(test_nan,:) = [];

pre_choice_correct_error(test_nan,:) = [];
p_pre_choice_correct_error(test_nan,:) = [];

SameDif_preTone_correct_error(test_nan,:) = [];
p_SameDif_preTone_correct_error(test_nan,:) = [];

preCorrectTone_index_correct_correct_error(test_nan,:) = [];
p_preCorrectTone_index_correct_correct_error(test_nan,:) = [];

preCorrectChoice_index_correct_correct_error(test_nan,:) = [];
p_preCorrectChoice_index_correct_correct_error(test_nan,:) = [];

pre_same_pLcL_pLcH_index(test_nan,:) = [];
p_pre_same_pLcL_pLcH_index(test_nan,:) = [];
pre_same_pHcL_pHcH_index(test_nan,:) = [];
p_pre_same_pHcL_pHcH_index(test_nan,:) = [];
cur_same_pLcL_pHcL_index(test_nan,:) = [];
p_cur_same_pLcL_pHcL_index(test_nan,:) = [];
cur_same_pLcH_pHcH_index(test_nan,:) = [];
p_cur_same_pLcH_pHcH_index(test_nan,:) = [];

use_sig_neuron(test_nan) = [];

sound_spike = get_spike_trace('HMM_ephys1_20230315_make_ave*', use_sig_neuron);
choice_spike = get_spike_trace('HMM_choice_ephys1_20230315_make_ave*', use_sig_neuron);


return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sound_spike = get_spike_trace(filename, use_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Activity aligned to sound onset
%temp = dir('HMM_ephys1_20230315_make_ave*');
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

% trial_group(25).matrix = PostLowCorrect_left;
% trial_group(26).matrix = PostHighCorrect_left;
% trial_group(27).matrix = PostLowCorrect_right;
% trial_group(28).matrix = PostHighCorrect_right;
% trial_group(29).matrix = PostLowCorrect_low;
% trial_group(30).matrix = PostHighCorrect_low;
% trial_group(31).matrix = PostLowCorrect_high;
% trial_group(32).matrix = PostHighCorrect_high;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_neuron = get_sig_neuron_time_window(p_task,new_p_thre,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_p_sound = min(p_task(:,window),[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_prefer_trace(left_sig_neuron, right_sig_neuron, all_sound_trace,brainarea,length_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if contains (brainarea,'repeat')
    left_order =  [7 9 8 10; 11 12 13 14]; 
    right_order = [10 8 9 7; 14 13 12 11]; 
else
    left_order = [10 8 9 7; 14 13 12 11]; 
    right_order =  [7 9 8 10; 11 12 13 14]; 
end

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


