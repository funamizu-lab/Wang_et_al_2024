
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
function Figure8D_S7B_HMM_20230504_scatter_ver2_depth_correct

%repeat
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('repeat_OFC_20230427');
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('repeat_Hippo_20230427');
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('repeat_AC_20230427');
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('repeat_PPC_20230427');

%zigzag
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('zigzag_OFC_20230427');
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('zigzag_Hippo_20230427');
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('zigzag_AC_20230427');
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('zigzag_PPC_20230427');
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('zigzag_M1_20230427');
process_HMM_20230504_surprise_at_choice_ver2_depth_correct('zigzag_STR_20230427');


function [prefer_sabunCorrect, nonprefer_sabunCorrect, p_sabunCorrect, length_neuron] = ...
   process_HMM_20230504_surprise_at_choice_ver2_depth_correct(folders, kaiseki_number)

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

[analysis_dir,depth_def] = eval(folders);

brainarea = folders(1:10);

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
        HMM_ephys_20230504_surprise_at_choice_depth(temp_dir, kaiseki_number,depth_def);

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
% right_neuron = find(all_neuron_choice_category == 1);
% left_neuron = find(all_neuron_choice_category == 0);
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

%%%%%%%%%%%%%%%%%%%%%%%%We first need to know whether they are choice or
%%%%%%%%%%%%%%%%%%%%%%%%sound representing

%PostLowCorrect_leftCorrect PostHighCorrect_leftCorrect PostLowCorrect_leftError PostHighCorrect_leftError
left_choice_matrix = [33 37 34 38];
left_surprise = [3,4];
%PostLowCorrect_rightCorrect PostHighCorrect_rightCorrect PostLowCorrect_rightError PostHighCorrect_rightError
right_choice_matrix = [35 39 36 40];
right_surprise = [5,6];

[Lsig_PostLow_leftCorrect, Lsig_PostHigh_leftCorrect, ~, ~, ~] = ...
    test_surprise_at_choice(all_norm_spike, left_sig_neuron, left_choice_matrix, all_p_surprise(:,left_surprise));
[Lsig_PostLow_rightCorrect, Lsig_PostHigh_rightCorrect, ~, ~, ~] = ...
    test_surprise_at_choice(all_norm_spike, left_sig_neuron, right_choice_matrix, all_p_surprise(:,right_surprise));

[Rsig_PostLow_leftCorrect, Rsig_PostHigh_leftCorrect, ~, ~, ~] = ...
    test_surprise_at_choice(all_norm_spike, right_sig_neuron, left_choice_matrix, all_p_surprise(:,left_surprise));
[Rsig_PostLow_rightCorrect, Rsig_PostHigh_rightCorrect, ~, ~, ~] = ...
    test_surprise_at_choice(all_norm_spike, right_sig_neuron, right_choice_matrix, all_p_surprise(:,right_surprise));

%Activity is already flipped
Repeat_preferCorrect = [Lsig_PostLow_leftCorrect; Rsig_PostHigh_rightCorrect];
Repeat_nonpreferCorrect = [Lsig_PostHigh_rightCorrect; Rsig_PostLow_leftCorrect];
Switch_preferCorrect = [Lsig_PostHigh_leftCorrect; Rsig_PostLow_rightCorrect];
Switch_nonpreferCorrect = [Lsig_PostLow_rightCorrect; Rsig_PostHigh_leftCorrect];

clear P_preferCorrect P_nonpreferCorrect P_preferError P_nonpreferError

%Correct trials

[prefer_sabunCorrect, nonprefer_sabunCorrect, p_sabunCorrect] = ...
    plot_surprise_choice(Repeat_preferCorrect, Switch_preferCorrect, Switch_nonpreferCorrect, Repeat_nonpreferCorrect, [0 0 0]);

temp1 = length(prefer_sabunCorrect);
temp2 = length(nonprefer_sabunCorrect);

if temp1 ~= temp2
    hoge
else
    length_neuron = temp1;
end


figure
plot_surprise_choice_prefer(Lsig_PostLow_leftCorrect, Lsig_PostHigh_leftCorrect,[0 0 1]);
hold on
plot_surprise_choice_prefer(Rsig_PostHigh_rightCorrect, Rsig_PostLow_rightCorrect, [1 0 0]);
text(6,0,['p = ' string(p_sabunCorrect(1))])
text(6,2,strcat(string(length_neuron),' neurons'))
set(gca,'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])
title({'Preferred choice';brainarea(8:10)})


figure
plot_surprise_choice_nonprefer(Lsig_PostHigh_rightCorrect,Lsig_PostLow_rightCorrect, [0 0 1]);
hold on
plot_surprise_choice_nonprefer(Rsig_PostLow_leftCorrect, Rsig_PostHigh_leftCorrect, [1 0 0]);
text(6,0,['p = ' string(p_sabunCorrect(2))])
text(6,2,strcat(string(length_neuron),' neurons'))
set(gca,'fontname','Arial','Box','off')
set(gcf,'Position',[584,652,295,263])
title({'Non-Preferred choice';brainarea(8:10)})

return


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefer_sabun, nonprefer_sabun, p_sabun] = plot_surprise_choice(Repeat_prefer, Switch_prefer, Switch_nonprefer, Repeat_nonprefer, ~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


prefer_sabun = Repeat_prefer-Switch_prefer;
nonprefer_sabun = Repeat_nonprefer - Switch_nonprefer;
if ~isempty(prefer_sabun)
    p_sabun(1) = signrank(prefer_sabun);
    p_sabun(2) = signrank(nonprefer_sabun);
  
    signrank(Repeat_prefer, Switch_prefer)

    signrank(Switch_nonprefer, Repeat_nonprefer)
else
    p_sabun = nan(1,2);
end

disp([length(Repeat_prefer),length(Switch_prefer),length(Switch_nonprefer), length(Repeat_nonprefer)])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


