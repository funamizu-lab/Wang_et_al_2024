
function Figure8F_G_HMM_20230504_all_region_depth_correct

%repeat
[repeat_OFC_prefer, repeat_OFC_nonprefer, repeat_p(1,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('repeat_OFC_20240427');
[repeat_Hippo_prefer, repeat_Hippo_nonprefer, repeat_p(2,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('repeat_Hippo_20240427');
[repeat_AC_prefer, repeat_AC_nonprefer, repeat_p(3,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('repeat_AC_20240427');
[repeat_PPC_prefer, repeat_PPC_nonprefer, repeat_p(4,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('repeat_PPC_20240427');
%zigzag
[zigzag_OFC_prefer, zigzag_OFC_nonprefer, zigzag_p(1,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('altern_OFC_20240427');
[zigzag_Hippo_prefer, zigzag_Hippo_nonprefer, zigzag_p(2,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('altern_Hippo_20240427');
[zigzag_AC_prefer, zigzag_AC_nonprefer, zigzag_p(3,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('altern_AC_20240427');
[zigzag_PPC_prefer, zigzag_PPC_nonprefer, zigzag_p(4,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('altern_PPC_20240427');
[zigzag_M1_prefer, zigzag_M1_nonprefer, zigzag_p(5,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('altern_M1_20240427');
[zigzag_STR_prefer, zigzag_STR_nonprefer, zigzag_p(6,:), ~] = process_HMM_20230504_surprise_at_choice_ver2_depth_correct('altern_STR_20240427');

group_repeat_OFC = ones(length(repeat_OFC_prefer),1) * 1;
group_repeat_HPC = ones(length(repeat_Hippo_prefer),1) * 2;
group_repeat_AC = ones(length(repeat_AC_prefer),1) * 3;
group_repeat_PPC = ones(length(repeat_PPC_prefer),1) * 4;

group_zigzag_OFC = ones(length(zigzag_OFC_prefer),1) * 5;
group_zigzag_Hippo = ones(length(zigzag_Hippo_prefer),1) * 6;
group_zigzag_AC = ones(length(zigzag_AC_prefer),1) * 7;
group_zigzag_PPC = ones(length(zigzag_PPC_prefer),1) * 8;
group_zigzag_M1 = ones(length(zigzag_M1_prefer),1) * 9;
group_zigzag_STR = ones(length(zigzag_STR_prefer),1) * 10;

repeat_prefer = [repeat_OFC_prefer; repeat_Hippo_prefer; repeat_AC_prefer; repeat_PPC_prefer];
repeat_nonprefer = [repeat_OFC_nonprefer; repeat_Hippo_nonprefer; repeat_AC_nonprefer; repeat_PPC_nonprefer];
group_repeat = [group_repeat_OFC; group_repeat_HPC; group_repeat_AC; group_repeat_PPC];
group_repeat2 = group_repeat + 4;
zigzag_prefer = [zigzag_OFC_prefer; zigzag_Hippo_prefer; zigzag_AC_prefer; zigzag_PPC_prefer; zigzag_M1_prefer; zigzag_STR_prefer];
zigzag_nonprefer = [zigzag_OFC_nonprefer; zigzag_Hippo_nonprefer; zigzag_AC_nonprefer; zigzag_PPC_nonprefer; zigzag_M1_nonprefer; zigzag_STR_nonprefer];
group_zigzag = [group_zigzag_OFC; group_zigzag_Hippo; group_zigzag_AC; group_zigzag_PPC; group_zigzag_M1; group_zigzag_STR];
group_zigzag2 = group_zigzag + 6;

figure
boxplot(cat(1, repeat_prefer, repeat_nonprefer), cat(1, group_repeat, group_repeat2),'symbol', '',...
    'Labels',{'OFC','HPC','AC','PPC','OFC','HPC','AC','PPC'},'Colors', [0 0.5 0])
yline(0,':k')
text(1.5, 2.3, 'Preferred Sound')
text(5.5, 2.3, 'Non-preferred Sound')
ylim([-2.5 2.5])
ylabel('Repeat - Switch choice')
title('During choice (Repeating condition)')
set(gca,'FontName', 'Arial')
set(gcf,'Position',[441,351,616,360])



figure
boxplot(cat(1, zigzag_prefer, zigzag_nonprefer), cat(1, group_zigzag, group_zigzag2),'symbol', '',...
      'Labels',{'OFC','HPC','AC','PPC','M1','STR','OFC','HPC','AC','PPC','M1','STR'},'Colors', 'm')
text(2.5, 2.3, 'Preferred Sound')
text(8.5, 2.3, 'Non-preferred Sound')
yline(0,':k')
ylim([-2.5 2.5])
ylabel('Repeat - Switch choice')
title('During choice (Alternating condition)')
set(gca,'FontName', 'Arial')
set(gcf,'Position',[389,433,1044,360])


return


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

[prefer_sabunCorrect, nonprefer_sabunCorrect, p_sabunCorrect] = ...
    plot_surprise_choice(Repeat_preferCorrect, Switch_preferCorrect, Switch_nonpreferCorrect, Repeat_nonpreferCorrect, [0 0 0]);

temp1 = length(prefer_sabunCorrect);
temp2 = length(nonprefer_sabunCorrect);

if temp1 ~= temp2
    hoge
else
    length_neuron = temp1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PostLow_leftCorrect, PostHigh_leftCorrect, PostLow_leftError, PostHigh_leftError, p_surprise...
    ] = test_surprise_at_choice(all_norm_spike, low_sig_neuron, use_number, p_surprise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Surprise at choice:
PostLow_leftCorrect = all_norm_spike(low_sig_neuron,use_number(1)); %PostLowCorrect_leftCorrect
PostHigh_leftCorrect = all_norm_spike(low_sig_neuron,use_number(2)); %PostHighCorrect_leftCorrect
PostLow_leftError = all_norm_spike(low_sig_neuron,use_number(3)); %PostLowCorrect_leftError
PostHigh_leftError = all_norm_spike(low_sig_neuron,use_number(4)); %PostHighCorrect_leftError

p_surprise = p_surprise(low_sig_neuron,:);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefer_sabun, nonprefer_sabun, p_sabun] = plot_surprise_choice(Repeat_prefer, Switch_prefer, Switch_nonprefer, Repeat_nonprefer, plot_color)
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
return













