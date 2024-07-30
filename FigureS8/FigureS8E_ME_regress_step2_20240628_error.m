function FigureS8E_ME_regress_step2_20240628_error

tic


% during choice, error
[Coef_OFC_choice_error,~,p_prefer_choice_error(1,:),p_nonprefer_choice_error(1,:)] = ...
    ME_regress_step1_choice_error_20240628('altern_OFC_20240427');
[Coef_HPC_choice_error,~,p_prefer_choice_error(2,:),p_nonprefer_choice_error(2,:)] = ...
    ME_regress_step1_choice_error_20240628('altern_Hippo_20240427');
[Coef_AC_choice_error,~,p_prefer_choice_error(3,:),p_nonprefer_choice_error(3,:)] = ...
    ME_regress_step1_choice_error_20240628('altern_AC_20240427');
[Coef_PPC_choice_error,~,p_prefer_choice_error(4,:),p_nonprefer_choice_error(4,:)] = ...
    ME_regress_step1_choice_error_20240628('altern_PPC_20240427');
[Coef_M1_choice_error,~,p_prefer_choice_error(5,:),p_nonprefer_choice_error(5,:)]  = ...
    ME_regress_step1_choice_error_20240628('altern_M1_20240427');
[Coef_STR_choice_error,~,p_prefer_choice_error(6,:),p_nonprefer_choice_error(6,:)]  = ...
    ME_regress_step1_choice_error_20240628('altern_STR_20240427');

toc

group_zigzag_OFC3 = ones(length(Coef_OFC_choice_error),1) * 1;
group_zigzag_HPC3 = ones(length(Coef_HPC_choice_error),1) * 2;
group_zigzag_AC3 = ones(length(Coef_AC_choice_error),1) * 3;
group_zigzag_PPC3 = ones(length(Coef_PPC_choice_error),1) * 4;
group_zigzag_M13 = ones(length(Coef_M1_choice_error),1) * 5;
group_zigzag_STR3 = ones(length(Coef_STR_choice_error),1) * 6;


group_choice = [group_zigzag_OFC3; group_zigzag_HPC3;group_zigzag_AC3;group_zigzag_PPC3; group_zigzag_M13; group_zigzag_STR3];
group_choice2 = group_choice+6;


prefer_choice_error = [Coef_OFC_choice_error(:,1); Coef_HPC_choice_error(:,1);Coef_AC_choice_error(:,1);Coef_PPC_choice_error(:,1); Coef_M1_choice_error(:,1); Coef_STR_choice_error(:,1)];
nonprefer_choice_error = [Coef_OFC_choice_error(:,2); Coef_HPC_choice_error(:,2);Coef_AC_choice_error(:,2);Coef_PPC_choice_error(:,2); Coef_M1_choice_error(:,2); Coef_STR_choice_error(:,2)];

figure
boxplot(cat(1, prefer_choice_error, nonprefer_choice_error), cat(1, group_choice, group_choice2),'symbol', '',...
    'Labels',{'OFC','HPC','AC','PPC','M1','STR','OFC','HPC','AC','PPC','M1','STR'})
yline(0,':k')
ylim([-0.8 0.8])
ylabel('Regression coefficient')
title('During choice, Error trials (Alternating condition with Motion Energy recording)')
set(gca,'FontName', 'Arial')

return


function [b_all, p_all, p_prefer_coef,p_nonprefer_coef] = ME_regress_step1_choice_error_20240628(folders)

switch nargin
    case 0
        folders = pwd;
    case 1
        disp('OK to analyze')
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end

b = [];
p = [];


[analysis_dir,depth_def] = eval(folders);

DLC_root = '/Volumes/Extreme SSD/DLC_data';
warning off

for i = 1:length(analysis_dir)
    temp_dir=analysis_dir{i};

    n = strfind(temp_dir,'/20');

    if contains(temp_dir,'AC')
        mouseID = temp_dir(n-6:n-4);
    else
        mouseID = temp_dir(n-7:n-5);
    end
    RecDate = [temp_dir(n+1:n+4),temp_dir(n+6:n+7),temp_dir(n+9:n+10)];
    DLC_test = [DLC_root,'/',mouseID,'/',RecDate];

    if ~isfolder(DLC_test)
        disp('not a dir')
        continue
    else
        disp(temp_dir)
        disp([i,length(analysis_dir)])

        cd(temp_dir)
        temp = dir('ME_neurons_20240530*');
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

        new_p_thre = 10;


        sig_during_choice = get_sig_neuron_time_window(p_task2,new_p_thre,6:15); %During choice 1000ms

        use_sig_neuron = intersect(sig_during_choice,depth_neuron);
        [Lpre,Rpre] = getPreferSide(neuron_index.choice2, p_index.choice2, use_sig_neuron);
    end

    if isempty(Lpre) && isempty(Rpre)
        disp('No significant L R neurons')
        continue
    else
        temp = dir('Bpod_mat*');
        if length(temp) ~= 1
            disp(string(temp))
            hoge
        end
        load(temp.name);
        Bpod_file = temp.name;

        [Choice_trial,~,~,use_trial_remove_first,...
            low,high,correct,error,~,~,...
            ~,~,~,...
            ~,~,~] = ...
            HMM_get_basic_task_structure_20210514(Bpod_file);

        left = find(Chosen_side == 0);
        right = find(Chosen_side == 1);

        temp = dir('task_frame_tokyo_ephys_*');
        if length(temp) ~= 1
            disp(string(temp))
            hoge
        end
        load(temp.name);

        frame_choice_select = nan(length(frame_choice),1);
        frame_choice_select(left) = frame_choice(left,1);
        frame_choice_select(right) = frame_choice(right,2);
        temp = frame_choice_select(Choice_trial);
        if max(isnan(temp)) == 1
            frame_choice_select(Choice_trial)
            length(Choice_trial)
            hoge
            test_nan = isnan(temp);
            %sum(isnan(temp))
            test_nan = find(test_nan == 1);
            Choice_trial(test_nan)
            frame_choice(test_nan,:)
            Chosen_side(test_nan,:)
            hoge
        end

        [MotionEnergy, use_trial] = getMEdata(DLC_test,frame_sound,frame_choice_select,use_trial_remove_first);


        PostLow  = low+1;
        PostHigh = high+1;
        PostLow  = intersect(PostLow, use_trial); %make sure the trials
        PostHigh = intersect(PostHigh, use_trial); %make sure the trials

        Postcorrect = correct+1;
        Postcorrect = intersect(Postcorrect,use_trial);

        PostLowCorrect = intersect(PostLow, Postcorrect);
        PostHighCorrect = intersect(PostHigh, Postcorrect);

        left_error = intersect(left,error);
        right_error = intersect(right,error);

        PostLowCorrect_lefterror  = intersect(PostLowCorrect, left_error);
        PostLowCorrect_righterror = intersect(PostLowCorrect, right_error);
        PostHighCorrect_lefterror  = intersect(PostHighCorrect, left_error);
        PostHighCorrect_righterror = intersect(PostHighCorrect, right_error);

        ChoseData(1).matrix = PostLowCorrect_lefterror;
        ChoseData(2).matrix = PostLowCorrect_righterror;
        ChoseData(3).matrix = PostHighCorrect_lefterror;
        ChoseData(4).matrix = PostHighCorrect_righterror;

        cd(temp_dir)
        spike_dir = dir('spike_ch*');
        if length(spike_dir) ~= 1
            hoge
        end
        spike_dir = spike_dir.name;

        cd(spike_dir)

        tif_name = dir('task_spike_stripe*.mat'); %get all the tif files
        max_tif = length(tif_name);

        time_window_base = 200; %ms
        base_for_spout_on  = [frame_spout(:,1)-time_window_base, frame_spout(:,1)-1]; %Before moving spout

        % time_before_sound = [frame_sound - 600, frame_sound-1];
        %
        % time_all_sound = [frame_sound, frame_sound+600-1];

        time_choice2 = [frame_choice_select, frame_choice_select+1000-1]; %1000 ms

        % [length_trial,~] = size(time_before_sound);

        [y,x,~,~] = ...
            getRegressData(ChoseData,Lpre,Rpre,max_tif,time_choice2,base_for_spout_on,use_trial,MotionEnergy.choice);


        [Coef_1,Pval] = LinMDL(y,x,ChoseData);
        b = [b; Coef_1];
        p = [p; Pval];

    end
end

delta_prefer_coef = b(:,2)-b(:,3);

delta_nonprefer_coef = b(:,5)-b(:,4);

p_prefer_coef = signrank(delta_prefer_coef);

p_nonprefer_coef = signrank(delta_nonprefer_coef);

b_all = [delta_prefer_coef, delta_nonprefer_coef];
p_all = [p_prefer_coef,p_nonprefer_coef];

return


function [MotionEnergy, use_trial] = getMEdata(pathname,frame_sound,frame_choice_select,use_trial_remove_first)

ME = Motion_Energy_cal(pathname);
time = ME.time;

allME = ME.EyeRaw+ME.noseRaw+ME.tongueRaw+ME.whiskerRaw;
allME1 = [0; allME];

allME_0_1 = rescale(allME1);

pre_frame = 2000;
post_frame = 4000;
pre_frame2 = 600;
sound_time = 600;
choice_time = 1000;
trial_time=5000;

time_before_sound = [frame_sound - pre_frame2, frame_sound-1];
time_all_sound = [frame_sound, frame_sound+sound_time-1];
time_choice2 = [frame_choice_select, frame_choice_select+choice_time-1]; %1000 ms

test_choice2 = time_choice2(use_trial_remove_first(end),:);
if time(end) < test_choice2(1)
    use_trial_remove_first(end) = [];
    frame_sound(end) = [];
    disp('reduce_trial')
    disp(pathname)
end


ME_presound  = nan(length(frame_sound),1);
ME_sound  = nan(length(frame_sound),1);
ME_choice = nan(length(frame_sound),1);

for i=1:length(frame_sound)
    time_trial = frame_sound(i)-pre_frame:frame_sound(i)+post_frame;
    time_presound = time_before_sound(i,:);
    time_sound = time_all_sound(i,:);
    time_choice = time_choice2(i,:);

    time_presound = time_presound(1):time_presound(2);
    time_sound = time_sound(1):time_sound(2);
    time_choice = time_choice(1):time_choice(2);

    Frame_trial = ismember(time,time_trial);
    Frame_presound = ismember(time,time_presound);
    Frame_sound = ismember(time,time_sound);
    Frame_choice= ismember(time,time_choice);

    Time_whole_trial = time(Frame_trial);
    ME_presound(i) = mean(allME_0_1(Frame_presound));
    ME_sound(i) = mean(allME_0_1(Frame_sound));
    ME_choice(i)= mean(allME_0_1(Frame_choice));

    if isempty(Time_whole_trial)
        disp([i,length(frame_sound)])
        hoge
    end

    t=Time_whole_trial - Time_whole_trial(1);

    if t(end)< trial_time
        use_trial_remove_first(end) = [];
        continue
    end


end
use_trial = use_trial_remove_first;

MotionEnergy.presound = ME_presound;
MotionEnergy.sound = ME_sound;
MotionEnergy.choice = ME_choice;

return


function [y,x,prefer_diff,nonprefer_diff] = getRegressData(ChoseData,Lprefer,Rprefer,max_tif,time_frame,base_frame,use_trial_remove_first,MotionEnergy)


y_n = [];

% repeat prefer, Switch prefer, switch nonprefer, repeat nonprefer

% ChoseData(1).matrix = PostLowCorrect_left;
% ChoseData(2).matrix = PostLowCorrect_right;
% ChoseData(3).matrix = PostHighCorrect_left;
% ChoseData(4).matrix = PostHighCorrect_right;

TmpPara2 = [1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1];

use_para = TmpPara2;



% base_frame = base_frame(use_trial_remove_first,:);
temp = mean(mean(base_frame));
if isnan(temp)
    disp('nan_detected')
    hoge
end
count = 1;
for file_count = 1 : max_tif
    if ~ismember(file_count,Lprefer) && ~ismember(file_count,Rprefer)
        continue
    else
        temp_file = sprintf('task_spike_stripe20210520_%d',file_count);

        data = load(temp_file);
        spike_mark = data.spike_mark;

        spike_target = nan(length(MotionEnergy),1);
        norm_spike = nan(length(MotionEnergy),1);
        base_spike = nan(length(MotionEnergy),1);
        std_base = nan(length(MotionEnergy),1);

        for i = 1:length(MotionEnergy)
            temp = max(isnan(base_frame(i,:)));
            if temp ~= 1
                base_frame1 = base_frame(i,1) : base_frame(i,2);
                base_spike(i) = mean(spike_mark(base_frame1));
                std_base(i) = std(spike_mark(base_frame1));
            else
                base_spike(i) = nan;
            end
        end

        base_spike1 = base_spike(use_trial_remove_first);
        mean_base = mean(base_spike1);
        std_base1 = std(base_spike1);

        for i = 1:length(MotionEnergy) %trial
            temp = max(isnan(time_frame(i,:)));
            if temp ~= 1
                target_frame1 = time_frame(i,1) : time_frame(i,2);
                spike_target(i) = sum(spike_mark(target_frame1));
                temp_target = spike_target(i);
                norm_spike(i) = (temp_target-mean_base) ./ std_base1;
            else
                norm_spike(i) = nan;
            end

        end

        spike_group = nan(length(ChoseData),1);
        norm_group = nan(length(ChoseData),1);

        for i = 1:length(ChoseData)
            use_trial = ChoseData(i).matrix;
            temp_spike = spike_target(use_trial);
            spike_group(i) = mean(temp_spike);
            norm_group(i) = (mean(temp_spike)-mean_base) ./ std_base1; %Normalized
        end

        yData = spike_target;


        [~,s] = size(use_para);

        tempxPara = nan(length(MotionEnergy),s);

        if ismember(file_count, Lprefer)
            disp('Left preferred neuron')
            tempxPara(ChoseData(1).matrix,:) = repmat(use_para(1,:),length(ChoseData(1).matrix),1);
            tempxPara(ChoseData(2).matrix,:) = repmat(use_para(3,:),length(ChoseData(2).matrix),1);
            tempxPara(ChoseData(3).matrix,:) = repmat(use_para(2,:),length(ChoseData(3).matrix),1);
            tempxPara(ChoseData(4).matrix,:) = repmat(use_para(4,:),length(ChoseData(4).matrix),1);

            TmpRepeat_prefer = norm_group(1);
            TmpSwitch_prefer = norm_group(3);
            TmpSwitch_nonprefer = norm_group(2);
            TmpRepeat_nonprefer = norm_group(4);
        elseif ismember(file_count,Rprefer)
            disp('Right preferred neuron')
            tempxPara(ChoseData(1).matrix,:) = repmat(use_para(4,:),length(ChoseData(1).matrix),1);
            tempxPara(ChoseData(2).matrix,:) = repmat(use_para(2,:),length(ChoseData(2).matrix),1);
            tempxPara(ChoseData(3).matrix,:) = repmat(use_para(3,:),length(ChoseData(3).matrix),1);
            tempxPara(ChoseData(4).matrix,:) = repmat(use_para(1,:),length(ChoseData(4).matrix),1);

            TmpRepeat_prefer = norm_group(4);
            TmpSwitch_prefer = norm_group(2);
            TmpSwitch_nonprefer = norm_group(3);
            TmpRepeat_nonprefer = norm_group(1);
        end

        Tmp_yData = [TmpRepeat_prefer;TmpSwitch_prefer;TmpSwitch_nonprefer;TmpRepeat_nonprefer];

        y_n = [y_n;Tmp_yData'];
    end

    xData = [tempxPara, MotionEnergy];

    xData = xData(use_trial_remove_first,:);
    yData = yData(use_trial_remove_first,:);

    test = mean(xData,2);
    test = isnan(test);
    test_nan = find(test == 1);
    yData(test_nan,:) = [];
    xData(test_nan,:) = [];


    y(count).matrix = yData;
    x(count).matrix = xData;
    count = count + 1;
end

if length(y) ~=length(x)
    disp('Error')
    HereisError
end

prefer_diff = y_n(:,1)-y_n(:,2);
nonprefer_diff = y_n(:,4)-y_n(:,3);

return

function sig_neuron = get_sig_neuron_time_window(p_task,new_p_thre,window)

use_p_sound = min(p_task(:,window),[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return

function [Lpre,Rpre] = getPreferSide(use_neuron_index, use_p_index, use_sig_neuron)

pre_sound(:,1) = use_neuron_index(:,1); %current_choice
pre_sound(:,2) = use_neuron_index(:,6); %pre sound
pre_sound(:,3) = use_neuron_index(:,7); %pre choice
pre_sound(:,4) = use_neuron_index(:,10); %pre SoundCorrect, current left
pre_sound(:,5) = use_neuron_index(:,11); %pre SoundCorrect, current right

pre_sound(:,6) = use_neuron_index(:,2); %current tone index


p_pre_sound(:,1) = use_p_index(:,1); %current choice index
p_pre_sound(:,2) = use_p_index(:,6); %pre tone index
p_pre_sound(:,3) = use_p_index(:,7); %pre choice index
p_pre_sound(:,4) = use_p_index(:,10); %pre SoundCorrect, current left
p_pre_sound(:,5) = use_p_index(:,11); %pre SoundCorrect, current right

p_pre_sound(:,6) = use_p_index(:,2); %current tone index

choice_correct_error(:,1) = use_neuron_index(:,4); %Correct choice index
choice_correct_error(:,2) = use_neuron_index(:,5); %Error choice index
% p_choice_correct_error(:,1) = use_p_index(:,4);
% p_choice_correct_error(:,2) = use_p_index(:,5);



right_neuron = find(choice_correct_error(:,1) >= 0);
left_neuron = find(choice_correct_error(:,1) < 0);
sig_sound = find(p_pre_sound(:,1) < 0.01);

right_sig_neuron = intersect(right_neuron, sig_sound);
left_sig_neuron = intersect(left_neuron, sig_sound);

Rpre = intersect(right_sig_neuron, use_sig_neuron);
Lpre = intersect(left_sig_neuron, use_sig_neuron);

return


function [b,p] = LinMDL(y,x,ChoseData)

length_trial = length(ChoseData(1).matrix)+length(ChoseData(2).matrix)+length(ChoseData(3).matrix)+length(ChoseData(4).matrix);

if length(y) == 1
    disp([length(y.matrix), length_trial, size(x.matrix,1)])
    if length(y.matrix)~=length_trial
        hoge
    elseif size(x.matrix,1)~=length_trial
        hoge
    end
    const = ones(length(y.matrix),1);
    x2 = [const, x.matrix];
    [b,~,stats] = glmfit(x2,y.matrix,'poisson','Constant','off');
    b = b';
    p = stats.p';
else
    [~,s2] = size(x(1).matrix);
    b = nan(length(y),s2+1);
    p = nan(length(y),s2+1);
    for i = 1 : length(y)
        disp([length(y(i).matrix), length_trial, size(x(i).matrix,1)])
        if length(y(i).matrix)~=length_trial
            hoge
        elseif size(x(i).matrix,1)~=length_trial
            hoge
        end
        const = ones(length(y(i).matrix),1);
        x2 = [const, x(i).matrix];
        [b(i,:),~,stats] = glmfit(x2,y(i).matrix,'poisson','Constant','off');
        p(i,:) = stats.p;
    end
end

return

