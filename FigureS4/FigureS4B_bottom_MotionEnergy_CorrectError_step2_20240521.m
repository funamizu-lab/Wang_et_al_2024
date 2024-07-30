function FigureS4B_bottom_MotionEnergy_CorrectError_step2_20240521
ME_correct_session_ave = [];ME_error_session_ave = [];

all_session_num = length(W25_MotionEnergy_path)+length(W27_MotionEnergy_path)+...
    length(W28_MotionEnergy_path)+length(W29_MotionEnergy_path)+...
    length(W30_MotionEnergy_path)+length(W31_MotionEnergy_path);
disp([string(all_session_num) 'sessions'])


[~,allCorrect,allError,~] = MotionEnergy_CorrectError_step1_20240521_process(W25_MotionEnergy_path);
ME_correct_session_ave = [ME_correct_session_ave;mean(allCorrect)];
ME_error_session_ave = [ME_error_session_ave;mean(allError)];

[~,allCorrect,allError,~] = MotionEnergy_CorrectError_step1_20240521_process(W27_MotionEnergy_path);
ME_correct_session_ave = [ME_correct_session_ave;mean(allCorrect)];
ME_error_session_ave = [ME_error_session_ave;mean(allError)];


[~,allCorrect,allError,~] = MotionEnergy_CorrectError_step1_20240521_process(W28_MotionEnergy_path);
ME_correct_session_ave = [ME_correct_session_ave;mean(allCorrect)];
ME_error_session_ave = [ME_error_session_ave;mean(allError)];


[~,allCorrect,allError,~] = MotionEnergy_CorrectError_step1_20240521_process(W29_MotionEnergy_path);
ME_correct_session_ave = [ME_correct_session_ave;mean(allCorrect)];
ME_error_session_ave = [ME_error_session_ave;mean(allError)];


[~,allCorrect,allError,~] = MotionEnergy_CorrectError_step1_20240521_process(W30_MotionEnergy_path);
ME_correct_session_ave = [ME_correct_session_ave;mean(allCorrect)];
ME_error_session_ave = [ME_error_session_ave;mean(allError)];

[~,allCorrect,allError,~] = MotionEnergy_CorrectError_step1_20240521_process(W31_MotionEnergy_path);
ME_correct_session_ave = [ME_correct_session_ave;mean(allCorrect)];
ME_error_session_ave = [ME_error_session_ave;mean(allError)];

all_session_num = length(W25_MotionEnergy_path)+length(W27_MotionEnergy_path)+...
    length(W28_MotionEnergy_path)+length(W29_MotionEnergy_path)+...
    length(W30_MotionEnergy_path)+length(W31_MotionEnergy_path);
disp([string(all_session_num), 'sessions'])

grayColor = [.7 .7 .7];
figure(1)
for i = 1:size(ME_correct_session_ave,1)
    plot(ME_correct_session_ave(i,:),'Color', grayColor)
    hold on
end
hold on
plot_mean_se_moto(ME_correct_session_ave,'r',0)
xticks(0:1000:5000)
xticklabels({'-2000','-1000','0','1000','2000','3000'})
ylim([0.15,0.6])
xlabel('Trial time')
ylabel('Nomalized total motion energy')
title('Correct trials')

figure(2)
for i = 1:size(ME_error_session_ave,1)
    plot(ME_error_session_ave(i,:),'Color', grayColor)
    hold on
end
hold on 
plot_mean_se_moto(ME_error_session_ave,'b',0)
xticks(0:1000:5000)
xticklabels({'-2000','-1000','0','1000','2000','3000'})
ylim([0.15,0.6])
xlabel('Trial time')
ylabel('Nomalized total motion energy')
title('Error trials')


end


function [MotionEnergy,allCorrect,allError,allME] = MotionEnergy_CorrectError_step1_20240521_process(analysis_dir)


for path = 1 : length(analysis_dir)
    disp([path, length(analysis_dir)])
    pathname = (analysis_dir{path});
    cd(pathname);

    pre_frame = 2000;
    post_frame = 4000;
    pre_frame2 = 600;
    sound_time = 600;
    choice_time = 1500;
    trial_time = 5000;

    temp = dir('Bpod_mat*.mat');
    Bpod_file=temp.name;
    load(Bpod_file);

    [~,~,~,use_trial_remove_first,...
        ~,~,correct,error,~,~,~,...
        ~,~,~,~,~] ...
        = HMM_get_basic_task_structure_20210514(Bpod_file);

    left  = find(Chosen_side == 0);
    right = find(Chosen_side == 1);

    temp = dir('task_frame_tokyo_ephys_20220210*');
    load(temp.name)

    frame_choice_select = nan(length(frame_choice),1);
    frame_choice_select(left) = frame_choice(left,1);
    frame_choice_select(right) = frame_choice(right,2);

    time_before_sound = [frame_sound - pre_frame2, frame_sound-1];
    time_all_sound = [frame_sound, frame_sound+sound_time-1];
    time_choice2 = [frame_choice_select, frame_choice_select+choice_time-1]; %1000 ms

    ME = Motion_Energy_cal(pathname);
    time = ME.time;

    allME = ME.EyeRaw+ME.noseRaw+ME.tongueRaw+ME.whiskerRaw;
    allME1 = [0; allME];
    allNormME = ME.EyeNorm+ME.noseNorm+ME.tongueNorm+ME.whiskerNorm;
    allNormME1 = [0; allNormME];

    allNormME = rescale(allNormME1);
   
    use_ME = allNormME;

    test_choice2 = time_choice2(use_trial_remove_first(end),:);
    if time(end) < test_choice2(1)
        reduce_trial = 1;
        disp({'reduce_trial =' string(reduce_trial)})
        use_trial_remove_first(end) = [];
        disp('reduce_trial')
        disp(pathname)
    end


    Correct = ismember(use_trial_remove_first,correct);
    Error = ismember(use_trial_remove_first,error);


    count1=1;
    count2=1;
    ME_presound  = nan(length(use_trial_remove_first),1);
    ME_sound  = nan(length(use_trial_remove_first),1);
    ME_choice = nan(length(use_trial_remove_first),1);

    ME_all = cell(length(use_trial_remove_first),1);
    ME_Correct = cell(length(find(Correct)),1);
    ME_Error = cell(length(find(Error)),1);

    % length(use_trial_remove_first)


    for i=1:length(use_trial_remove_first)
        time_trial = frame_sound(use_trial_remove_first(i),1)-pre_frame:frame_sound(use_trial_remove_first(i),1)+post_frame;
        time_presound = time_before_sound(use_trial_remove_first(i),:);
        time_sound = time_all_sound(use_trial_remove_first(i),:);
        time_choice = time_choice2(use_trial_remove_first(i),:);

        time_presound = time_presound(1):time_presound(2);
        time_sound = time_sound(1):time_sound(2);
        time_choice = time_choice(1):time_choice(2);

        Frame_trial = ismember(time,time_trial);
        Frame_presound = ismember(time,time_presound);
        Frame_sound = ismember(time,time_sound);
        Frame_choice= ismember(time,time_choice);

        Time_whole_trial = time(Frame_trial);
        ME_presound(i) = mean(allME1(Frame_presound));
        ME_sound(i) = mean(allME1(Frame_sound));
        ME_choice(i)= mean(allME1(Frame_choice));

        t=Time_whole_trial - Time_whole_trial(1);

        if t(end)< trial_time
            disp([use_trial_remove_first(i), use_trial_remove_first(end)])
            continue
        else

        %check the 0 frames
        if max(Frame_presound) == 0
            disp([i, use_trial_remove_first(i), length(frame_sound)])
            disp([time(1), time(end)])
            disp([time_presound(1), time_presound(end)])
            hoge
        end
        if max(Frame_sound) == 0
            disp([i, use_trial_remove_first(i), length(frame_sound)])
            disp([time(1), time(end)])
            disp([time_sound(1), time_sound(end)])
            hoge
        end
        if max(Frame_choice) == 0
            disp([i, use_trial_remove_first(i), length(frame_sound)])
            disp([time(1), time(end)])
            disp([time_choice(1), time_choice(end)])
            hoge
        end

        ME_all{i} = interp_ME(trial_time,t,use_ME(Frame_trial));
        if Correct(i)
            ME_Correct{count1} = interp_ME(trial_time,t,use_ME(Frame_trial));
            count1=count1+1;
        end
        if Error(i)
            ME_Error{count2} = interp_ME(trial_time,t,use_ME(Frame_trial));
            count2=count2+1;
        end
        end


    end

    MotionEnergy.presound = ME_presound;
    MotionEnergy.sound = ME_sound;
    MotionEnergy.choice = ME_choice;

    allCorrect= cell2mat(ME_Correct);
    allError= cell2mat(ME_Error);
    allME = cell2mat(ME_all);


end

end


function  chazhi_data= interp_ME(sec,time,ME)
TEMP = griddedInterpolant(time,ME);
Time = 0:1:time(end);
MEvalue = TEMP(Time);
chazhi_data= MEvalue(1:sec);
end
