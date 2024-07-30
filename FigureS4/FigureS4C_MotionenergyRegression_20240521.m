function FigureS4C_MotionenergyRegression_20240521


% repeat : W31
% zigzag : W25, W27, W28, W29, W30, W31
tic

analysis_dir = [W25_MotionEnergy_path;W27_MotionEnergy_path;W28_MotionEnergy_path;...
    W29_MotionEnergy_path;W30_MotionEnergy_path;W31_MotionEnergy_path];


for path = 1 : length(analysis_dir)
    disp([path, length(analysis_dir)])
    [Coef_all(path,:),p_all(path,:)] = get_regress_data(analysis_dir{path});
end

toc

Coef_rep = Coef_all(78:88);
Coef_zig = Coef_all(1:77);

p_rep = p_all(78:88);
p_zig = p_all(1:77);

Coef_presound = nan(length(Coef_all),6);Coef_sound = nan(length(Coef_all),6);Coef_choice = nan(length(Coef_all),6);
Coef_rep_presound = nan(length(Coef_rep),6);Coef_rep_sound= nan(length(Coef_rep),6);Coef_rep_choice= nan(length(Coef_rep),6);
Coef_zig_presound = nan(length(Coef_zig),6);Coef_zig_sound= nan(length(Coef_zig),6);Coef_zig_choice= nan(length(Coef_zig),6);

p_zig_presound = nan(length(Coef_zig),6);p_zig_sound = nan(length(Coef_zig),6);p_zig_choice = nan(length(Coef_zig),6);
p_rep_presound = nan(length(Coef_rep),6);p_rep_sound = nan(length(Coef_rep),6);p_rep_choice = nan(length(Coef_rep),6);

for i = 1 : length(Coef_all)
    Coef_presound(i,:) = Coef_all(i).presound;
    Coef_sound(i,:) = Coef_all(i).sound;
    Coef_choice(i,:) = Coef_all(i).choice;
end

for i = 1 : length(Coef_rep)
    Coef_rep_presound(i,:) = Coef_rep(i).presound;
    Coef_rep_sound(i,:) = Coef_rep(i).sound;
    Coef_rep_choice(i,:) = Coef_rep(i).choice;
    p_rep_presound(i,:) = p_rep(i).presound;
    p_rep_sound(i,:) = p_rep(i).sound;
    p_rep_choice(i,:) = p_rep(i).choice;
end

for i = 1 : length(Coef_zig)
    Coef_zig_presound(i,:) = Coef_zig(i).presound;
    Coef_zig_sound(i,:) = Coef_zig(i).sound;
    Coef_zig_choice(i,:) = Coef_zig(i).choice;
    p_zig_presound(i,:) = p_zig(i).presound;
    p_zig_sound(i,:) = p_zig(i).sound;
    p_zig_choice(i,:) = p_zig(i).choice;
end


for i = 1 : size(Coef_choice,2)
    p_pre_all(:,i) = signrank(Coef_presound(:,i));
    p_sound_all(:,i)= signrank(Coef_sound(:,i));
    p_choice_all(:,i)= signrank(Coef_choice(:,i));
    p_pre_rep(:,i)= signrank(Coef_rep_presound(:,i));
    p_sound_rep(:,i)= signrank(Coef_rep_sound(:,i));
    p_choice_rep(:,i)= signrank(Coef_rep_choice(:,i));
    p_pre_zig(:,i)= signrank(Coef_zig_presound(:,i));
    p_sound_zig(:,i)= signrank(Coef_zig_sound(:,i));
    p_choice_zig(:,i)= signrank(Coef_zig_choice(:,i));
end

figure
plot_regress_result(Coef_rep_presound, Coef_rep_sound, Coef_rep_choice,1)

figure
plot_regress_result(Coef_zig_presound, Coef_zig_sound, Coef_zig_choice,2)




end


function [Coef,p] = get_regress_data(pathname)

cd(pathname)

pre_frame = 2000;
post_frame = 4000;
pre_frame2 = 600;
sound_time = 600;
choice_time = 1000;
trial_time=5000;

temp = dir('Bpod_mat*.mat');
Bpod_file=temp.name;
load(Bpod_file);

[~,~,~,use_trial_remove_first,...
    low,high,correct,error,~,~,~,...
    ~,~,~,~,~] ...
    = HMM_get_basic_task_structure_20210514(Bpod_file);

use_trial = use_trial_remove_first;

% trial_evidence = trial_evidence(use_trial);
% Choice = Chosen_side(use_trial,:);

left  = find(Chosen_side == 0);
right = find(Chosen_side == 1);

PostLow  = low+1;
PostHigh = high+1;
Postleft  = left+1;
Postright = right+1;
Postcorrect  = correct+1;
Posterror = error+1;

PostLow  = intersect(PostLow,use_trial);
PostHigh = intersect(PostHigh,use_trial);
Postleft  = intersect(Postleft,use_trial);
Postright = intersect(Postright,use_trial);
Postcorrect  = intersect(Postcorrect,use_trial);
Posterror = intersect(Posterror,use_trial);

low  = intersect(low,use_trial);
high = intersect(high,use_trial);
left  = intersect(left,use_trial);
right = intersect(right,use_trial);
correct  = intersect(correct,use_trial);
error = intersect(error,use_trial);

PredCurr = ones(length(Chosen_side), 3);
PredCurr(left,1) = -1;
PredCurr(low,2) = -1;
PredCurr(error,3) = -1;
PredCurr = PredCurr(use_trial,:);

PredPre = ones(length(Chosen_side), 3);
PredPre(Postleft,1) = -1;
PredPre(PostLow,2) = -1;
PredPre(Posterror,3) = -1;
PredPre = PredPre(use_trial,:);

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

% allME = ME.EyeRaw+ME.noseRaw+ME.whiskerRaw;
allME = ME.EyeRaw+ME.noseRaw+ME.tongueRaw+ME.whiskerRaw;
allME1 = [0; allME];
allME_0_1 = rescale(allME1);

test_choice2 = time_choice2(use_trial_remove_first(end),:);
if time(end) < test_choice2(1)
    reduce_trial = 1;
    disp(strcat('reduce trial = ', string(reduce_trial)))
    use_trial_remove_first(end) = [];
    PredPre(end,:) = [];
    PredCurr(end,:) = [];
    disp('reduce_trial')
    disp(pathname)
else
    reduce_trial = 0;
    disp(strcat('reduce trial = ', string(reduce_trial)))
end


ME_presound  = nan(length(use_trial_remove_first),1);
ME_sound  = nan(length(use_trial_remove_first),1);
ME_choice = nan(length(use_trial_remove_first),1);



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
    ME_presound(i) = mean(allME_0_1(Frame_presound));
    ME_sound(i) = mean(allME_0_1(Frame_sound));
    ME_choice(i)= mean(allME_0_1(Frame_choice));

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
    end


end


[Coef.presound, p.presound] = LinMdl(ME_presound,PredCurr,PredPre);
[Coef.sound, p.sound] = LinMdl(ME_sound,PredCurr,PredPre);
[Coef.choice,p.choice]= LinMdl(ME_choice,PredCurr,PredPre);
end


function  [Coef,Pval]= LinMdl(ME,Curr,Pre)

X = [Curr,Pre];

y = ME;

data=[X,y];

tbl= table(data(:,1),data(:,2),data(:,3),data(:,4),...
    data(:,5),data(:,6),data(:,7),'VariableNames',...
    {'CurrentChoice','CurrentSound','CurrentOutcome',...
    'PreChoice','PreSound','PreOutcome','MotionEnergy'});

mdl = fitlm(tbl,'MotionEnergy ~ CurrentChoice + CurrentSound + CurrentOutcome + PreChoice + PreSound + PreOutcome ');
Coef = mdl.Coefficients.Estimate(2:end);
Pval = mdl.Coefficients.pValue(2:end);
end


function  chazhi_data= interp_ME(sec,time,ME)
TEMP = griddedInterpolant(time,ME);
Time = 0:1:time(end);
MEvalue = TEMP(Time);
chazhi_data= MEvalue(1:sec);
end

function plot_regress_result(Coef_presound, Coef_sound, Coef_choice,flag)

subplot(1,3,1)
plot_mean_se_moto(Coef_presound,'k',2)
hold on
yline(0,':')
xlim([0,7])
xticks(0:7)
xticklabels({[],'Current choice','Current sound','Current outcome','Pre choice','Pre sound','Pre outcome',[]})
xtickangle(45)
ylabel('Regression Coefficient')
title('Previous Sound')
subplot(1,3,2)
plot_mean_se_moto(Coef_sound,'k',2)
hold on
yline(0,':')
xlim([0,7])
xticks(0:7)
xticklabels({[],'Current choice','Current sound','Current outcome','Pre choice','Pre sound','Pre outcome',[]})
xtickangle(45)
ylabel('Regression Coefficient')
title('During Sound')
subplot(1,3,3)
plot_mean_se_moto(Coef_choice,'k',2)
hold on
yline(0,':')
xlim([0,7])
xticks(0:7)
xticklabels({[],'Current choice','Current sound','Current outcome','Pre choice','Pre sound','Pre outcome',[]})
xtickangle(45)
ylabel('Regression Coefficient')
title('During Choice')

if flag == 1
    sgtitle('Repeating condition')
elseif flag == 2
    sgtitle('Alternating condition')
else
    sgtitle('All sessions')

end

set(gca,'fontname','Arial')

end


