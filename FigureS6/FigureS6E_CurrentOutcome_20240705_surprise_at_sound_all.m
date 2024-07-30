
function FigureS6E_CurrentOutcome_20240705_surprise_at_sound_all

% %  repeat
[~,~,repeat_OFC_prefer_correct, repeat_OFC_nonprefer_correct, repeat_OFC_prefer_error, repeat_OFC_nonprefer_error, ...
    repeat_p(1,:), repeat_neuron(1)] = process_CurrentOutcome_20240705_surprise_during_sound_depth('repeat_OFC_20240427');

[~,~,repeat_AC_prefer_correct, repeat_AC_nonprefer_correct, repeat_AC_prefer_error, repeat_AC_nonprefer_error, ...
    repeat_p(2,:), repeat_neuron(2)] = process_CurrentOutcome_20240705_surprise_during_sound_depth('repeat_AC_20240427');
% 
[~,~,repeat_PPC_prefer_correct, repeat_PPC_nonprefer_correct, repeat_PPC_prefer_error, repeat_PPC_nonprefer_error, ...
    repeat_p(3,:), repeat_neuron(3)] = process_CurrentOutcome_20240705_surprise_during_sound_depth('repeat_PPC_20240427');


% %  zigzag

[~,~,zigzag_OFC_prefer_correct, zigzag_OFC_nonprefer_correct, zigzag_OFC_prefer_error, zigzag_OFC_nonprefer_error,...
    zigzag_p(1,:), zigzag_neuron(1)] = process_CurrentOutcome_20240705_surprise_during_sound_depth('altern_OFC_20240427');

[~,~,zigzag_Hippo_prefer_correct, zigzag_Hippo_nonprefer_correct, zigzag_Hippo_prefer_error, zigzag_Hippo_nonprefer_error,...
    zigzag_p(2,:), zigzag_neuron(2)] = process_CurrentOutcome_20240705_surprise_during_sound_depth('altern_Hippo_20240427');

[~,~,zigzag_AC_prefer_correct, zigzag_AC_nonprefer_correct, zigzag_AC_prefer_error, zigzag_AC_nonprefer_error,...
    zigzag_p(3,:), zigzag_neuron(3)] = process_CurrentOutcome_20240705_surprise_during_sound_depth('altern_AC_20240427');

[~,~,zigzag_PPC_prefer_correct, zigzag_PPC_nonprefer_correct, zigzag_PPC_prefer_error, zigzag_PPC_nonprefer_error,...
    zigzag_p(4,:), zigzag_neuron(4)] = process_CurrentOutcome_20240705_surprise_during_sound_depth('altern_PPC_20240427');

[~,~,zigzag_M1_prefer_correct, zigzag_M1_nonprefer_correct, zigzag_M1_prefer_error, zigzag_M1_nonprefer_error,...
    zigzag_p(5,:), zigzag_neuron(5)] = process_CurrentOutcome_20240705_surprise_during_sound_depth('altern_M1_20240427');

[~,~,zigzag_STR_prefer_correct, zigzag_STR_nonprefer_correct, zigzag_STR_prefer_error, zigzag_STR_nonprefer_error,...
    zigzag_p(6,:), zigzag_neuron(6)] = process_CurrentOutcome_20240705_surprise_during_sound_depth('altern_STR_20240427');

close all

repeat_prefer_OFC =[repeat_OFC_prefer_correct,repeat_OFC_prefer_error];
repeat_nonprefer_OFC =[repeat_OFC_nonprefer_correct,repeat_OFC_nonprefer_error];
repeat_prefer_AC =[repeat_AC_prefer_correct,repeat_AC_prefer_error];
repeat_nonprefer_AC =[repeat_AC_nonprefer_correct,repeat_AC_nonprefer_error];
repeat_prefer_PPC =[repeat_PPC_prefer_correct,repeat_PPC_prefer_error];
repeat_nonprefer_PPC =[repeat_PPC_nonprefer_correct,repeat_PPC_nonprefer_error];

zigzag_prefer_OFC =[zigzag_OFC_prefer_correct,zigzag_OFC_prefer_error];
zigzag_nonprefer_OFC =[zigzag_OFC_nonprefer_correct,zigzag_OFC_nonprefer_error];
zigzag_prefer_HPC =[zigzag_Hippo_prefer_correct,zigzag_Hippo_prefer_error];
zigzag_nonprefer_HPC =[zigzag_Hippo_nonprefer_correct,zigzag_Hippo_nonprefer_error];
zigzag_prefer_AC =[zigzag_AC_prefer_correct,zigzag_AC_prefer_error];
zigzag_nonprefer_AC =[zigzag_AC_nonprefer_correct,zigzag_AC_nonprefer_error];
zigzag_prefer_PPC =[zigzag_PPC_prefer_correct,zigzag_PPC_prefer_error];
zigzag_nonprefer_PPC =[zigzag_PPC_nonprefer_correct,zigzag_PPC_nonprefer_error];
zigzag_prefer_M1 =[zigzag_M1_prefer_correct,zigzag_M1_prefer_error];
zigzag_nonprefer_M1 =[zigzag_M1_nonprefer_correct,zigzag_M1_nonprefer_error];
zigzag_prefer_STR =[zigzag_STR_prefer_correct,zigzag_STR_prefer_error];
zigzag_nonprefer_STR =[zigzag_STR_nonprefer_correct,zigzag_STR_nonprefer_error];


AArep = cell(1,2);BBrep = cell(1,2);CCrep = cell(1,2);DDrep = cell(1,2);
EErep = cell(1,2);FFrep = cell(1,2);
temp1 = cell(2,6);
temp2 = cell(2,12);
for ii=1:size(temp1,1)
    AArep{ii}=repeat_prefer_OFC(:,ii);
    BBrep{ii}=repeat_prefer_AC(:,ii);
    CCrep{ii}=repeat_prefer_PPC(:,ii);
    DDrep{ii}=repeat_nonprefer_OFC(:,ii);
    EErep{ii}=repeat_nonprefer_AC(:,ii);
    FFrep{ii}=repeat_nonprefer_PPC(:,ii);
end
AAzig = cell(1,2);BBzig = cell(1,2);CCzig = cell(1,2);DDzig = cell(1,2);EEzig = cell(1,2);
FFzig = cell(1,2);GGzig = cell(1,2);HHzig = cell(1,2);IIzig = cell(1,2);JJzig = cell(1,2);
KKzig = cell(1,2);LLzig = cell(1,2);
for ii=1:size(temp2,1)
    AAzig{ii}=zigzag_prefer_OFC(:,ii);
    BBzig{ii}=zigzag_prefer_HPC(:,ii);
    CCzig{ii}=zigzag_prefer_AC(:,ii);
    DDzig{ii}=zigzag_prefer_PPC(:,ii);
    EEzig{ii}=zigzag_prefer_M1(:,ii);
    FFzig{ii}=zigzag_prefer_STR(:,ii);
    GGzig{ii}=zigzag_nonprefer_OFC(:,ii);
    HHzig{ii}=zigzag_nonprefer_HPC(:,ii);
    IIzig{ii}=zigzag_nonprefer_AC(:,ii);
    JJzig{ii}=zigzag_nonprefer_PPC(:,ii);
    KKzig{ii}=zigzag_nonprefer_M1(:,ii);
    LLzig{ii}=zigzag_nonprefer_STR(:,ii);
end

data_repeat=vertcat(AArep,BBrep,CCrep,DDrep,EErep,FFrep);
data_zigzag = vertcat(AAzig,BBzig,CCzig,DDzig,EEzig,FFzig,GGzig,HHzig,IIzig,JJzig,KKzig,LLzig);

Mlab = {'Current correct', 'Current error'};
colorcorrect = [0.4660 0.6740 0.1880]';
colorerror = [0.6350 0.0780 0.1840]';
color = [colorerror,colorcorrect];


figure
multiple_boxplot(data_repeat,{'OFC','AC','PPC','OFC','AC','PPC'},Mlab,color)
hold on
yline(0,':k')
ylim([-3.5 3.5])
legend(fliplr(Mlab),'Location','southeast');
title('Repeating condition (During sound)')
text(1.7, 3.2, 'Preferred Choice')
text(3.7, 3.2, 'Non-preferred Choice')
set(gcf,'Position',[441,351,616,360])
box off


figure
multiple_boxplot(data_zigzag,{'OFC','HPC','AC','PPC','M1','STR','OFC','HPC','AC','PPC','M1','STR'},Mlab,color)
hold on
yline(0,':k')
ylim([-6.5 6.5])
legend(fliplr(Mlab),'Location','southeast');
title('Alternating condition (During sound)')
text(2.5, 5, 'Preferred Choice')
text(7.5, 5, 'Non-preferred Choice')
set(gcf,'Position',[389,433,1044,360])
box off

return


function multiple_boxplot(data,xlab,~,colors)

if ~iscell(data)
    error('Input data is not even a cell array!');
end

% Get sizes
M=size(data,2);
L=size(data,1);
if nargin>=4
    if size(colors,2)~=M
        error('Wrong amount of colors!');
    end
end
if nargin>=2
    if length(xlab)~=L
        error('Wrong amount of X labels given');
    end
end

% Calculate the positions of the boxes
positions=1:0.25:M*L*0.25+1+0.25*L;
positions(1:M+1:end)=[];

% Extract data and label it in the group correctly
x=[];
group=[];
for ii=1:L
    for jj=1:M
        aux=data{ii,jj};
        x=vertcat(x,aux(:));
        group=vertcat(group,ones(size(aux(:)))*jj+(ii-1)*M);
    end
end
% Plot it

boxplot(x,group, 'positions', positions,'Symbol',' ','BoxStyle','outline');

% Set the Xlabels
aux=reshape(positions,M,[]);
labelpos = sum(aux,1)./M;

set(gca,'xtick',labelpos)
if nargin>=2
    set(gca,'xticklabel',xlab);
else
    idx=1:L;
    set(gca,'xticklabel',strsplit(num2str(idx),' '));
end
    

% Get some colors
if nargin>=4
    cmap=colors;
else
    cmap = hsv(M);
    cmap=vertcat(cmap',ones(1,M)*0.5);
end
color=repmat(cmap, 1, L);

% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),color(1:3,jj)','EdgeColor',color(1:3,jj),'FaceColor','none');
end

% if nargin>=3
%     legend(fliplr(Mlab));
% end
return



function [prefer_sabun,nonprefer_sabun,prefer_sabun_correct, nonprefer_sabun_correct, prefer_sabun_error, nonprefer_sabun_error,...
    p_sabun, length_neuron] = process_CurrentOutcome_20240705_surprise_during_sound_depth(folders, kaiseki_number)

switch nargin
    case 0
        hoge
    case 1
        kaiseki_number = 2;
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end

[analysis_dir,depth_def] = eval(folders);


Brainarea = [folders(1:6),' ', folders(8:10)];

all_norm_spike = [];
all_p_surprise = [];
all_sound_correct_error = [];
all_tone_correct_error = [];
all_p_sound_correct_error = [];
all_p_tone_correct_error = [];
all_neuron_choice_category = [];

all_pre_sound = [];
all_p_pre_sound = [];
for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category,tone_correct_error,p_tone_correct_error] = ...
        HMM_sound_trans20240705_surprise_during_depth(temp_dir, kaiseki_number,depth_def);
    
    all_pre_sound = [all_pre_sound; sound_neuron];
    all_p_pre_sound = [all_p_pre_sound; p_sound_neuron];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];

    all_sound_correct_error = [all_sound_correct_error; sound_correct_error];
    all_p_sound_correct_error = [all_p_sound_correct_error; p_sound_correct_error];

    all_tone_correct_error = [all_tone_correct_error;tone_correct_error];
    all_p_tone_correct_error = [all_p_tone_correct_error;p_tone_correct_error];

    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];
end
delete(gcp('nocreate'))

disp([length(all_pre_sound),length(all_norm_spike),length(all_p_surprise)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Based on focusing on the sound or choice, should change the analysis way

%Focus on choice
target_sig = 1;
target_choice = [1,2]; %left, right
target_surprise_both = [25,26,27,28];

target_surprise_correct = [33,37,35,39];
target_surprise_error = [34,38,36,40];


[prefer_sabun,nonprefer_sabun,prefer_sabun_correct, nonprefer_sabun_correct, prefer_sabun_error, nonprefer_sabun_error,p_sabun, length_neuron] = ...
    get_analysis_sound_or_choice(target_sig, target_choice, target_surprise_both,target_surprise_correct, target_surprise_error,...
    all_neuron_choice_category,all_pre_sound,all_p_pre_sound,all_sound_correct_error,all_tone_correct_error,all_norm_spike,Brainarea);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefer_sabun,nonprefer_sabun,prefer_sabun_correct, nonprefer_sabun_correct, prefer_sabun_error, nonprefer_sabun_error,p_sabun, length_neuron] = ...
    get_analysis_sound_or_choice(target_sig, ~, target_surprise_both,target_surprise_correct, target_surprise_error,...
    ~,all_pre_sound,all_p_pre_sound,all_sound_correct_error,~,all_norm_spike,Brainarea)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

right_neuron = find(all_sound_correct_error(:,1) > 0);
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

[L_sig_PostLow_left, L_sig_PostHigh_left, L_sig_PostLow_right, L_sig_PostHigh_right] = ...
    test_surprise_sound_choice(all_norm_spike, left_sig_neuron, target_surprise_both);
%Use only on high sig neuron
[R_sig_PostLow_left, R_sig_PostHigh_left, R_sig_PostLow_right, R_sig_PostHigh_right] = ...
    test_surprise_sound_choice(all_norm_spike, right_sig_neuron, target_surprise_both);

[L_sig_PostLow_leftcorrect, L_sig_PostHigh_leftcorrect, L_sig_PostLow_rightcorrect, L_sig_PostHigh_rightcorrect] = ...
    test_surprise_sound_choice(all_norm_spike, left_sig_neuron, target_surprise_correct);
%Use only on high sig neuron
[R_sig_PostLow_leftcorrect, R_sig_PostHigh_leftcorrect, R_sig_PostLow_rightcorrect, R_sig_PostHigh_rightcorrect] = ...
    test_surprise_sound_choice(all_norm_spike, right_sig_neuron, target_surprise_correct);

[L_sig_PostLow_lefterror, L_sig_PostHigh_lefterror, L_sig_PostLow_righterror, L_sig_PostHigh_righterror] = ...
    test_surprise_sound_choice(all_norm_spike, left_sig_neuron, target_surprise_error);
%Use only on high sig neuron
[R_sig_PostLow_lefterror, R_sig_PostHigh_lefterror, R_sig_PostLow_righterror, R_sig_PostHigh_righterror] = ...
    test_surprise_sound_choice(all_norm_spike, right_sig_neuron, target_surprise_error);


%Activity is already flipped

Repeat_prefer = [L_sig_PostLow_left; R_sig_PostHigh_right];
Repeat_nonprefer = [L_sig_PostHigh_right; R_sig_PostLow_left];
Switch_prefer = [L_sig_PostHigh_left; R_sig_PostLow_right];
Switch_nonprefer = [L_sig_PostLow_right; R_sig_PostHigh_left];

Repeat_prefer_correct = [L_sig_PostLow_leftcorrect; R_sig_PostHigh_rightcorrect];
Repeat_nonprefer_correct = [L_sig_PostHigh_rightcorrect; R_sig_PostLow_leftcorrect];
Switch_prefer_correct = [L_sig_PostHigh_leftcorrect; R_sig_PostLow_rightcorrect];
Switch_nonprefer_correct = [L_sig_PostLow_rightcorrect; R_sig_PostHigh_leftcorrect];

Repeat_prefer_error = [L_sig_PostLow_lefterror; R_sig_PostHigh_righterror];
Repeat_nonprefer_error = [L_sig_PostHigh_righterror; R_sig_PostLow_lefterror];
Switch_prefer_error = [L_sig_PostHigh_lefterror; R_sig_PostLow_righterror];
Switch_nonprefer_error = [L_sig_PostLow_righterror; R_sig_PostHigh_lefterror];



p_prefer_correct = signrank(Repeat_prefer_correct,Switch_prefer_correct);
p_nonprefer_correct = signrank(Repeat_nonprefer_correct,Switch_nonprefer_correct);



p_prefer_error = signrank(Repeat_prefer_error,Switch_prefer_error);
p_nonprefer_error = signrank(Repeat_nonprefer_error,Switch_nonprefer_error);

prefer_sabun = Repeat_prefer - Switch_prefer;
nonprefer_sabun = Repeat_nonprefer - Switch_nonprefer;

prefer_sabun_correct = Repeat_prefer_correct - Switch_prefer_correct;
prefer_sabun_error = Repeat_prefer_error - Switch_prefer_error;

nonprefer_sabun_correct = Repeat_nonprefer_correct - Switch_nonprefer_correct;
nonprefer_sabun_error = Repeat_nonprefer_error - Switch_nonprefer_error;

p_sabun = [p_prefer_correct, p_nonprefer_correct, p_prefer_error, p_nonprefer_error];



temp1 = length(prefer_sabun_correct);
temp2 = length(nonprefer_sabun_correct);

if temp1 ~= temp2
    hoge
else
    length_neuron = temp1;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PostLow_left, PostHigh_left, PostLow_right, PostHigh_right] = ...
    test_surprise_sound_choice(all_norm_spike, low_sig_neuron, target_surprise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PostLow_left = all_norm_spike(low_sig_neuron,target_surprise(1)); %PostLowCorrect_left
PostHigh_left = all_norm_spike(low_sig_neuron,target_surprise(2)); %PostHighCorrect_left
PostLow_right = all_norm_spike(low_sig_neuron,target_surprise(3)); %PostLowCorrect_right
PostHigh_right = all_norm_spike(low_sig_neuron,target_surprise(4)); %PostHighCorrect_right

function plot_surprise_choice_prefer(Repeat_prefer, Switch_prefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(Repeat_prefer, Switch_prefer, '.','color', plot_color) %same choice, X_repeat or Y_zigzag
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])


function plot_surprise_choice_nonprefer(Switch_nonprefer, Repeat_nonprefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(Switch_nonprefer, Repeat_nonprefer, '.','color', plot_color)  %same choice, X_zigzag or Y_repeat
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])

function [sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category,tone_correct_error,p_tone_correct_error] = ...
    HMM_sound_trans20240705_surprise_during_depth(pathname, kaiseki_number,depth_def)

switch nargin
    case 0
        pathname = pwd;
    case 3
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

temp = dir('HMM_spike_count_neurons_sound_trans_20240326*');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADD depth definition
use_sig_neuron = intersect(use_sig_neuron,depth_neuron);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% tone_correct_error(:,1) = use_neuron_index(use_sig_neuron,12); %Correct Low high
% tone_correct_error(:,2) = use_neuron_index(use_sig_neuron,13); %Error Low high
% p_tone_correct_error(:,1) = use_p_index(use_sig_neuron,12);
% p_tone_correct_error(:,2) = use_p_index(use_sig_neuron,13);

tone_correct_error = [];
p_tone_correct_error = [];

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

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_neuron = get_sig_neuron_time_window(p_task,new_p_thre,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_p_sound = min(p_task(:,window),[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return




