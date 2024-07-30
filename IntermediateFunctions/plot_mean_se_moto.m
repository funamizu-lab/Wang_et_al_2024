%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mean_trace,std_trace,se_trace] = plot_mean_se_moto(trace, trace_color,std_se)
%std_se: 0 -> plot only mean
%      : 1 -> plot with std
%      : 2 -> plot with se
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%trace
%y_axis: trials
%x_axis: time

[trial,time_trace] = size(trace);

%Extract nan trial
nan_check = mean(trace,2);
nan_check = find(isnan(nan_check) == 0); %detect_non_nan
if trial ~= length(nan_check)
    disp('detect nan trial')
    [trial, length(nan_check)]
end
trace = trace(nan_check,:);

mean_trace = mean(trace,1);
std_trace  = std(trace,1);
se_trace = std_trace ./ (sqrt(trial));

std_plus  = mean_trace + std_trace;
std_minus = mean_trace - std_trace;
se_plus  = mean_trace + se_trace;
se_minus = mean_trace - se_trace;

temp_x1 = [1:time_trace];
temp_x2 = [time_trace:-1:1];
temp_x = [temp_x1, temp_x2];

if std_se == 0 %mean only
plot(mean_trace,'color',trace_color,'LineWidth',1)
box off
elseif std_se == 1 %std
fill(temp_x,[std_plus, fliplr(std_minus)],trace_color,'edgecolor','none')
alpha(0.1)
hold on
plot(mean_trace,'color',trace_color,'LineWidth',1)
box off

elseif std_se == 2 %se
fill(temp_x,[se_plus, fliplr(se_minus)],trace_color,'edgecolor','none')
alpha(0.1)
hold on
plot(mean_trace,'color',trace_color,'LineWidth',1)
box off

end