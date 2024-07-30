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
function analysis_dir = repeat_overtrained_20240420

analysis_dir(1).matrix = repeat_a17_20240420;
analysis_dir(2).matrix = repeat_w03_20240420;
analysis_dir(3).matrix = repeat_w06_20240420;
analysis_dir(4).matrix = repeat_w22_20240420;
analysis_dir(5).matrix = repeat_w31_20240420;

%analysis_dir = [analysis_dir1;analysis_dir2;analysis_dir3;analysis_dir4;analysis_dir5];

