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
function analysis_dir = zigzag_overtrain_20240420

analysis_dir(1).matrix = zigzag_w02_20240420;
analysis_dir(2).matrix = zigzag_w13_20240420;
analysis_dir(3).matrix = zigzag_w14_20240420;
analysis_dir(4).matrix = zigzag_w23_20240420;
analysis_dir(5).matrix = zigzag_w25_20240420;
analysis_dir(6).matrix = zigzag_w27_20240420;
analysis_dir(7).matrix = zigzag_w28_20240420;
analysis_dir(8).matrix = zigzag_w29_20240420;
analysis_dir(9).matrix = zigzag_w30_20240420;

%analysis_dir = [analysis_dir1;analysis_dir2;analysis_dir3;analysis_dir4;analysis_dir5];

