function make_neuron_spike_file

%repeat
express_spike_data_each_region('repeat_OFC_20240428');
express_spike_data_each_region('repeat_Hippo_20240428');
express_spike_data_each_region('repeat_AC_20240428');

%alterante
express_spike_data_each_region('altern_OFC_20240428');
express_spike_data_each_region('altern_AC_20240428');
express_spike_data_each_region('altern_PPC_20240428');
express_spike_data_each_region('altern_M1_20240428');

end

function express_spike_data_each_region(folders)

[analysis_dir,~] = eval(folders);

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    process_make_neuron_file(analysis_dir{i});
end
end

function process_make_neuron_file(folder)

cd(folder)
mkdir spike_ch1/
load("spikedata.mat")
max_tif = length(spikedata.FireTiming);

for i = 1 : max_tif
    filename = append('task_spike_stripe20210520_', string(i), '.mat');
    spike_mark = zeros(1, spikedata.TimeLength);
    spike_mark(1, spikedata.FireTiming{i}) = 1;
    cd([folder, '/spike_ch1'])
    save(filename, 'spike_mark')
end

end