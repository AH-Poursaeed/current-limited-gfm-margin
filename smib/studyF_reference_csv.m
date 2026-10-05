function p = studyF_reference_csv()
% Path of the Table XII results in data/, which the fine sweep is compared
% with.

p = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'smib', ...
             'studyC_aggregation_error.csv');
if ~isfile(p)
    error('studyF_reference_csv:missing', 'cannot find %s', p);
end
end
