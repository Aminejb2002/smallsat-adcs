function rows = log_rows(a)
% To Workspace (Array) stores matrix signals as rows x cols x samples; return one row per sample.
if ndims(a) == 3
    rows = reshape(a, [], size(a, 3))';
else
    rows = a;
end
end
