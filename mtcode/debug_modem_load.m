% DEBUG SCRIPT for load_data_modem parsing issues
% Run in MATLAB from the mtcode folder:
%   cd('/home/ramonvanoli/MasterThesis/mtcode')
%   debug_modem_load('path/to/your/datafile.dat')

function debug_modem_load(filename)

if nargin < 1
    [f, p] = uigetfile('*.dat', 'Pick ModEM data file');
    if isequal(f, 0), return; end
    filename = fullfile(p, f);
end

fprintf('\n=== Debugging: %s ===\n\n', filename);

% 1) Count lines in file by type
fid = fopen(filename, 'r');
n_comment = 0; n_header = 0; n_data = 0; n_other = 0;
data_lines = {};
while true
    line = fgetl(fid);
    if line == -1, break; end
    if isempty(line), continue; end
    tline = strtrim(line);
    if tline(1) == '#'
        n_comment = n_comment + 1;
    elseif tline(1) == '>'
        n_header = n_header + 1;
    elseif ~isempty(sscanf(line, '%f', 1))
        n_data = n_data + 1;
        data_lines{end+1} = line; %#ok<AGROW>
    else
        n_other = n_other + 1;
        fprintf('  Unclassified line: %s\n', line(1:min(80,end)));
    end
end
fclose(fid);

fprintf('Line counts:\n');
fprintf('  # comment lines : %d\n', n_comment);
fprintf('  > header lines  : %d (expect multiple of 6)\n', n_header);
fprintf('  data lines      : %d\n', n_data);
fprintf('  other lines     : %d\n\n', n_other);

% 2) Old vs new data-line detection on first 4 data lines
if ~isempty(data_lines)
    fprintf('First-char check on first 4 data lines:\n');
    for k = 1:min(4, numel(data_lines))
        line = data_lines{k};
        old_ok = ~isempty(str2num(line(1))); %#ok<ST2NM>
        new_ok = ~isempty(sscanf(line, '%f', 1));
        data = textscan(line, '%f %s %f %f %f %f %f %s %f %f %f %f %f %f %f');
        comp = data{8}{1};
        fprintf('  line %d: first char = ''%c'' | old check = %d | new check = %d | component = %s\n', ...
            k, line(1), old_ok, new_ok, comp);
    end
    fprintf('\n');
end

% 3) Component breakdown in file
comps = {'ZXX','ZXY','ZYX','ZYY','TX','TY'};
comp_counts = zeros(1, numel(comps));
for k = 1:numel(data_lines)
    data = textscan(data_lines{k}, '%f %s %f %f %f %f %f %s %f %f %f %f %f %f %f');
    comp = data{8}{1};
    idx = find(strcmp(comp, comps), 1);
    if ~isempty(idx)
        comp_counts(idx) = comp_counts(idx) + 1;
    else
        fprintf('  Unknown component on data line %d: %s\n', k, comp);
    end
end
fprintf('Components found in file:\n');
for k = 1:numel(comps)
    fprintf('  %s : %d lines\n', comps{k}, comp_counts(k));
end
fprintf('\n');

% 4) Load with fixed function and inspect result
addpath('Load_Functions');
d = load_data_modem(filename);

fprintf('After load_data_modem:\n');
fprintf('  ns = %d, nf = %d, nlines parsed = %d\n', d.ns, d.nf, n_data);
fprintf('  d.responses = '); disp(d.responses);

for is = 1:min(d.ns, 2)
    fprintf('\n  Site %d (%s): non-NaN Z components per column\n', is, d.site{is});
    for ic = 1:4
        n_ok = sum(~isnan(real(squeeze(d.Z(:,ic,is)))));
        fprintf('    %s (col %d): %d / %d periods have data\n', comps{ic}, ic, n_ok, d.nf);
    end
    fprintf('  Z tensor for site %s (first 5 periods):\n', d.site{is});
    disp(squeeze(d.Z(1:min(5,d.nf),:,is)));
end

fprintf('\nDone. If file has all 4 components but load shows only 2,\n');
fprintf('re-run after the load_data_modem fix (pull latest version).\n\n');

end
