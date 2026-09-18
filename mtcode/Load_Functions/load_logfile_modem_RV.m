function invlog = load_logfile_modem_RV(filename)
% Robust parser for ModEM NLCG .log files
% Parses "with:" lines which contain the per-iteration summary

fid = fopen(filename, 'r');
if fid < 0
    error('Cannot open file: %s', filename);
end

iter_num = [];
f_vals   = [];
m2_vals  = [];
rms_vals = [];
lam_vals = [];
alp_vals = [];

iter_counter = 0;

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    
    % Only parse the final "with:" line per iteration (the confirmed result)
    if ~isempty(regexp(line, 'with:', 'once'))
        iter_counter = iter_counter + 1;
        
        % Extract each named value using regex
        f   = extract_val(line, 'f');
        m2  = extract_val(line, 'm2');
        rms = extract_val(line, 'rms');
        lam = extract_val(line, 'lambda');
        alp = extract_val(line, 'alpha');  % NaN if ****
        
        iter_num(end+1) = iter_counter;
        f_vals(end+1)   = f;
        m2_vals(end+1)  = m2;
        rms_vals(end+1) = rms;
        lam_vals(end+1) = lam;
        alp_vals(end+1) = alp;
    end
end
fclose(fid);

invlog.niter  = numel(iter_num);
invlog.rms    = rms_vals(:);
invlog.lambda = lam_vals(:);
invlog.alpha  = alp_vals(:);
invlog.m2     = m2_vals(:);
invlog.f      = f_vals(:);

fprintf('Parsed %d iterations from %s\n', invlog.niter, filename);
end

function val = extract_val(line, name)
% Extract a named value like "rms=  61.398872" or "alpha=************"
    pattern = [name, '=\s*([\d.E+\-]+|\*+)'];
    tok = regexp(line, pattern, 'tokens', 'once');
    if isempty(tok) || ~isempty(regexp(tok{1}, '\*', 'once'))
        val = NaN;
    else
        val = str2double(tok{1});
    end
end