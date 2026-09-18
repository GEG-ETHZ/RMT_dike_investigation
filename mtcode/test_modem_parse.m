% Quick test of load_data_modem parsing logic
testfile = 'test_modem_sample.dat';
fid = fopen(testfile,'w');
fprintf(fid,'# written by Matlab function write_data_modem\n');
fprintf(fid,'# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n');
fprintf(fid,'> Full_Impedance\n');
fprintf(fid,'> exp(+i\\omega t)\n');
fprintf(fid,'> [mV/km]/[nT]\n');
fprintf(fid,'> 0.00\n');
fprintf(fid,'> 51.751 4.633\n');
fprintf(fid,'> 2 1\n');
lines = {
    ' 1.152339E-05             1201_edited_01-Jun-2   51.750    4.634        -27.000         21.000          0.000     ZXX   2.514005E-08   2.605264E-08   1.000000E+13    0.000   90.000    0.000   90.000'
    ' 1.152339E-05             1201_edited_01-Jun-2   51.750    4.634        -27.000         21.000          0.000     ZXY   4.649860E+03   4.818650E+03   1.000000E+13    0.000   90.000    0.000   90.000'
    ' 1.152339E-05             1201_edited_01-Jun-2   51.750    4.634        -27.000         21.000          0.000     ZYX  -4.649860E+03  -4.818650E+03   1.000000E+13    0.000   90.000    0.000   90.000'
    ' 1.152339E-05             1201_edited_01-Jun-2   51.750    4.634        -27.000         21.000          0.000     ZYY  -2.514005E-08  -2.605264E-08   1.000000E+13    0.000   90.000    0.000   90.000'
};
for k = 1:numel(lines)
    fprintf(fid,'%s\n', lines{k});
end
fclose(fid);

% Test first-char check
for k = 1:numel(lines)
    line = lines{k};
    fprintf('Line %d: line(1)=%c str2num(line(1))=%s sscanf=%g\n', k, line(1), mat2str(str2num(line(1))), sscanf(line,'%f',1));
end

addpath('Load_Functions');
d = load_data_modem(testfile);
fprintf('Loaded Z for site 1:\n');
disp(squeeze(d.Z(:,:,1)));
