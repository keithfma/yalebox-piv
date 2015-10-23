% Exploration of peak finding using the fourier shift theorum.

load(data_file, 'xcr');

[rpk0, cpk0] = find(xcr == max(xcr(:)));

shift = FourierShift2D(xcr, [cpk0, rpk0]-1);

if abs(shift(1,1)-xcr(rpk0, cpk0)) < 10*eps
    fprintf('peak is shifted to 1,1\n');
else
    fprintf('something went wrong\n');
end

offset = 0.05;

shift = FourierShift2D(xcr, [cpk0, rpk0]-1+0.05);

if shift(1,1) > xcr(rpk0, cpk0)
    fprintf('shifted peak is better\n');
else
    fprintf('original peak is better\n');
end
