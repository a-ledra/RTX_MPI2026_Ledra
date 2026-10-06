%{
main_faithful_select.m

One main code that finds the faithful (feasible) set and then selects a
codebook from it with whichever algorithm you choose.  Set the parameters
and METHOD below, then run.

Algorithms (all in this folder):
  'greedy_first'      Dr. Edwards, Greedy First      select_greedy_first.m
  'greedy_last'       Dr. Edwards, Greedy Last       select_greedy_last.m
  'maxgap_backward'   Lawan, Backward Max-Gap        select_maxgap_backward.m
  'maxgap_forward'    Lawan, Forward Max-Gap         select_maxgap_forward.m
  'greedy_rows'       Adrianna, greedy farthest-pt   select_greedy_rows.m
  'quantile_rows'     Adrianna, quantile spacing     select_quantile_rows.m
  'all'               run all six and compare them

Faithful-set code (getfaithful, bestmatch, dae_DAD) is from daeintraplot.m

%}

clear all

%% Choose the algorithm here
method = 'greedy_first';

%% Parameters
n = 512;          % Size of Hadamard matrix
omegaw = 10;      % Word transmission rate
fmax = 512;       % Desired analog bandwidth (Hz): also called omegaa
snr = -3;         % Signal-to-noise ratio (dB); Inf means no noise
maxrun = 50;      % Number of noisy runs (use 1 if snr = Inf)
w = 4;            % Word length (bits)
esize = 2^w;      % Size of codebook

methodList = {'greedy_first','greedy_last','maxgap_backward', ...
              'maxgap_forward','greedy_rows','quantile_rows'};
if strcmpi(method,'all')
    methods = methodList;
else
    methods = {method};
end

%% Find the faithful set
% Compute the Walsh matrix of size n.
h = n*fwht(eye(n));

faithmat = cell(maxrun,1);
for run = 1:maxrun
    faithmat{run,1} = getfaithful(h,omegaw,fmax,snr);
end

% The feasible set is every row that is faithful in EVERY run.
onemat = vertcat(faithmat{:});
feasible = find(all(onemat==1, 1))';
fprintf('Feasible set: %d of %d rows (omegaw=%d, fmax=%d, snr=%g, %d runs)\n', ...
    length(feasible),n,omegaw,fmax,snr,maxrun);

%% Select a codebook with each requested algorithm
codebooks = NaN(esize,length(methods));
mindist = NaN(1,length(methods));
for k = 1:length(methods)
    codebooks(:,k) = select_codebook(feasible,esize,methods{k});
    mindist(k) = min(diff(codebooks(:,k)));
end

%% Report
results = array2table(codebooks,'VariableNames',methods);
disp('Selected codebook rows (Walsh scores):');
disp(results);

summary = table(string(methods(:)),mindist(:), ...
    'VariableNames',{'Method','MinDistance'});
disp('Minimum distance between neighboring codewords:');
disp(summary);

%%
function f = bestmatch(i,h,n,omegaw,fmax,snr)
% Computes the index of the closest transmitted encoding to the received
% encoding i, or 0 if the minimum is not unique.

fs_sym = n*omegaw;  % symbol/sample rate (Hz)
fs_analog = fs_sym; % analog sampling frequency (Hz)

% The input string is the ith row of the Walsh matrix:
x = h(i,:);
x_rec = dae_DAD(x,fs_analog,fs_sym,n,fmax,snr);

% Find the encoding with minimum Hamming distance from x_rec.
d = sum(h ~= x_rec, 2);
f = find(d == min(d));
% To be faithful, the minimum must be unique.
if length(f)~=1
    f = 0;
end

end

function x_rec = dae_DAD(x,fs_analog,fs_sym,n,fmax,snr)
% Passes x through the digital-analog-digital process: bandlimit, add
% noise, sample, and make a hard decision.  Written by David A. Edwards
% and Copilot.

% Create analog (bandlimited) waveform by resampling.
analog = resample(x,fs_analog,fs_sym);

% Enforce fmax with a 6th-order Butterworth lowpass.
Wn = fmax/(fs_analog/2);
if Wn < 1
    [b,a] = butter(6, Wn);
    analog = filtfilt(b,a,double(analog));
end

% Add white noise with SNR level snr.
analog = awgn(analog, snr, 'measured', 'db');

% Recover digital by sampling analog at symbol instants.
L = fs_analog / fs_sym;
if abs(L - round(L)) > 1e-10
    error('fs_analog must be an integer multiple of fs_sym for simple downsampling. Use resample for arbitrary ratios.');
end
L = round(L);
recovered_samples = analog(1:L:end);

% Decision device (hard decision to +-1)
x_rec = sign(recovered_samples);
x_rec(x_rec==0) = 1;

end

function faithful = getfaithful(h,omegaw,fmax,snr)
% Calculates the faithful set of Walsh matrix rows given a signal-noise
% ratio.

n = size(h,2);
% bestvec must be 1xn so the final calculation of faithful is correct.
bestvec = zeros(1,n);

for i = 1:n
    bestvec(i) = bestmatch(i,h,n,omegaw,fmax,snr);
end

% Faithful set: rows where the best match equals the input row.
faithful = (bestvec == (1:length(bestvec)));

end
