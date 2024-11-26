function analytic_signal = lfp_wavelet_decomposition_seq(samples, srate, frex, ...
                                     numcycles,num_frex)

%% function analytic_signal = wavelet_decomposition(samples, srate, frex, ...
%                                     numcycles)
% samples are either expected to a row vector, or a 2D array with
% samples(segment, samplepoint);

nsegs = size(samples,1);
nsamples = size(samples,2);
if nsegs > 1
    samples = reshape(samples',1,nsegs*nsamples);
end

time  = -2:1/srate:2; % time, from -1 to 1 second in steps of 1/sampling-rate

if ~exist('frex', 'var') || isempty(frex)
    min_freq = 2^(6/8);
    max_freq = 2^(30/8);
    frex = logspace(log10(min_freq),log10(max_freq),num_frex);
end

num_frex = numel(frex);

if ~exist('numcycles', 'var') || isempty(numcycles)
    numcycles = 4; % sets the trade-off between precision in temporal
                   % and frequency domain
end
s = numcycles./(2*pi*frex);

tic
% definte convolution parameters 
n_wavelet = length(time);
n_data = numel(samples);
n_convolution = n_wavelet+n_data-1;
n_conv_pow2 = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% init output
if nsegs > 1
    analytic_signal = zeros(num_frex, nsegs, nsamples);
else    
    analytic_signal = zeros(num_frex,n_data);
end

% get fft of signal
eegfft = fft(samples, n_conv_pow2);
for fi = 1:numel(frex)

    wavelet = sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2)));
    wavelet = fft(wavelet , n_conv_pow2);    
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end- ...
                      half_of_wavelet_size);
    
    % save output
    if nsegs == 1
        analytic_signal(fi,1:n_data) = eegconv;
    else
       
        analytic_signal(fi,1:nsegs,1:nsamples) = reshape(eegconv,nsamples,nsegs)';
    end
    

end

toc