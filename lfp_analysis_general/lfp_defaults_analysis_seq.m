
function anaparams=lfp_defaults_analysis_seq(analysis)

switch analysis
    
    case 'fft_time'
        anaparams.Fs = 1000;
        anaparams.analysis=analysis;
        anaparams.window_size=250;
        anaparams.step=10;
        anaparams.fpass=[1 100];
        anaparams.trialave=0;
        anaparams.norm=1;
        anaparams.freqs=[1:1.5:100];
%         num_freqs=100;
%         min_freq = 2^(6/8);
%         max_freq = 2^(54/8);
%         anaparams.freqs = logspace(log10(min_freq),log10(max_freq),num_freqs);

    case 'wavelet'
        anaparams.Fs = 1000;
        anaparams.analysis=analysis;
        
        num_freqs=100;
        min_freq = 2^(6/8);
        max_freq = 2^(54/8);
        
        
        anaparams.freqs = logspace(log10(min_freq),log10(max_freq),num_freqs);
        anaparams.Zscoringbaseline=1;
        % remove mean signal from raw data before doing wavelet
        % transformation?
        anaparams.rmmean=0;
        
        
        

        
        
end