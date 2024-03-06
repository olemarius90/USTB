classdef refocus < preprocess
    %REFOCUS
    %
    %   authors: Anders E. Vrålstad <anders.e.vralstad@ntnu.no>
    %   
    %   Code adapted from 
    %   github.com/nbottenus/REFoCUS
    %   
    %   $Last updated: 2024/03/06$
    
    %% constructor
    methods (Access = public)
        function h = refocus()
            h.name='REFoCUS implemented in MATLAB';
            h.reference='www.github.com/nbottenus/REFoCUS';
            h.implemented_by={'Anders E. Vrålstad <anders.e.vralstad@ntnu.no>'};
            h.version='v1.0.0'; 
        end
    end
    
    properties (Access = public)
        post_pad_samples = 200;
        use_filter = false;
        regularization = @Hinv_adjoint;
        decode_parameter = [];
    end
    
    properties (Dependent)
                                                                 
    end
    
    methods (Access = public)
        function output=go(h)  
            % Check if we can skip calculation
            if h.check_hash()
                output= h.output;
                return;
            end

            % If IQ
            if abs(h.input.modulation_frequency)>eps
                error('Only implemented for RF channel data, not IQ.')
            end
            
            N_channels = h.input.N_channels;
            N_waves = h.input.N_waves;
            N_frames = h.input.N_frames;
            
            tx_delays = zeros(N_channels,N_waves);
            tx_apod = zeros(N_channels,N_waves);
            for wave = 1:N_waves
                tx_delays(:,wave) = h.input.sequence(wave).delay_values - h.input.sequence(wave).delay;
                tx_apod(:,wave) = h.input.sequence(wave).apodization_values;
            end
         
            rxdata_multiTx = padarray(h.input.data,h.post_pad_samples,'post');
            normalized_rxdata_multiTx = double(rxdata_multiTx / max(rxdata_multiTx(:)));
            N_samples_output = size(normalized_rxdata_multiTx,1);

            %%
            full_synth_data = zeros(N_samples_output, N_channels, N_channels,N_frames);
            for fr = 1:N_frames
            % Decode Multistatic Data Using REFoCUS
                full_synth_data(:,:,:,fr) = refocus_decode(normalized_rxdata_multiTx(:,:,:,fr),tx_delays.'*h.input.sampling_frequency,...
                'fun',h.regularization,'apod',tx_apod.','param',h.decode_parameter);
            end
            
            % Passband Filter Channel Data
            if h.use_filter
                N = 10; % Filter Order
                Wn = [0.05,0.3];
                [b, a] = butter(N, Wn); % Filter
                full_synth_data = filtfilt(b, a, double(full_synth_data));
            end

            % Create output channel data object
            h.output = uff.channel_data();
            h.output.initial_time = h.input.initial_time;
            h.output.modulation_frequency = h.input.modulation_frequency;
            h.output.sampling_frequency = h.input.sampling_frequency;
            h.output.data = single(full_synth_data);
            h.output.sound_speed = h.input.sound_speed;
            h.output.probe = uff.probe; h.output.probe.geometry = h.input.probe.geometry;
            h.output.sequence = uff.wave;

            for wave = 1:N_channels
                h.output.sequence(wave) = uff.wave;
                h.output.sequence(wave).probe = uff.probe; 
                h.output.sequence(wave).probe.geometry = h.input.probe.geometry; 
                h.output.sequence(wave).sound_speed = h.input.sound_speed;
                h.output.sequence(wave).source = uff.point('xyz', h.input.probe.geometry(wave,1:3));
                h.output.sequence(wave).delay = sqrt(sum(h.output.sequence(wave).source.xyz.^2)) / h.output.sequence(wave).sound_speed; 
                h.output.sequence(wave).origin = uff.point('xyz', h.input.probe.geometry(wave,1:3));
            end
            
            % Pass reference
            output = h.output;

            % Update hash
            h.save_hash();
        end
    end
end