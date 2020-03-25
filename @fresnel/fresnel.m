classdef fresnel < handle
%fresnel   fresnel definition
%
%   See also PULSE, BEAM, PHANTOM, PROBE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/02/24 $

    %% public properties
    properties  (Access = public)
        phantom             % phantom class
        pulse               % pulse class
        probe               % probe class
        sequence            % collection of wave classes
        sampling_frequency  % sampling frequency [Hz]
        PRF                 % pulse repetition frequency [Hz]
    end
    
    %% dependent properties
    properties  (Dependent)   
        N_elements         % number of elements in the probe
        N_points           % number of points in the phantom
        N_waves            % number of waves
        N_events           % number of events (waves*frames)
        N_frames           % number of frames
    end
    
    %% private properties
    properties  (Access = private)   
        version='v1.0.7';  % fresnel version
    end
    
    %% constructor
    methods (Access = public)
        function h=fresnel()
            %fresnel   Constructor of fresnel class
            %
            %   Syntax:
            %   h = fresnel()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE                      
            
        end
    end
    
    %% set methods
    methods  
        function out_dataset=go(h)
            disp(sprintf('USTB''s Fresnel impulse response simulator (%s)',h.version));
            disp('---------------------------------------------------------------');
            
            %% checking we have all we need
            assert(~isempty(h.probe),'The PROBE parameter is not set.');
            assert(~isempty(h.phantom),'The PHANTOM parameter is not set.');
            assert(length(h.phantom)==1,'Only static phantoms are supported');
            assert(~isempty(h.pulse),'The PULSE parameter is not set.');
            assert(~isempty(h.sequence),'The SEQUENCE parameter is not set.');
            assert(~isempty(h.sampling_frequency),'The SAMPLING_FREQUENCY parameter is not set.');
            assert(length(h.phantom)==1,'Only single frames are supported');

            % checking number of elements
            assert(h.probe.N_elements==h.sequence(1).N_elements,'Mismatch in the number of elements in probe and the size of delay and apodization vectors in beam');
           
            %% unwrapping the signal
            focusing_delay=zeros(1,h.N_elements,h.N_waves);
            apodization=zeros(1,h.N_elements,h.N_waves);
            for n_wave=1:h.N_waves 
                focusing_delay(1,:,n_wave)=h.sequence(n_wave).delay_values;
                apodization(1,:,n_wave)=h.sequence(n_wave).apodization_values;
            end
            
            c0 = h.phantom.sound_speed;
            f0 = h.pulse.center_frequency;
            w0 = 2*pi*f0;
            k0 = w0/c0;
            fs = h.sampling_frequency;
            bw = h.pulse.fractional_bandwidth;
            
            %% minimum distance for including geometric dispersion
            delta0=4*pi*0.1e-3;

            % save the data into a CHANNEL_DATA structure
            out_dataset=uff.channel_data();
            out_dataset.probe=h.probe();
            out_dataset.pulse=h.pulse();
            out_dataset.phantom=h.phantom();
            out_dataset.sequence=h.sequence();
            out_dataset.sampling_frequency=h.sampling_frequency();
            out_dataset.sound_speed=h.phantom.sound_speed;
          
            
            % computing geometry relations to the point
            distance  = sqrt((h.phantom.x.'-h.probe.x).^2+(h.phantom.y.'-h.probe.y).^2+(h.phantom.z.'-h.probe.z).^2);
            theta = atan2(h.phantom.x.'-h.probe.x, h.phantom.z.'-h.probe.z)-h.probe.theta;
            phi = asin((h.phantom.y.'-h.probe.y)./distance)-h.probe.phi;
            
            % directivity between probe and the point
            directivity = sinc(k0*h.probe.width/2/pi.*tan(theta)).*sinc(k0*h.probe.height/2/pi.*tan(phi)./cos(theta));
            
            % delay between probe and the point
            propagation_delay = permute(distance/c0, [3,1,2]);
            
            % attenuation (absorption & geometrical dispersion)
            attenuation = permute(10.^(-h.phantom.alpha*(distance*1e2)*(f0*1e-6)).*directivity.*delta0./(4*pi*distance), [3,1,2]);
            
            min_range = min(distance(:));
            max_range = max(distance(:));
            min_delay = min(focusing_delay(:));
            max_delay = max(focusing_delay(:));
            
            time_1w = ((min_range/c0 - 8/f0/bw + min_delay):1/fs:(max_range/c0 + 8/f0/bw + max_delay)).';                                                  % time vector [s]
            time_2w = ((2*min_range/c0 - 8/f0/bw + min_delay):1/fs:(2*max_range/c0 + 8/f0/bw + max_delay)).';                                               % time vector [s]
            N_samples = length(time_2w);                                                                              % number of time samples
            
            % check if there wave delays are imposed in the sequence
            % definition
            if any(abs([h.sequence.delay])>0)
                out_dataset.initial_time = 0;
                wave_delays=true;
            else
                out_dataset.initial_time = time_2w(1);
                wave_delays  = false;
            end
            
            out_dataset.data = zeros([N_samples,h.N_elements,h.N_waves]);
            out_dataset.PRF = h.PRF;
            
            F = griddedInterpolant();
            F.Method = 'linear';
            F.ExtrapolationMethod = 'none';
            F.GridVectors = {time_1w};

            % the wave loop
            receive_signal = zeros([N_samples,h.N_elements,h.phantom.N_points]);
            for n_wave=1:h.N_waves
                
                % computing the transmit signal
                transmit_delay = time_1w - (propagation_delay + h.sequence(n_wave).delay_values.');
                transmit_signal = sum(h.pulse.signal(transmit_delay).*apodization(:,:,n_wave).*attenuation, 2);
                
                % computing the receive signal
                receive_delay = time_2w - propagation_delay;
                if wave_delays
                    receive_delay = receive_delay + h.sequence(n_wave).delay-time_2w(1);
                end
                
                for n_point = 1:h.phantom.N_points
                    F.Values = transmit_signal(:,1,n_point);
                    receive_signal(:,:,n_point) = F(receive_delay(:,:,n_point));
                end
                out_dataset.data(:,:,n_wave) = sum(receive_signal.*attenuation, 3, 'omitnan');
                
                % computing second order scattering
                %                         extra_distance = sqrt(sum((bsxfun(@minus,current_phantom.points([1:n_p-1 n_p+1:h.N_points],1:3),current_phantom.points(n_p,1:3))).^2,2));
                %                         extra_delay = extra_distance/current_phantom.sound_speed;
                %                         extra_attenuation= 10.^(-current_phantom.alpha*(extra_distance/1e-2)*(h.pulse.center_frequency/1e6)).*delta0./(4*pi*extra_distance);
                %                         for nnp=1:length(extra_distance)
                %                             h.reverb(:,:,n_w,n_f)=h.reverb(:,:,n_w,n_f)+bsxfun(@times,interp1(time_1w,transmit_signal,delayed_time-extra_delay(nnp),'linear',0),attenuation.*extra_attenuation(nnp)).';
                %                         end
            end
        end
    end
    
    %% set methods
    methods  
        function set.phantom(h,in_phantom)
            assert(isa(in_phantom,'uff.phantom'), 'The phantom is not a PHANTOM class. Check HELP PHANTOM.');
            h.phantom=in_phantom;
        end
        function set.pulse(h,in_pulse)
            assert(isa(in_pulse,'uff.pulse'), 'The pulse is not a PULSE class. Check HELP PULSE.');
            h.pulse=in_pulse;
        end
        function set.probe(h,in_probe)
            assert(isa(in_probe,'uff.probe'), 'The probe is not a PROBE class. Check HELP PROBE.');
            h.probe=in_probe;
        end
        function set.sequence(h,in_sequence)
            assert(isa(in_sequence,'uff.wave'), 'The sequence members are not a WAVE class. Check HELP WAVE.');
            h.sequence=in_sequence;
        end
        function set.sampling_frequency(h,in_sampling_frequency)
            assert(numel(in_sampling_frequency)==1, 'The sampling frequency must be a scalar');
            h.sampling_frequency=in_sampling_frequency;
        end       
        function set.PRF(h,in_PRF)
            assert(numel(in_PRF)==1, 'The PRF must be a scalar');
            h.PRF=in_PRF;
        end       
    end
    
    %% get methods
    methods  
        function value=get.N_elements(h)
            value=h.probe.N_elements;
        end
        function value=get.N_points(h)
            value=h.phantom.N_points;
        end
        function value=get.N_waves(h)
            value=numel(h.sequence);
        end
        function value=get.N_events(h)
            value=length(h.phantom);
        end
        function value=get.N_frames(h)
            value=ceil(h.N_events/h.N_waves);
        end
    end
    
end