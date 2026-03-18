classdef integration_coherence_factor_test < matlab.unittest.TestCase

    methods (Test)
        function test_coherence_factor_pipeline(testCase)
            pha = uff.phantom();
            pha.sound_speed = 1540;
            pha.points = [0, 0, 40e-3, 1];

            prb = uff.linear_array();
            prb.N = 128;
            prb.pitch = 300e-6;
            prb.element_width = 270e-6;
            prb.element_height = 5000e-6;

            pul = uff.pulse();
            pul.center_frequency = 5.2e6;
            pul.fractional_bandwidth = 0.6;

            N = 5;
            angles = linspace(-0.3, 0.3, N);
            seq = uff.wave();
            for n = 1:N
                seq(n) = uff.wave();
                seq(n).source.azimuth = angles(n);
                seq(n).source.distance = Inf;
                seq(n).probe = prb;
                seq(n).sound_speed = pha.sound_speed;
            end

            sim = fresnel();
            sim.phantom = pha;
            sim.pulse = pul;
            sim.probe = prb;
            sim.sequence = seq;
            sim.sampling_frequency = 41.6e6;
            channel_data = sim.go();

            scan = uff.linear_scan();
            scan.x_axis = linspace(-3e-3, 3e-3, 100).';
            scan.z_axis = linspace(38e-3, 42e-3, 100).';

            mid = midprocess.das();
            mid.code = code.matlab;
            mid.dimension = dimension.none;
            mid.channel_data = channel_data;
            mid.scan = scan;
            mid.receive_apodization.window = uff.window.boxcar;
            mid.receive_apodization.f_number = 1.7;
            mid.transmit_apodization.window = uff.window.none;

            b_data = mid.go();

            cf = postprocess.coherence_factor();
            cf.input = b_data;
            cf.transmit_apodization = mid.transmit_apodization;
            cf.receive_apodization = mid.receive_apodization;
            b_data_cf = cf.go();

            testCase.verifyEqual(size(b_data_cf.data, 1), 100*100);
            testCase.verifyTrue(all(isfinite(b_data_cf.data(:))));

            cf_values = abs(cf.CF.data);
            testCase.verifyGreaterThanOrEqual(min(cf_values(:)), 0);
            testCase.verifyLessThanOrEqual(max(cf_values(:)), 1 + 1e-6);
        end

        function test_coherence_factor_high_at_target(testCase)
            pha = uff.phantom();
            pha.sound_speed = 1540;
            pha.points = [0, 0, 20e-3, 1];

            prb = uff.linear_array();
            prb.N = 64;
            prb.pitch = 300e-6;
            prb.element_width = 270e-6;
            prb.element_height = 5000e-6;

            pul = uff.pulse();
            pul.center_frequency = 5.2e6;
            pul.fractional_bandwidth = 0.6;

            N = 5;
            angles = linspace(-0.3, 0.3, N);
            seq = uff.wave();
            for n = 1:N
                seq(n) = uff.wave();
                seq(n).source.azimuth = angles(n);
                seq(n).source.distance = Inf;
                seq(n).probe = prb;
                seq(n).sound_speed = pha.sound_speed;
            end

            sim = fresnel();
            sim.phantom = pha;
            sim.pulse = pul;
            sim.probe = prb;
            sim.sequence = seq;
            sim.sampling_frequency = 41.6e6;
            channel_data = sim.go();

            scan = uff.linear_scan();
            scan.x_axis = linspace(-3e-3, 3e-3, 64).';
            scan.z_axis = linspace(18e-3, 22e-3, 64).';

            mid = midprocess.das();
            mid.code = code.matlab;
            mid.dimension = dimension.none;
            mid.channel_data = channel_data;
            mid.scan = scan;
            mid.receive_apodization.window = uff.window.boxcar;
            mid.receive_apodization.f_number = 1.7;
            mid.transmit_apodization.window = uff.window.none;

            b_data = mid.go();

            cf = postprocess.coherence_factor();
            cf.input = b_data;
            cf.transmit_apodization = mid.transmit_apodization;
            cf.receive_apodization = mid.receive_apodization;
            cf.go();

            cf_values = abs(cf.CF.data);
            [~, peak_idx] = max(cf_values(:));
            testCase.verifyGreaterThan(cf_values(peak_idx), 0.3);
        end
    end

end
