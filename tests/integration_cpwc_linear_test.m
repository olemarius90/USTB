classdef integration_cpwc_linear_test < matlab.unittest.TestCase

    methods (Test)
        function test_cpwc_pipeline_produces_image(testCase)
            pha = uff.phantom();
            pha.sound_speed = 1540;
            pha.points = [0, 0, 20e-3, 1];

            prb = uff.linear_array();
            prb.N = 128;
            prb.pitch = 300e-6;
            prb.element_width = 270e-6;
            prb.element_height = 5000e-6;

            pul = uff.pulse();
            pul.center_frequency = 5.2e6;
            pul.fractional_bandwidth = 0.6;

            N = 11;
            angles = linspace(-0.3, 0.3, N);
            seq = uff.wave();
            for n = 1:N
                seq(n) = uff.wave();
                seq(n).wavefront = uff.wavefront.plane;
                seq(n).source.azimuth = angles(n);
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

            scan = uff.linear_scan('x_axis', linspace(-10e-3, 10e-3, 128).', ...
                                   'z_axis', linspace(5e-3, 35e-3, 128).');

            mid = midprocess.das();
            mid.code = code.matlab;
            mid.dimension = dimension.both;
            mid.channel_data = channel_data;
            mid.scan = scan;

            mid.receive_apodization.window = uff.window.hanning;
            mid.receive_apodization.f_number = 1.7;

            mid.transmit_apodization.window = uff.window.hanning;
            mid.transmit_apodization.f_number = 1.7;

            b_data = mid.go();

            testCase.verifyEqual(numel(b_data.data(:)), 128*128);
            testCase.verifyTrue(all(isfinite(b_data.data(:))));

            [~, idx] = max(abs(b_data.data));
            peak_x = scan.x(idx);
            peak_z = scan.z(idx);
            testCase.verifyEqual(peak_x, 0, 'AbsTol', 2e-3);
            testCase.verifyEqual(peak_z, 20e-3, 'AbsTol', 2e-3);
        end

        function test_cpwc_output_has_expected_shape(testCase)
            pha = uff.phantom();
            pha.sound_speed = 1540;
            pha.points = [0, 0, 15e-3, 1];

            prb = uff.linear_array();
            prb.N = 64;
            prb.pitch = 300e-6;
            prb.element_width = 270e-6;
            prb.element_height = 5000e-6;

            pul = uff.pulse();
            pul.center_frequency = 5.2e6;
            pul.fractional_bandwidth = 0.6;

            N = 5;
            angles = linspace(-0.2, 0.2, N);
            seq = uff.wave();
            for n = 1:N
                seq(n) = uff.wave();
                seq(n).wavefront = uff.wavefront.plane;
                seq(n).source.azimuth = angles(n);
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

            Nx = 64; Nz = 64;
            scan = uff.linear_scan('x_axis', linspace(-5e-3, 5e-3, Nx).', ...
                                   'z_axis', linspace(10e-3, 20e-3, Nz).');

            mid = midprocess.das();
            mid.code = code.matlab;
            mid.dimension = dimension.both;
            mid.channel_data = channel_data;
            mid.scan = scan;
            mid.receive_apodization.window = uff.window.boxcar;
            mid.receive_apodization.f_number = 1.7;
            mid.transmit_apodization.window = uff.window.boxcar;
            mid.transmit_apodization.f_number = 1.7;

            b_data = mid.go();

            testCase.verifyClass(b_data, ?uff.beamformed_data);
            testCase.verifyEqual(size(b_data.data, 1), Nx * Nz);
            testCase.verifyTrue(max(abs(b_data.data(:))) > 0);
        end
    end

end
