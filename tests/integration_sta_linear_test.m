classdef integration_sta_linear_test < matlab.unittest.TestCase

    methods (Test)
        function test_sta_pipeline_produces_image(testCase)
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

            N = prb.N_elements;
            seq = uff.wave();
            for n = 1:N
                seq(n) = uff.wave();
                seq(n).probe = prb;
                seq(n).wavefront = uff.wavefront.spherical;
                seq(n).source.xyz = [prb.x(n) prb.y(n) prb.z(n)];
                seq(n).sound_speed = pha.sound_speed;
            end

            sim = fresnel();
            sim.phantom = pha;
            sim.pulse = pul;
            sim.probe = prb;
            sim.sequence = seq;
            sim.sampling_frequency = 41.6e6;
            channel_data = sim.go();

            scan = uff.linear_scan('x_axis', linspace(-3e-3, 3e-3, 100).', ...
                                   'z_axis', linspace(18e-3, 22e-3, 100).');

            mid = midprocess.das();
            mid.code = code.matlab;
            mid.dimension = dimension.both;
            mid.channel_data = channel_data;
            mid.scan = scan;

            mid.receive_apodization.window = uff.window.tukey50;
            mid.receive_apodization.f_number = 1.7;

            mid.transmit_apodization.window = uff.window.tukey50;
            mid.transmit_apodization.f_number = 1.7;

            b_data = mid.go();

            testCase.verifyEqual(numel(b_data.data(:)), 100*100);
            testCase.verifyTrue(all(isfinite(b_data.data(:))));

            [~, idx] = max(abs(b_data.data));
            peak_x = scan.x(idx);
            peak_z = scan.z(idx);
            testCase.verifyEqual(peak_x, 0, 'AbsTol', 2e-3);
            testCase.verifyEqual(peak_z, 20e-3, 'AbsTol', 2e-3);
        end
    end

end
