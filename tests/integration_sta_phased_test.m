classdef integration_sta_phased_test < matlab.unittest.TestCase

    methods (Test)
        function test_sta_phased_with_coherent_compounding(testCase)
            pha = uff.phantom();
            pha.sound_speed = 1540;
            pha.points = [0, 0, 40e-3, 1];

            prb = uff.linear_array();
            prb.N = 64;
            prb.pitch = 300e-6;
            prb.element_width = 270e-6;
            prb.element_height = 7000e-6;

            pul = uff.pulse();
            pul.center_frequency = 3e6;
            pul.fractional_bandwidth = 0.6;

            seq = uff.wave();
            for n = 1:prb.N_elements
                seq(n) = uff.wave();
                seq(n).probe = prb;
                seq(n).source.xyz = [prb.x(n) prb.y(n) prb.z(n)];
                seq(n).sound_speed = pha.sound_speed;
            end

            sim = fresnel();
            sim.phantom = pha;
            sim.pulse = pul;
            sim.probe = prb;
            sim.sequence = seq;
            sim.sampling_frequency = 4 * pul.center_frequency;
            channel_data = sim.go();

            scan = uff.sector_scan('azimuth_axis', linspace(-10*pi/180, 10*pi/180, 100).', ...
                                   'depth_axis', linspace(35e-3, 45e-3, 100).');

            pipe = pipeline();
            pipe.channel_data = channel_data;
            pipe.scan = scan;

            pipe.receive_apodization.window = uff.window.tukey50;
            pipe.receive_apodization.f_number = 1.7;

            pipe.transmit_apodization.window = uff.window.tukey50;
            pipe.transmit_apodization.f_number = 1.7;

            das = midprocess.das();
            das.code = code.matlab;

            b_data = pipe.go({das postprocess.coherent_compounding()});

            testCase.verifyEqual(size(b_data.data, 2), 1, ...
                'Coherent compounding should reduce wave dimension to 1');
            testCase.verifyTrue(all(isfinite(b_data.data(:))));

            [~, idx] = max(abs(b_data.data));
            peak_x = scan.x(idx);
            peak_z = scan.z(idx);
            testCase.verifyEqual(peak_x, 0, 'AbsTol', 3e-3);
            testCase.verifyEqual(peak_z, 40e-3, 'AbsTol', 3e-3);
        end
    end

end
