classdef integration_coherent_compounding_test < matlab.unittest.TestCase

    methods (Test)
        function test_coherent_compounding_reduces_waves(testCase)
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

            Nx = 64; Nz = 64;
            scan = uff.linear_scan('x_axis', linspace(-5e-3, 5e-3, Nx).', ...
                                   'z_axis', linspace(15e-3, 25e-3, Nz).');

            pipe = pipeline();
            pipe.channel_data = channel_data;
            pipe.scan = scan;

            pipe.receive_apodization.window = uff.window.hanning;
            pipe.receive_apodization.f_number = 1.7;

            pipe.transmit_apodization.window = uff.window.hanning;
            pipe.transmit_apodization.f_number = 1.7;

            das = midprocess.das();
            das.code = code.matlab;

            b_data_cc = pipe.go({das postprocess.coherent_compounding()});

            testCase.verifyEqual(size(b_data_cc.data, 2), 1, ...
                'After compounding, wave dimension should be 1');
            testCase.verifyEqual(size(b_data_cc.data, 1), Nx * Nz);
            testCase.verifyTrue(all(isfinite(b_data_cc.data(:))));

            [~, idx] = max(abs(b_data_cc.data));
            peak_x = scan.x(idx);
            peak_z = scan.z(idx);
            testCase.verifyEqual(peak_x, 0, 'AbsTol', 2e-3);
            testCase.verifyEqual(peak_z, 20e-3, 'AbsTol', 2e-3);
        end
    end

end
