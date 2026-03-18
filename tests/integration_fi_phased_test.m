classdef integration_fi_phased_test < matlab.unittest.TestCase

    methods (Test)
        function test_fi_phased_with_sector_scan(testCase)
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

            N = 50;
            azimuth_axis = linspace(-10*pi/180, 10*pi/180, N).';
            depth = 40e-3;
            seq = uff.wave();
            for n = 1:N
                seq(n) = uff.wave();
                seq(n).probe = prb;
                seq(n).source = uff.point();
                seq(n).source.azimuth = azimuth_axis(n);
                seq(n).source.distance = depth;
                seq(n).apodization = uff.apodization();
                seq(n).apodization.window = uff.window.tukey50;
                seq(n).apodization.f_number = 1.7;
                seq(n).apodization.focus = uff.sector_scan('xyz', seq(n).source.xyz);
                seq(n).sound_speed = pha.sound_speed;
            end

            sim = fresnel();
            sim.phantom = pha;
            sim.pulse = pul;
            sim.probe = prb;
            sim.sequence = seq;
            sim.sampling_frequency = 41.6e6;
            channel_data = sim.go();

            depth_axis = linspace(35e-3, 45e-3, 100).';
            scan = uff.sector_scan('azimuth_axis', azimuth_axis, ...
                                   'depth_axis', depth_axis);

            mid = midprocess.das();
            mid.code = code.matlab;
            mid.dimension = dimension.both;
            mid.channel_data = channel_data;
            mid.scan = scan;

            mid.transmit_apodization.window = uff.window.scanline;

            mid.receive_apodization.window = uff.window.tukey50;
            mid.receive_apodization.f_number = 1.7;

            b_data = mid.go();

            testCase.verifyEqual(numel(b_data.data(:)), N * 100);
            testCase.verifyTrue(all(isfinite(b_data.data(:))));

            [~, idx] = max(abs(b_data.data));
            peak_x = scan.x(idx);
            peak_z = scan.z(idx);
            testCase.verifyEqual(peak_x, 0, 'AbsTol', 3e-3);
            testCase.verifyEqual(peak_z, 40e-3, 'AbsTol', 3e-3);
        end
    end

end
