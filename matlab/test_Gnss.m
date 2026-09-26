% Regression tests for the Gnss class.

classdef test_Gnss < matlab.unittest.TestCase
    methods (Test)
        function test_explicitStartTimeNotIgnored(testCase)
            % B10: an explicit gpsStartTime_s must not be silently
            % clobbered by the max-TOC auto-detect, even when it equals
            % the old value.
            rng(1);
            explicitT = g_time(datetime(2019,03,02,13,00,00));
            G         = Gnss('gpsStartTime_s', explicitT);
            testCase.verifyTrue(G.gpsStartTime_s == explicitT);
        end

        function test_defaultStartTimeAutoDetected(testCase)
            % B10 / regression (a): default (no explicit start) still
            % auto-detects from max TOC to the known-good value.
            rng(1);
            G        = Gnss();
            expected = g_time(datetime(2019,03,02,14,00,00));
            testCase.verifyTrue(G.gpsStartTime_s == expected);
            testCase.verifyEqual(G.gpsStartTime_s.g_sow, 568800, 'AbsTol', 1e-6);
        end

        function test_fieldSetConsistentAcrossVisibility(testCase)
            % B11: rawxSoln's field set must not change with visibility.
            rng(1);
            G     = Gnss();
            [posEcef,...
             velEcef,...
             Rpy]  = test_Gnss.samplePose();
            tStart = G.gpsStartTime_s.g_sow;

            % default mask angle
            solnVisible  = G.update(tStart, posEcef, velEcef, Rpy);

            % mask everything
            G.mask_angle = 90;
            solnNone     = G.update(tStart+1, posEcef, velEcef, Rpy);

            testCase.verifyEqual(sort(fieldnames(solnVisible)),     sort(fieldnames(solnNone)));
            testCase.verifyEqual(sort(fieldnames(solnVisible.svg)), sort(fieldnames(solnNone.svg)));
            testCase.verifyNotEmpty(solnNone.SOW);
            testCase.verifyNotEmpty(solnNone.bias_m);
            testCase.verifyNotEmpty(solnNone.pos_ecef);
        end

        function test_satelliteSettingGoesBackToEmpty(testCase)
            % Regression (b): the obs-staleness bug - a satellite that
            % sets must not keep reporting last epoch's values.
            rng(1);
            G      = Gnss();
            [posEcef,...
             velEcef,...
             Rpy]  = test_Gnss.samplePose();
            tStart = G.gpsStartTime_s.g_sow;

            soln1       = G.update(tStart, posEcef, velEcef, Rpy);
            visiblePrns = find(arrayfun(@(s) ~isempty(s.C1C), soln1.svg));
            testCase.assumeNotEmpty(visiblePrns);
            prn         = visiblePrns(1);

            G.mask_angle = 90;
            soln2        = G.update(tStart+1, posEcef, velEcef, Rpy);

            testCase.verifyEmpty(soln2.svg(prn).C1C);
        end

        function test_incIonoFlagActuallyGatesIono(testCase)
            % Regression (c): guard against the config-wiring bug - the
            % inc* flags must actually switch models off, not be
            % silently ignored.
            rng(1);
            Gon  = Gnss('incIono', true);
            rng(1);
            Goff = Gnss('incIono', false);

            [posEcef, velEcef, Rpy] = test_Gnss.samplePose();
            tStart = Gon.gpsStartTime_s.g_sow;

            rng(2);
            solnOn  = Gon.update(tStart, posEcef, velEcef, Rpy);
            rng(2);
            solnOff = Goff.update(tStart, posEcef, velEcef, Rpy);

            visible = find(arrayfun(@(s) ~isempty(s.C1C), solnOn.svg));
            testCase.assumeNotEmpty(visible);
            prn = visible(1);

            testCase.verifyNotEqual(solnOn.svg(prn).C1C, solnOff.svg(prn).C1C);
        end
    end

    methods (Static, Access = private)
        function [posEcef, velEcef, Rpy] = samplePose()
            lat     = deg2rad(52.9519816);
            lon     = deg2rad(-1.1907585);
            hd      = 100;
            [xEcef,...
             yEcef,...
             zEcef] = geodetic2ecef(lat,lon,hd,'wgs84');
            posEcef = [xEcef;yEcef;zEcef];
            velEcef = [0;0;0];
            Rpy     = [0;0;0];
        end
    end
end
