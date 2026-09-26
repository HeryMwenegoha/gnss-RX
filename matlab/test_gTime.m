% Regression test for gTime's power/mpower fix.

classdef test_gTime < matlab.unittest.TestCase
    methods (Test)
        function test_powerCombinesIntegerAndFractionalParts(testCase)
            % power/mpower used to return obj.t^b + obj.f^b instead of
            % (obj.t + obj.f)^b.
            d = gTime(3600,0) - gTime(0,0.07);
            expected = 12959496.0049;
            testCase.verifyEqual(d.^2, expected, 'AbsTol', 1e-6);
            testCase.verifyEqual(d^2,  expected, 'AbsTol', 1e-6);
        end
    end
end
