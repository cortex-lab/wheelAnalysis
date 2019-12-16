classdef test_wheel < matlab.unittest.TestCase
  % TEST_WHEEL Unit tests for findWheelMoves3
  %  Tests the results of wheel.findWheelMoves3 algorithm against previous
  %  results.  Also checks memory usage as this function used to be memory
  %  intensive.
  %
  %  How to add a test to the test data structure
  %   load('test_wheel.mat', 'test')
  %
  %   m=length(test)+1;
  %   data = load(block_files{103}); data = data.block;
  %   wval = data.inputs.wheelValues;
  %   wt = data.inputs.wheelTimes-data.events.expStartTimes;
  %
  %   test(m).wval = wval;
  %   test(m).wt = wt;
  %   [mon,mof] = wheel.findWheelMoves3(test(m).wval, test(m).wt, 1000);
  %   test(m).mon = mon;
  %   test(m).mof = mof;
  %
  %   save('test_wheel.mat', 'test')
  
  properties
    % A structure containing the test data and expected wheelmoves output
    tdata
    % The encoder resolution used for the test dataset
    res = 1024
    % The wheel radius in cm used for the test dataset
    r = 3.1
  end
  
  
  methods(TestClassSetup)
    function read_data(testCase)
      % Load the test data.  This should be a MAT file with the same name
      % and location as this class
      testCase.tdata = load(mfilename('fullpath'));
      testCase.tdata = testCase.tdata.test;
    end
  end
  
  methods(Test)
    function test_values(testCase)
      % Loads test data and runs algorithm, then checks the onset and
      % offset values match as well as the movement amplitudes and peak
      % velocities.  The data are small so this test is quick.
      data = testCase.tdata(1);
      [mon, mof, amp, vel] = wheel.findWheelMoves3(data.wval, data.wt, 1000);
      testCase.verifyEqual(mon, data.mon, 'Unexpected movement onset times')
      testCase.verifyEqual(mof, data.mof, 'Unexpected movement offset times')
      testCase.verifyEqual(amp, data.amp, 'Unexpected movement amplitudes')
      testCase.verifyEqual(vel, data.vel, 'Unexpected peak velocity times')
    end
    
    function test_memory_use(testCase)
      % This test uses a long-ish session of 40 min.  With the previous
      % code this usually fails due to lack of memory.  This test is
      % even stricter, requiring a very little peak memory value to
      % pass.
      data = testCase.tdata(2);
      % Ensure profiler turned off on error
      testCase.addTeardown(@profile, 'off')
      profile('-memory', 'on') % Record memory usage
      [mon, mof, amp, vel] = wheel.findWheelMoves3(data.wval, data.wt, 1000);
      info = profile('info'); % Stop profiler and gather info
      I = strcmp({info.FunctionTable.FunctionName}, 'findWheelMoves3');
      maxPeak = 25e6; % Expected peak memory should be < 25MB
      testCase.verifyLessThan(info.FunctionTable(I).PeakMem, maxPeak);
      testCase.verifyEqual(mon, data.mon, 'Unexpected movement onset times')
      testCase.verifyEqual(mof, data.mof, 'Unexpected movement offset times')
      testCase.verifyEqual(amp, data.amp, 'Unexpected movement amplitudes')
      testCase.verifyEqual(vel, data.vel, 'Unexpected peak velocity times')
    end
    
    function test_empty_dataset(testCase)
      % Tests the algorithm on a nonsense dataset
      data = testCase.tdata(3);
      [mon, mof, amp, vel] = wheel.findWheelMoves3(data.wval, data.wt, 1000);
      testCase.verifyEqual(mon, data.mon, 'Unexpected movement onset times')
      testCase.verifyEqual(mof, data.mof, 'Unexpected movement offset times')
      testCase.verifyEqual(amp, data.amp, 'Unexpected movement amplitudes')
      testCase.verifyEqual(vel, data.vel, 'Unexpected peak velocity times')
    end
    
    function test_displacement_units(testCase)
      % Tests the movement detection on data that has been extracted in
      % units of cm linear displacement
      data = testCase.tdata(4); % Data extracted from IBL FPGA setup
      cm = @(v) v./(4*testCase.res)*2*pi*testCase.r; % convert to cm
      params.posThresh = cm(8);
      params.posThreshOnset = cm(1.5);
      [mon, mof, amp, vel] = ...
        wheel.findWheelMoves3(data.wval, data.wt, 1000, params);
      testCase.verifyEqual(mon, data.mon, 'Unexpected movement onset times')
      testCase.verifyEqual(mof, data.mof, 'Unexpected movement offset times')
      testCase.verifyEqual(amp, data.amp, 'Unexpected movement amplitudes')
      testCase.verifyEqual(vel, data.vel, 'Unexpected peak velocity times')
    end
  end
end