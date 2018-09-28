classdef test_wheel < matlab.unittest.TestCase
    
    properties
        tdata
    end
    
    
    methods(TestMethodSetup)
        function read_data(testCase)
            [chem, nam] = fileparts(which('test_wheel'));
            testCase.tdata  = load([chem filesep nam '.mat' ]);
            testCase.tdata   = testCase.tdata.test;
        end
    end
    
    methods(Test)
        function test_values(testCase)
            tdata = testCase.tdata(1);
            [mon,mof] = wheel.findWheelMoves3(tdata.wval, tdata.wt, 1000, []);
            testCase.assertEqual(mon, tdata.mon)
            testCase.assertEqual(mof, tdata.mof)
        end
        % not really testing on a big computer...
        function test_out_of_memory(testCase)
            tdata = testCase.tdata(2);
            [mon,mof] = wheel.findWheelMoves3(tdata.wval, tdata.wt, 1000, []);
            testCase.assertEqual(mon, tdata.mon)
            testCase.assertEqual(mof, tdata.mof)
        end
        function test_empty_dataset(testCase)
            tdata = testCase.tdata(2);
            [mon,mof] = wheel.findWheelMoves3(tdata.wval, tdata.wt, 1000, []);
            testCase.assertEqual(mon, tdata.mon)
            testCase.assertEqual(mof, tdata.mof)
        end
    end
end

%% How to add a test to the test data structure
% load('/home/owinter/MATLAB/Rigbox/wheelAnalysis/test_wheel.mat', 'test')
% 
% m=length(test)+1;
% data = load(block_files{103}); data = data.block;
% wval = data.inputs.wheelValues;
% wt = data.inputs.wheelTimes-data.events.expStartTimes;
% 
% test(m).wval = wval;
% test(m).wt = wt;
% [mon,mof] = wheel.findWheelMoves3(test(m).wval, test(m).wt, 1000, []);
% test(m).mon = mon;
% test(m).mof = mof;
% 
% %
% save('/home/owinter/MATLAB/Rigbox/wheelAnalysis/test_wheel.mat', 'test')