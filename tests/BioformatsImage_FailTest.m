classdef BioformatsImage_FailTest < matlab.unittest.TestCase
    
    properties
        
        testfile = 'test.nd2';
        
    end
    
    methods(TestClassSetup)
        function addTbxToPath(TestCase)
            %Adds the toolbox path to the path
            
            p = path;   %Store current path
            TestCase.addTeardown(@path,p);  %Restore it on teardown
            addpath('../tbx/bfimage/');  %Add the toolbox path
        end
    end
    
    methods (Test)

        function testSetUnknownFilename(TestCase)
           
            to = BioformatsImage;
            
            TestCase.verifyClass(to,'BioformatsImage');
            
            try
                to.filename = 'noSuchFILENAME_anfi29nsi3.nd2';
            catch ME
                                
            end
            TestCase.verifyEqual(ME.identifier,'BioformatsImage:CannotFindFile');
            
        end

    end
end