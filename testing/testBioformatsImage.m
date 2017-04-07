classdef testBioformatsImage < matlab.unittest.TestCase
    
    properties
        
        testfile = '3-29-17_(640nm-8percent)_Cy5-20-1frame_GFP-8-100ms.nd2';
        
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
                to.filename = 'test.nd2';
            catch ME
                                
            end
            
            TestCase.verifyEqual(ME.identifier,'BioformatsImage:CannotFindFile');
            
        end
        
        function testSetFilename(TestCase)
            
            to = BioformatsImage;
            
            to.filename = TestCase.testfile;
            
            TestCase.verifyEqual(to.filename, which(TestCase.testfile))
            
        end
    end
end