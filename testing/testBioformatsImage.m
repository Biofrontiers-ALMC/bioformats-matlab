classdef testBioformatsImage < matlab.unittest.TestCase
    
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
        
        function testSetFilename(TestCase)
            
            to = BioformatsImage;
            
            to.filename = TestCase.testfile;
            
            TestCase.verifyEqual(to.filename, which(TestCase.testfile))
            
        end
        
        function testGetPlane(TestCase)
            
            BIobj = BioformatsImage(TestCase.testfile);
            
            imgBFI = BIobj.getPlane([1 1 1]);
            
            nd2r = bfGetReader(TestCase.testfile);
            
            imgND2r = bfGetPlane(nd2r,nd2r.getIndex(0,0,0) + 1);
            
            TestCase.verifyEqual(imgBFI, imgND2r);
        end
        
        function testGetTile(TestCase)
            
            %[Z C T]
            iZ = 1;
            iC = 1;
            iT = 10;
            
            BIobj = BioformatsImage(TestCase.testfile);
            
            [tileImg, tileROI] = BIobj.getTile([iZ, iC, iT],[4 5], 5);
                        
            nd2r = bfGetReader(TestCase.testfile);
            
            imgND2r = bfGetPlane(nd2r,nd2r.getIndex(iZ - 1, iC - 1, iT - 1) + 1,...
                tileROI(1), tileROI(2), tileROI(3), tileROI(4));
            
            TestCase.verifyEqual(tileImg, imgND2r);
                    
        end
        
    end
end