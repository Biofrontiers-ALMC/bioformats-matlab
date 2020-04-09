classdef testBioformatsImage < matlab.unittest.TestCase
    
    properties
        
        testfile = 'test.nd2';
%         zStackTestfile = '../local/zStackTest.nd2';
    end
    
    methods(TestClassSetup)
        function addTbxToPath(TestCase)
            %Adds the toolbox path to the path
            
            p = path;   %Store current path
            TestCase.addTeardown(@path,p);  %Restore it on teardown
            addpath('../tbx/BioformatsImage/');  %Add the toolbox path
        end
    end
    
    methods (Test)

        function setUnknownFilename(TestCase)
           
            to = BioformatsImage;
            
            TestCase.verifyClass(to,'BioformatsImage');
            
            try
                to.filename = 'noSuchFILENAME_anfi29nsi3.nd2';
            catch ME
                                
            end
            TestCase.verifyEqual(ME.identifier,'BioformatsImage:CannotFindFile');
            
            clear to
        end
        
        function setFilename(TestCase)
            
            to = BioformatsImage;
            
            to.filename = TestCase.testfile;
            
            TestCase.verifyEqual(to.filename, which(TestCase.testfile))
            
        end
        
        function getPlane(TestCase)
            
            BIobj = BioformatsImage(TestCase.testfile);
            
            imgBFI = BIobj.getPlane(1, 1, 1);
            
            nd2r = bfGetReader(TestCase.testfile);
            
            imgND2r = bfGetPlane(nd2r,nd2r.getIndex(0,0,0) + 1);
            
            TestCase.verifyEqual(imgBFI, imgND2r);
        end
        
        function getPlaneByChannelName (TestCase)
            BIobj = BioformatsImage(TestCase.testfile);
            
            imgBFI = BIobj.getPlane(1,'Mono',1);
            
            nd2r = bfGetReader(TestCase.testfile);
            
            imgND2r = bfGetPlane(nd2r,nd2r.getIndex(0,0,0) + 1);
            
            TestCase.verifyEqual(imgBFI, imgND2r);
            
        end
        
        function getTile(TestCase)
            
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
    
    methods (Test)
        
        function checkDelete (TestCase)
            %Related to Issue #259: Image remains "in use" in MATLAB
            %
            %After clearing the BioformatsImage object, the file remains in
            %memory (cannot rename or delete).
            %
            %Solution: Call the close() function on the bfReader object
            %before deleting the object.
            
            bfObj = BioformatsImage(TestCase.testfile);
            
            %Try renaming the file. This should cause an error since the
            %image file should be locked by MATLAB.
            [fpath, fname] = fileparts(TestCase.testfile);
            mvSuccess = movefile(fullfile(fpath,fname), fullfile(fpath,'new.nd2'));
            TestCase.assertEqual(mvSuccess,false);
            
            %Delete the obj
            clear bfObj
            
            mvSuccess = movefile(fullfile(fpath,fname), fullfile(fpath,'new.nd2'));
            TestCase.verifyEqual(mvSuccess,true);
            
            %Rename the file back to the original
            mvSuccess = movefile(fullfile(fpath,'new.nd2'),fullfile(fpath,fname));
        end
        
    end
    
end