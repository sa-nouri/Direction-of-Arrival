%SETUP  Compile the SPGL1 MEX interfaces

%   setup.m
%   $Id: spgsetup.m

root = pwd;
try
    cd('private')
    mex oneProjectorMex.c oneProjectorCore.c heap.c -output oneProjectorMex -DNDEBUG
    fprintf('Successfully compiled oneProjector.\n');
    cd(root)
catch
    cd(root)
    fprintf('Could not compile oneProjector.');
    fprintf('You can still use the slower ".m" version.');
    rethrow(lasterr);
end
    
