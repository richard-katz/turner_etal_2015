function A = loadMorOutput(filename)
%
% USAGE: A = loadMorOutput(<filename>)
%
%   input: where <filename> is be the name of the PETSc
%      binary file to be loaded.  <filename>.info
%      MUST be present in the directory also.
%
%   output: a structure containing parameters, fields
%      and coordinates extracted from the input file.
%
   
   % if no filename is given, load the continuation file
   if (nargin==0 || isempty(filename)); 
       filename = ls('*_cont'); 
       filename = filename(1:end-1); 
   end
   
   try
       % for output generated with a PetscViewerBinaryOpen
       A = PetscReadBinaryMatlab([filename,'.info']);
   catch
      errstr = sprintf('\nloadMorOutput: cannot load file %s',filename);
      error(errstr); 
   end
   
   A = ProcessPetscOutput(A);
   
   return;

function A = ProcessPetscOutput(A)
    A = MakeCoordinateVectors(A);
    A.par.sec_per_year = 365.25*24*60^2;
   
function A = MakeCoordinateVectors(A)
    A.par.is = 3; A.par.ie = A.par.ni-2;
    A.par.js = 3; A.par.je = A.par.nj-2;
    A.par.hdim  = A.par.height/(A.par.nj-4); % grid spacing in kilometres
    A.coord.zc = ([0:A.par.nj-1]-1.5)*A.par.hdim;
    A.coord.xc = ([0:A.par.ni-1]-1.5)*A.par.hdim - A.par.xridge;
    A.coord.iin = [3:A.par.ni-2];
    A.coord.xe = A.coord.xc + A.par.hdim/2;
    A.coord.ze = A.coord.zc + A.par.hdim/2;
    A.coord.jin = [3:A.par.nj-2];   