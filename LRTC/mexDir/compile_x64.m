%Compile a mex file.
%
% Inputs :
%   - filename : name of the c or c++ source file to compile.
%   - extIDirs : (Optional) list of include directories containing include files of external dependencies.
%   - extLDirs0 : (Optional) list of library directories containing lib files of external dependencies.
%   - libs0 :  (Optional) list of lib filenames of external dependencies.
%   - debugMode : (Optional:default=false), compile in debug mode.
% Note : the lists should be either cells of strings or a single string (or empty).
% See the function AddLibPath to generate the lists.
%
% Dependencies (.h and .lib) can be added to the ./include and ./lib/x64 directories in the
% same folder as this script. Only the list of lib filenames is necessary then.
%
% The compiled mex file is generated in the ./x64 directory in the same folder as this script.

function compile_x64(filename,extIDirs,extLDirs,libs,debugMode)
clear mex;
if(~exist('extIDirs','var')), extIDirs='';end
if(~exist('extLDirs','var')), extLDirs='';end
if(~exist('libs','var')), libs='';end

debugOpt = '';
if(exist('debugMode','var') && debugMode==true)
    debugOpt = '-g';
end

%Initialize paths:
src_dir = [fileparts(mfilename('fullpath')) '/'];
output_dir =    ['"' src_dir 'x64/"'];
include_dir =   ['-I"' src_dir 'include/"'];
lib_dir =       ['-L"' src_dir 'libs/x64/"'];
lib_flag='';

%Generate source filename with absolute path from the input filename which
%can be given either with absolute or relative path.
if(exist([src_dir filename],'file'))
    src_file = ['"' src_dir filename '"'];
elseif(exist(filename,'file'))
    src_file = ['"' filename '"'];
else
    error(['Can''t find the source file ' filename]);
end

%Setup external library include folders:
if ~isempty(extIDirs) && ischar(extIDirs)     %ischar => single string (array of char).
    extIDirs = {extIDirs};
elseif ~isempty(extIDirs) && ~iscellstr(extIDirs)	%iscellstr() => cell array of strings.
    error('The argument ''extIDirs'' should be either a string, a cell array of strings or empty!')
end

for i=1:length(extIDirs)
    include_dir = [include_dir ' -I"' extIDirs{i} '"'];
end

%Setup external library lib folders:
if ~isempty(extLDirs) && ischar(extLDirs)       %ischar => single string (array of char).
    extLDirs = {extLDirs};
elseif ~isempty(extLDirs) && ~iscellstr(extLDirs)	%iscellstr() => cell array of strings.
    error('The argument ''extLDirs'' should be either a string, cell array of strings or empty!')
end

for i=1:length(extLDirs)
    lib_dir = [lib_dir ' -L"' extLDirs{i} '"'];
end

%Setup external library libs:
if ~isempty(libs) && ischar(libs)     %ischar => single string (array of char).
    libs = {libs};
elseif ~isempty(libs) && ~iscellstr(libs)	%iscellstr() => cell array of strings.
    error('The argument ''libs'' should be either a string, cell array of strings or empty!')
end

for i=1:length(libs)
    if(i>1) lib_flag = [lib_flag ' '];end
    lib_flag = [lib_flag '-l' libs{i}];
end


%Do compilation:
%That is the correct way to call the mex build function with variables as argument!!!
%(Why is this not mentioned anywhere in the mex documentation????)
mexCommand = [debugOpt ' -largeArrayDims -outdir ' output_dir ' ' include_dir  ' ' lib_dir ' ' lib_flag ' ' src_file];
args = textscan( mexCommand, '%s' );
fprintf('compiling with arguments:\n')
fprintf( '%s\n', args{1}{:} )
mex(args{1}{:});


%Those are incorrect ways that have worked in some cases but not always:
%Incorrect way #1:
%options=[debugOpt ' -largeArrayDims -outdir ' output_dir ' ' include_dir  ' ' lib_dir ' ' lib_flag];
%mex(['COMPFLAGS="$COMPFLAGS ' options '"'], [src_dir '\' filename]);
%
%Incorrect way #2:
%mex(debugOpt, '-largeArrayDims', include_dir, lib_dir, lib_flag, [src_dir '\' filename], '-outdir', output_dir);


clear mex;