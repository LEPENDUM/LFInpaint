% AddLibPath
% Adds an include directory, a lib directory and a list of lib filenames to
% the list.

% Inputs:
%   - include_dir : include directory to add to the list.
%   - lib_dir : library directories to add to the list.
%   - libList : cell array of string of the library file names to add to the list.
%   - extIDirs0 : (Optional) input list of include directories (starts empty if not given).
%   - extLDirs0 : (Optional) input list of library directories (starts empty if not given).
%   - libs0 :  (Optional) input list of lib filenames (starts empty if not given).
% Note : the lists should be either cells of strings or a single string (or empty).
%
% Outputs:
%   -extIDirs,extLDirs,libs : new lists after adding the input data.

function [extIDirs,extLDirs,libs] = AddLibPath(include_dir,lib_dir,libList,extIDirs0,extLDirs0,libs0)

if(~exist('extIDirs0','var')), extIDirs0='';end
if(~exist('extLDirs0','var')), extLDirs0='';end
if(~exist('libs0','var')), libs0='';end

%Add include folder:
if ~isempty(extIDirs0) && ischar(extIDirs0)
    extIDirs0 = {extIDirs0};
elseif ~isempty(extIDirs0) && ~iscellstr(extIDirs0)
    error('The argument ''extIDirs0'' should be either a string, a cell array of strings or empty!')
end
extIDirs = [extIDirs0, include_dir];

%Add lib folder:
if ~isempty(extLDirs0) && ischar(extLDirs0)
    extLDirs0 = {extLDirs0};
elseif ~isempty(extLDirs0) && ~iscellstr(extLDirs0)
    error('The argument ''extLDirs0'' should be either a string, a cell array of strings or empty!')
end
extLDirs = [extLDirs0, lib_dir];

%Add list of lib files:
if ~isempty(libs0) && ischar(libs0)
    libs0 = {libs0};
elseif ~isempty(libs0) && ~iscellstr(libs0)
    error('The argument ''libs0'' should be either a string, a cell array of strings or empty!')
end
libs = [libs0, libList];