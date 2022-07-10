function ntp_slope_codegen(nk, ni, nj, Zvec, OPTS)
%NTP_SLOPE_ERROR_CODEGEN Create MEX function for ntp_slope_error
%
%
% ntp_slope_error_codegen(nk, ni, nj, false)
% runs codegen on ntp_slope_error.m, appropriate for a grid
% of ni by nj points in the horizontal and nk points in the vertical.
%
% ntp_slope_error_codegen(nk, ni, nj, true)
% specifies that Z in ntp_slope_error.m is just a vector: Z(k)
% specifies the pressure or depth of all grid points having vertical index
% k. Use this for simple Z-level models (not hybrid coordinate models).
%
% ntp_slope_error_codegen(..., OPTS)
% overrides default verbosity by OPTS.VERBOSE and file output by
% OPTS.FILE_ID.

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% Set defaults
VERBOSE = 1; % verbose mode
FILE_ID = 1; % output to MATLAB terminal

% Override defaults
if nargin == 2 && isstruct(OPTS)
    if isfield(OPTS, 'VERBOSE')
        VERBOSE = OPTS.VERBOSE;
    end
    if isfield(OPTS, 'FILE_ID')
        FILE_ID = OPTS.FILE_ID;
    end
end

V = filesep();
folder_start = pwd();

try
    
    % Get info about which functions are on the path
    name = 'ntp_slope';
    name_mex = [name '_mex'];
    file_mex = dir(which(name_mex));
    file_mat = dir(which(name));
    assert(~isempty(file_mat), ['Cannot locate ' name '.m']);
    which_eos   = which('eos');
    file_eos    = dir(which_eos);
    
    % Test values
    s = 34.5;
    t = 3;
    p = 1000;
    m = eos(s, t, p);
    
    % Create textual identifier for this build of the MEX function.
    build_text = sprintf('%s_k%d_m=%.59e', name, nk, m);
    fileName_build_text = [file_mat.folder V name '_info.txt'];
    fileID_build_text = fopen(fileName_build_text, 'rt');
    
    
    if isempty(file_mex) ... % No mex file yet
            || file_mex.datenum < max([file_mat.datenum, file_eos.datenum]) ... % MEX is too old
            || (fileID_build_text < 0) ... % No text file specification
            || (fileID_build_text >= 0 && ~strcmp(fgetl(fileID_build_text), build_text) && ~fclose(fileID_build_text)) % Requested parameters do not match those recorded in text file. Also close the text file.
        
        if VERBOSE
            mytic = tic;
            fprintf(FILE_ID, 'Compiling MEX for %s, with\n', name);
            fprintf(FILE_ID, ' %s in %s\n', read_function_name(which_eos), which_eos);
            fprintf(FILE_ID, ' eos(%g,%g,%g) = %e\n', s, t, p, m);
        end

        vs = true;
        t_SppZ = coder.typeof(0, [8, nk-1, ni, nj], [true, vs, vs, vs]);
%         if Zvec
%             t_Z  = coder.typeof(0, [nk, 1], [vs, false]);
%         else
        t_Z  = coder.typeof(0, [nk, ni, nj], [vs, vs, vs]);
%         end
        t_z = coder.typeof(0, [ni, nj], [vs, vs]);
        
        % (SppZ, TppZ, Z, z, tolz, dx, dy)
        args = {t_SppZ, t_SppZ, t_Z, t_z, 1e-6, t_z, t_z};
        
        % Configure MEX for speed.
        mexconfig = coder.config('mex');
        mexconfig.ExtrinsicCalls = false;
        mexconfig.ResponsivenessChecks = false;
        mexconfig.IntegrityChecks = false;
        
        % Compile the MEX function
        cd(file_mat.folder);
        codegen(name, '-args', args, '-config', mexconfig, '-o', [file_mat.folder V name_mex], '-d', [file_mat.folder V 'codegen' V 'mex' V name]);
        clear(name_mex)  % Ensure new mex file gets used
        
        % Save textual identifier for this MEX function
        fileID_build_text = fopen(fileName_build_text, 'wt'); % overwrite
        fprintf(fileID_build_text, build_text);
        fclose(fileID_build_text);
        
        if VERBOSE
            fprintf(FILE_ID, ' ...done compiling, time %.2f\n', toc(mytic));
        end
    end
    
    cd(folder_start);
    
catch err
    cd(folder_start);
    rethrow(err);
end




