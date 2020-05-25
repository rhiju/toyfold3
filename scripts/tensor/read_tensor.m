function T = read_tensor( filename )
% read_tensor( filename )
%
%  Read tensor from file with name like "my_tensor.bin" or
%  "my_tensor.bin.gz"
%
%  and look for accompanying .json file with n_bins information.
%
% (C) Rhiju Das, Stanford University

filename = strrep(filename,'//','/'); % causing an issue with gunzip

T.tensor = [];
[filename, json_file, is_gzip, use_binary ] = process_tensor_filename( filename );

% MATLAB claims to be able to read ascii from .gz, but this doesn't seem to work:
%if ( is_gzip & ~use_binary ) filename = [filename, '.gz' ]; end; 

if ( ~check_exists_or_gunzip( filename ) ) return; end;
if ( ~check_exists_or_gunzip( json_file ) ) return; end;

json = loadjson( json_file );
assert( isfield( json, 'n_bins' ) );
assert( isfield( json, 'type' ) );

fid = fopen( filename, 'r' );
if ( use_binary )
    [X, n_data] = fread( fid, json.type );
else
   [X, n_data] = fscanf(fid, '%f');
end
if ( prod( json.n_bins ) ~= n_data )
    fprintf( 'Mismatch between n_data %d and expected value based on product of n_bins %d\n', n_data, prod( json.n_bins) );
    return;
end
fclose( fid );
if is_gzip & exist( [filename,'.gz'], 'file' ) 
    fprintf( 'Deleting (since there is a .gz version): %s\n\n', filename );
    delete( filename );
end


% the ordering of MathNTensor output in Rosetta requires reshaping and
% permuting to get into 'reasonable' order.
Ndim = length( json.n_bins );
T.tensor = reshape( X, json.n_bins(Ndim:-1:1) );
T.tensor = permute( T.tensor, [Ndim:-1:1] );
T.json = json;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [success, did_gunzip] = check_exists_or_gunzip( filename )
success = 1;
if ~exist( filename, 'file' );
    filename_gz =  [filename, '.gz' ];
    if exist( filename_gz, 'file' )
        fprintf( 'Gunzipping: %s\n', filename_gz );
        gunzip( filename_gz );
    else
        fprintf( 'Could not find: %s or %s\n', filename, filename_gz );
        success = 0;
    end
end
