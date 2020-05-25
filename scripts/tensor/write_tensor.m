function write_tensor( T, filename, format )
% write_tensor( T, filename )
%
%  Create tensor with name like "my_tensor.bin.gz"
%  and create accompanying .json file with n_bins information.
%
% INPUTS
%   T = struct with .tensor and .json fields
%   filename = output filename (should be .bin, .txt, .bin.gz, .txt.gz)
%   format   = for txt output, the format (default is '%8.3f' )
%
% (C) Rhiju Das, Stanford University 2017
if ~exist( 'format', 'var' ) format = '%8.3f'; end;
assert( strcmp( T.json.type, class( T.tensor ) ) );
[filename, json_file, do_gzip, output_binary ] = process_tensor_filename( filename );


% the ordering of MathNTensor output in Rosetta requires reshaping and
% permuting to get into 'reasonable' order.
Ndim = length( size( T.tensor ) );
F = permute( T.tensor, [Ndim:-1:1] );

fid = fopen( filename, 'w');
if (output_binary )
    count = fwrite( fid, F(:), class(F));
else
    count = fprintf(fid, format, F(:));
    count = count / 8;
end
fprintf( 'Put %d values into: %s \n', count, filename );
fclose( fid );
assert( count == prod( size( F ) ) );

if do_gzip
    gzip( filename );
    delete( filename )
    fprintf( 'Gzipped into %s.gz\n', filename );
end

json = savejson( '', T.json, json_file );
fprintf( 'Created companion JSON file: %s\n', json_file );
