function [filename, json_file, is_gzip, is_binary ] = process_tensor_filename( filename );
%  [filename, json_file, is_gzip, is_binary ] = process_tensor_filename( filename );
%
% helper function for read_tensor & write_tensor.
%
is_gzip = 0;
if strcmp( filename(end-2:end), '.gz' )
    is_gzip = 1;
    filename = filename(1:end-3);
end
if (~strcmp(filename(end-3:end),'.bin') & ...
    ~strcmp(filename(end-3:end),'.txt') )
    fprintf( 'Filename should have .bin, .bin.gz, .txt, or txt.gz as suffix\n');
    return;
end
is_binary = strcmp( filename( end-3:end), '.bin' );
json_file = strrep( strrep( filename, '.bin', '.json' ),'.txt','.json' );
