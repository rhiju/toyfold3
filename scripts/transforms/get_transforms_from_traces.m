function pts = get_transforms_from_traces(traces,i,j);
% pts = get_transforms_from_traces(traces,i,j);
%
% (C) R. Das, Stanford University

NITER = length( traces );
for n = 1:NITER
    trace = traces{n};
    [t(:,n),R(:,:,n)] = get_transform(...
        trace.t(:,i),trace.R(:,:,i),...
        trace.t(:,j),trace.R(:,:,j));
end
TransformSet = struct( 't',t,'R',R);
pts = fill_T6_from_t_and_R( TransformSet );
