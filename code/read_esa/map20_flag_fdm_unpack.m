function [flags] = map20_flag_fdm_unpack(packed)

% function to unpack 20 bit 18hz 'mapped' flags (uint32) to 20 seperate flags
%  [flags] = map20_flag_unpack(map_flag)

if ~isa(packed,'uint32')
    error('Input must be uint32(N,1) or uint32(1,N) array');
end

n_recs=length(packed);

k = [1:20 32];
flags = cast(bitget(repmat(packed,1,21),repmat(k,n_recs,1)),'uint8');

end