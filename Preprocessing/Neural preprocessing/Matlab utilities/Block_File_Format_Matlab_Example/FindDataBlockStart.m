
function [ indicesOfSubarray ] = FindDataBlockStart( mainArray, subArray )
%FINDSUBARRAY finds all instances of subarray in mainArray

    dataAsDoubles = typecast(mainArray, 'uint64');
    constId = typecast(subArray, 'uint64');
    indicesOfSubarrayAsDouble = find(dataAsDoubles == constId);
    indicesOfSubarray = (indicesOfSubarrayAsDouble - 1) * length(subArray) + 1; %because MATLAB is 1 indexed, and switched to 8 byte uints

end

