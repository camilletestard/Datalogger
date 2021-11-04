function [ RTKData, RTKTimestamp ] = ExtractRTKData( RTKPartitionedData, headerStruct, rawData, blockStartIndices )
%EXTRACTRTKDATA extracts actual RTK data from the extracted partitioned
%data, which contains extraneous zeros as well as headers. Returns actual
%RTK data and timestamps of each event
    fun = @(x) headerStruct.PartitionInfo(x).DataType == DataTypeEnum.RTK;
    RTKPartitionSize = headerStruct.PartitionInfo(find(cell2mat(arrayfun(fun, 1:size(headerStruct.PartitionInfo, 2), 'un', 0)))).PartitionSize; % divide RTK data into blocks
    RTKBlockArray = reshape(RTKPartitionedData, RTKPartitionSize, []);
    RTKNumDataBytes =  RTKBlockArray(RTKHeaderStruct.NumDataBytesStartIndex:RTKHeaderStruct.NumDataBytesStartIndex + RTKHeaderStruct.NumDataBytesSizeInBytes - 1, :);
    RTKNumDataBytes = cell2mat(arrayfun(@(x) typecast(RTKNumDataBytes(:, x), 'uint32'), 1:size(RTKNumDataBytes, 2), 'un', 0)); % get number of bytes as uint32
    RTKData = cell(sum(RTKNumDataBytes ~= 0), 1);
    RTKTimestamp = cell(size(RTKData));
    count = 1;
    for blockIndex = 1:length(RTKNumDataBytes) % for each block, extract data bytes
        if RTKNumDataBytes(blockIndex) == 0 % no RTK data in this block
            continue;
        end
        startIndex = RTKHeaderStruct.DataStartIndex;
        endIndex = startIndex + RTKNumDataBytes(blockIndex) - 1;
        RTKData{count} = RTKBlockArray(startIndex:endIndex, blockIndex);
        startOfHeader = blockStartIndices(blockIndex);
        endOfHeader = startOfHeader + HeaderConstants.HeaderTotalBytes - 1;
        headerStruct = ExtractHeaderData(rawData(startOfHeader:endOfHeader));
        RTKTimestamp{count} = headerStruct.Timestamp;
        count = count + 1;

        
    end
    
    
end

