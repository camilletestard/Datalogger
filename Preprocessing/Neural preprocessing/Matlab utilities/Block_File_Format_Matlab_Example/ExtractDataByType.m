function [ dataArray ] = ExtractDataByType( rawData, HeaderStruct, dataTypeIndex, blockStartIndices, numberOfBlocks)
%EXTRACTDATABYTYPE extracts data of a single type from blocks and concatenates them into
%arrays
    dataPartitionSize = HeaderStruct.PartitionInfo(dataTypeIndex).PartitionSize;
    dataPartitionStart = HeaderStruct.PartitionInfo(dataTypeIndex).DataStart;
    dataArray = uint8(zeros(numberOfBlocks * dataPartitionSize, 1)); %preallocate       
    for blockIndex = 1:numberOfBlocks % for each block
        startIndexInSourceArray = blockStartIndices(blockIndex) + dataPartitionStart;
        startIndexInDestArray = (blockIndex - 1) * dataPartitionSize + 1;
        dataArray(startIndexInDestArray:startIndexInDestArray + dataPartitionSize - 1) = rawData(startIndexInSourceArray:startIndexInSourceArray + dataPartitionSize - 1);
    end

end

