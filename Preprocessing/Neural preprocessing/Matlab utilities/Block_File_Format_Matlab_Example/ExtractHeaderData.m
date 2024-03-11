function [ HeaderStruct ] = ExtractHeaderData( header )
%EXTRACTHEADERDATA extracts meta data from header and returns a header
%struct containing meta-data
      HeaderStruct.FileFormatId = GetValueFromHeader(header, HeaderConstants.FileFormatIdPosition, HeaderConstants.FileFormatIdBytes, 'uint32');
      HeaderStruct.SizeOfBlock = GetValueFromHeader(header, HeaderConstants.BlockSizePosition, HeaderConstants.BlockSizeBytes, 'uint32');
      HeaderStruct.Timestamp = GetValueFromHeader(header, HeaderConstants.TimestampPosition, HeaderConstants.TimestampBytes, 'uint32'); %ms since midnight
      PartitionStruct = struct('DataType', cell(1, PartitionConstants.NumberOfDataTypes), ...
      'DataStart', cell(1, PartitionConstants.NumberOfDataTypes), 'PartitionSize', cell(1,PartitionConstants.NumberOfDataTypes));

    for dataTypeIdx = 1:HeaderConstants.MaxNumberOfDataTypes
        dataStructOffset = HeaderConstants.PartitionStartPosition + PartitionConstants.StructSize * (dataTypeIdx - 1);
        thisDataType = GetValueFromHeader(header, dataStructOffset, PartitionConstants.SizeInBytes, 'uint32');
        switch thisDataType
            case uint32(DataTypeEnum.NoData)
                continue;
            case uint32(DataTypeEnum.EventData)
                PartitionStruct(dataTypeIdx).DataType = DataTypeEnum.EventData;
            case uint32(DataTypeEnum.NeuralData)
                PartitionStruct(dataTypeIdx).DataType = DataTypeEnum.NeuralData;
            case uint32(DataTypeEnum.MotionSensor)
                PartitionStruct(dataTypeIdx).DataType = DataTypeEnum.MotionSensor;
            case uint32(DataTypeEnum.Audio)
                PartitionStruct(dataTypeIdx).DataType = DataTypeEnum.Audio;
            case uint32(DataTypeEnum.RTK)
                PartitionStruct(dataTypeIdx).DataType = DataTypeEnum.RTK;
        end               
        partitionStartPosition = dataStructOffset + PartitionConstants.PartitionStartPosition;
        PartitionStruct(dataTypeIdx).DataStart = GetValueFromHeader(header, partitionStartPosition, PartitionConstants.SizeInBytes, 'uint32'); 
        partitionSizePosition = partitionStartPosition + PartitionConstants.PartitionStartPosition;
        PartitionStruct(dataTypeIdx).PartitionSize = GetValueFromHeader(header, partitionSizePosition, PartitionConstants.SizeInBytes, 'uint32'); 
    end
    HeaderStruct.PartitionInfo = PartitionStruct;
end
    


function [headerValue] = GetValueFromHeader(header, startPosition, numberOfBytes, dataTypeString)
    headerValueTemp = header(startPosition:startPosition + numberOfBytes - 1);
    headerValue = typecast(headerValueTemp, dataTypeString);
end


