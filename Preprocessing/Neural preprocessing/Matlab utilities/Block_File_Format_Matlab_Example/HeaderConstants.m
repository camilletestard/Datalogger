classdef HeaderConstants
    % HEADERCONSTANTS contains constants relating to block header
    %   Contains position (1 indexed) and size in bytes within header of each value and any
    %   constant values relating to header
    
    properties (Constant = true)
        HexConstId = '1234ABCD567890EF';
        HeaderTotalBytes = 108;
        ConstIdPosition = 1;
        ConstIdBytes = 8;
        FileFormatIdPosition = 9;
        FileFormatIdBytes = 4;
        BlockSizePosition = 13;
        BlockSizeBytes = 4
        TimestampPosition = 17;
        TimestampBytes = 4;
        PartitionStartPosition = 25;
        PartitionBytes = 84;
        MaxNumberOfDataTypes = 7;
    end
   
end



