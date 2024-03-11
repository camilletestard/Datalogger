classdef PartitionConstants
    %PARTITIONCONSTANTS Structure of the partition information
    %   
    
    properties (Constant = true)
        NumberOfDataTypes = 8;
        StructSize = 12;
        SizeInBytes = 4;
        DataTypePosition = 0;
        PartitionStartPosition = 4;
        PartitionSizePosition = 8;        
    end
    

end

