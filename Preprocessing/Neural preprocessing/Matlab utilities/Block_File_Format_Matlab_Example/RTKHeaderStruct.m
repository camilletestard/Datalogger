classdef RTKHeaderStruct
    %RTKPARTITIONCONSTANTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant = true)
        BlockHeaderConstantStart = 1;
        BlockHeaderConstantLength = 4;
        NumDataBytesStartIndex = 5;%bytes 5 - 8 inclusive are a uint32 stating the number of data bytes in this block
        NumDataBytesSizeInBytes = 4;
        DataStartIndex = 17;
        
    end
    

    
end

