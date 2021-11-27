classdef DataTypeEnum < uint32
    %DATATYPEENUM the enum for types of data stored here
    
    enumeration
        NoData (0),
        EventData (1),
        NeuralData (2),
        MotionSensor (3),
        Audio (4),
        RTK (8)
    end
        
end

