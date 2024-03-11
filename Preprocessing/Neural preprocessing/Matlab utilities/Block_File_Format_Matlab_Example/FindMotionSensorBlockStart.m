function [ startBlockIndices ] = FindMotionSensorBlockStart( audioData, constId )
%FINDMOTIONSENSORSTART Find start of each motion sensor block
startBlockIndices = zeros(10000, 1);
count = 0;
    for dataIdx = 1:length(audioData) - 1
        if (audioData(dataIdx) == constId(1) && audioData(dataIdx + 1) == constId(2))
            count = count + 1;
            startBlockIndices(count) = dataIdx;
        end
    end
    startBlockIndices = startBlockIndices(1:count);
            


end

