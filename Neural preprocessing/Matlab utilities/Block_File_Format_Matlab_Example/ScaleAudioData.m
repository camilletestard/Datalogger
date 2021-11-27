function [ scaledAudioData ] = ScaleAudioData( isAudioSigned, audioDataAsBytes, numberOfAudioBits )
%SCALEAUDIODATA Summary of this function goes here
%   Detailed explanation goes here

if (isAudioSigned)
        audioData = single(typecast(audioDataAsBytes, 'int16'));
        minValue = -2^(numberOfAudioBits - 1) + 1;
        maxValue = 2^(numberOfAudioBits - 1);
        range = maxValue - minValue;
        for dataPointIndex = 1:length(audioData)
            scaledAudioData(dataPointIndex) = (audioData(dataPointIndex) - minValue) / range;
        end   
    else
        audioData = single(typecast(audioDataAsBytes, 'uint16'));
        scaledAudioData = zeros(length(audioData), 1);
        minValue = 0;
        maxValue = 2^(numberOfAudioBits);
        range = maxValue - minValue;
        for dataPointIndex = 1:length(audioData)
            scaledAudioData(dataPointIndex) = (audioData(dataPointIndex) - minValue) / range;
        end
end
end

