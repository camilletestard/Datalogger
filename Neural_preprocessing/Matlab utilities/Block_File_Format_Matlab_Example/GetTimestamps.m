function timestamps = GetTimestamps(firstTimestamp, frequency, dataLength)

%GETTIMESTAMPS gets the timestamp in ms from midnight of all data points
    firstTimestamp = double(firstTimestamp);
    stepSizeMs = 1 / frequency * 1000;
    lastTimestamp = firstTimestamp + (dataLength - 1) * stepSizeMs;
    timestamps = firstTimestamp:stepSizeMs:lastTimestamp;


end

