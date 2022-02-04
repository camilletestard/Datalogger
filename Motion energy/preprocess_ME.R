# Script to pull together ME data

library(xlsx)

#Load data
file = file.choose() # choose the formatted behavior file
log <- xlsx::read.xlsx(file, 1)

LeftVideo_file = file.choose() # choose the left video ME file
ME_left_file = read.delim(LeftVideo_file, header = TRUE, sep = " ")
ME_left = ME_left_file[,1]

RightVideo_file = file.choose() # choose the right video ME file
ME_right_file = read.delim(RightVideo_file, header = TRUE, sep = " ")
ME_right = ME_right_file[,1]

#First important note: a few frames difference between the two videos.
diff_vids = nrow(ME_right_file) - nrow(ME_left_file)
print(paste('Difference in frames between videos is:', abs(diff_vids), 'frames.'))
if (diff_vids <0){
  ME_right[(length(ME_right)+1):(length(ME_right)+abs(diff_vids))]=0
} else {
  ME_left[(length(ME_left)+1):(length(ME_left)+abs(diff_vids))]=0
}

MotionEnergy = data.frame(cbind(ME_left, ME_right))

# #Visualize ME
# plot(1:3000, ME_right[1:3000], type = "S", frame = FALSE, 
#      col = "red", xlab = "x", ylab = "y")
# lines(1:3000, ME_left[1:3000], col = "blue", type = "S", lty = 2)

#Get recording length in sec
length_recording = log[nrow(log),"Time"]/1000
#Important note: the length of the recording session should not be more than 2sec different between the
#deuteron log and "unit rasters". That is because in the script "log_read_sorted_NEX_file_CombinedChannels.m"
#Unit_rasters is initialized by the line: Unit_rasters(unit,:) = zeros(1,round(length_recording*temp_resolution)+1);
#round and +1 can lead to up to 2sec difference with the length of recording from the log.
time_before_cameraSync = log[2,"Time"]/1000
time_after_cameraSync = log[nrow(log),"Time"]/1000 - log[nrow(log)-1,"Time"]/1000

#Estimate frame rate
length_video = (log[n row(log)-1,"Time"] - log[2,"Time"])/1000
frame_rate = 29.971
expected_num_frames = frame_rate*length_video
diff_frames = abs(expected_num_frames-length(ME_right))

#Pad ME with 0 in the extra recording spots
total_frames_needed = round(log[nrow(log),"Time"]/1000*frame_rate)
numFrames_before_cameraSync = round(time_before_cameraSync*frame_rate)
numFrames_after_cameraSync = round(time_after_cameraSync*frame_rate)
zero_pad_before = data.frame(matrix(data = 0, nrow = numFrames_before_cameraSync, ncol = 2)); names(zero_pad_before) = c("ME_left","ME_right")
zero_pad_after = data.frame(matrix(data = 0, nrow = numFrames_after_cameraSync, ncol = 2)); names(zero_pad_after) = c("ME_left","ME_right")

ME_final = rbind(zero_pad_before, MotionEnergy, zero_pad_after)

#Save ME file
mac =0
if (mac == 1) { home = '~'} else
{home = 'C:/Users/GENERAL'}

setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/', sep="" ))
write.csv(ME_final, 'Amos_2021-07-29/ME.csv',row.names=F)
