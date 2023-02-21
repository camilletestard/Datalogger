# Script to pull together ME data

library(xlsx)

#Load data
file = file.choose() # choose the formatted behavior file
log <- xlsx::read.xlsx(file, 1)
date = substr(file,nchar(file)-24, nchar(file)-15)
if (length(grep('Hooke',file))>0){monkey ="Hooke"}else{monkey="Amos"}
session=paste(monkey,date,sep="_")

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
} 
if (diff_vids >0){
  ME_left[(length(ME_left)+1):(length(ME_left)+abs(diff_vids))]=0
}

MotionEnergy = data.frame(cbind(ME_left, ME_right))

# #Visualize ME
# plot(1:3000, ME_right[1:3000], type = "S", frame = FALSE, 
#      col = "red", xlab = "x", ylab = "y")
# lines(1:3000, ME_left[1:3000], col = "blue", type = "S", lty = 2)

#Get indices for camera sync and stopped recording
idx_start_rec = which(log$Behavior=="Started recording")
idx_stop_rec = which(log$Behavior=="Stopped recording")
idx_camera_sync = which(log$Behavior=="Camera Sync")

#Get recording length in sec
length_recording = log[idx_stop_rec,"Time"]/1000
#Important note: the length of the recording session should not be more than 2sec different between the
#deuteron log and "unit rasters". That is because in the script "log_read_sorted_NEX_file_CombinedChannels.m"
#Unit_rasters is initialized by the line: Unit_rasters(unit,:) = zeros(1,round(length_recording*temp_resolution)+1);
#round and +1 can lead to up to 2sec difference with the length of recording from the log.
time_before_cameraSync = log[2,"Time"]/1000
if (session == "Hooke_2021-08-02"){
  time_after_cameraSync = 0 #Camera ended after start of recording...
}else{
  time_after_cameraSync = log[idx_stop_rec,"Time"]/1000 - log[idx_camera_sync[2],"Time"]/1000
}
#Special note for Hooke_2021-08-02: There was no end recording time due to a bug.

#Estimate frame rate
length_video = (log[idx_camera_sync[2],"Time"] - log[idx_camera_sync[1],"Time"])/1000
frame_rate = 29.973 
#For Amos_2021-07-29: 29.791
# Hooke_2021-08-02: 29.97
# Amos_2021-08-03: 29.973
# Hooke_2021-08-05: 29.972
expected_num_frames = frame_rate*length_video
diff_frames = expected_num_frames-length(ME_right)
print(paste('Difference in frames between real and expected is:', diff_frames, 'frames.'))

#Pad ME with 0 in the extra recording spots
total_frames_needed = round(length_recording*frame_rate)
numFrames_before_cameraSync = round(time_before_cameraSync*frame_rate)
numFrames_after_cameraSync = round(time_after_cameraSync*frame_rate)
zero_pad_before = data.frame(matrix(data = 0, nrow = numFrames_before_cameraSync, ncol = 2)); names(zero_pad_before) = c("ME_left","ME_right")
zero_pad_after = data.frame(matrix(data = 0, nrow = numFrames_after_cameraSync, ncol = 2)); names(zero_pad_after) = c("ME_left","ME_right")

ME_final = rbind(zero_pad_before, MotionEnergy, zero_pad_after)

if (session=="Hooke_2021-08-02"){
  #Hooke_2021-08-02 has a longer video length than recording. To rectify this mismatch,
  #remove the extra frames to match the rneural recording length.
  idx_to_cut = round(frame_rate*abs(log[nrow(log),"Time"]/1000 - log[nrow(log)-1,"Time"]/1000))
  ME_final=ME_final[1:(nrow(ME_final)-idx_to_cut+1),]
}

#Save ME file
mac =1
if (mac == 1) { home = '~'} else
{home = 'C:/Users/GENERAL'}

setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/', sep="" ))
write.csv(ME_final, paste(session,'/ME.csv',sep=""),row.names=F)
