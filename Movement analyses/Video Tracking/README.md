# Deep learning-based video tracking of natural monkey behavior

Recent advances in deep learning-based computer vision have enabled semi-automated tracking of animal behavior from video. We leverage these advances to conduct monkey detection, tracking, and pose estimation. We use the output of these models for neural analyses. 

The principal packages we use include: [mmdetection](https://mmdetection.readthedocs.io/en/latest/) for developing object detection and tracking models, [mmpose](https://mmpose.readthedocs.io/en/latest/) for developing pose estimation models, and [Roboflow ](https://roboflow.com/) for helpful labeling and visualization. We'd also like to acknowledge [DeepLabCut](http://www.mackenziemathislab.org/deeplabcut) and [SLEAP](https://sleap.ai/) as many analyses and visualizations were inspired by these prior works.

Finally, integrating computer vision-based tracking into (monkey) neuroscience experiments is a bit of pain. Please feel free reach out to the lead author, [Camille Testard](mailto:camille.testard94@gmail.com), or [Felipe Parodi](mailto:parodifelipe07@gmail.com), if you'd like help familiarizing yourself with this pipeline.

### Contents

Video Tracking/
  - README.md
  - video_demo.py
  - configs/
  - demo-vids/

### Usage
1. Create a conda environment following the [mmpose](https://mmpose.readthedocs.io/en/latest/) instructions.

2. Download the model configs from this [google drive link](https://drive.google.com/drive/folders/16greEMThJV7bKUbvTxQz4y6_BOXm91yl?usp=drive_link).

3. Activate your conda environment and run:
```
python video_demo.py \
  <detection-model-config> \
  <detection-model-weights> \
  <pose-model-config> \
  <pose-model-weights> \
  --input <input-video.mp4> \
  --device cuda:0 \
  --save-predictions \
  --output-root <output-root.json>
```
This will save a video with the visualized pose predictions and a json file containing the predictions for all detected instances.
