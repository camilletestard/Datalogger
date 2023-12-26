"""
Script for running the mmpose video demo on videos with two monkeys.

Usage:
    Run from the command line as follows:

    python video_demo.py <detection-model-config> <detection-model-weights> <pose-model-config> <pose-model-weights> --input <input-video.mp4>

Requirements:
    Create a conda environment following the mmpose instructions.

Authors: Testard & Tremblay et al., 2023
"""

import os
from argparse import ArgumentParser

import cv2
import json_tricks as json
import mmcv
import numpy as np
from tqdm import tqdm

from mmdet.apis import inference_detector, init_detector
from mmdet.structures import DetDataSample
from mmengine.registry import MODELS as MMENGINE_MODELS
from mmpose.apis import inference_topdown
from mmpose.apis import init_model as init_pose_estimator
from mmpose.evaluation.functional import nms
from mmpose.registry import VISUALIZERS
from mmpose.structures import merge_data_samples
from mmpose.utils import adapt_mmdet_pipeline


NUM_MONKEYS = 2
DET_CAT_ID = 0


def split_instance(instances):
    """Convert instances into a list where each element is a dict that contains
    information about one instance."""
    results = []

    # return an empty list if there is no instance detected by the model
    if instances is None:
        return results

    for i in range(len(instances.keypoints)):
        result = dict(
            keypoints=instances.keypoints[i].tolist(),
            keypoint_scores=instances.keypoint_scores[i].tolist(),
        )
        if "bboxes" in instances:
            result["bbox"] = (instances.bboxes[i].tolist(),)
            if "bbox_scores" in instances:
                result["bbox_score"] = instances.bbox_scores[i]
            if "instance_id" in instances:
                result["instance_id"] = instances.instance_id[i]
        results.append(result)

    return results


def process_one_image(
    args,
    img,
    detector,
    pose_estimator,
    visualizer=None,
    tracker=None,
    show_interval=0,
    frame_id=0,
):
    """Visualize predicted keypoints (and heatmaps) of one image."""

    # predict bbox
    det_result = inference_detector(detector, img)
    if len(det_result.pred_instances) > NUM_MONKEYS:
        det_result.pred_instances = det_result.pred_instances[:NUM_MONKEYS]
    pred_instances = det_result.pred_instances[
        det_result.pred_instances.scores > args.bbox_thr
    ]
    monkey_indices = (pred_instances.labels == 0).nonzero(as_tuple=True)[0]
    pred_instances = pred_instances[monkey_indices]

    # track bounding boxes
    img_data_sample = DetDataSample()
    img_data_sample.pred_instances = pred_instances
    img_data_sample.set_metainfo(dict(frame_id=frame_id))
    pred_track_instances = tracker.track(img_data_sample)
    img_data_sample.pred_track_instances = pred_track_instances

    # format for pose:
    bboxes_cpu = pred_track_instances.bboxes.cpu().numpy()
    scores_cpu = pred_track_instances.scores.cpu().numpy()
    labels_cpu = pred_track_instances.labels.cpu().numpy()
    instance_ids_cpu = pred_track_instances.instances_id.cpu().numpy()

    # Ensure the scores and labels are in the right shape for concatenation
    scores_expanded = np.expand_dims(scores_cpu, axis=-1)
    labels_expanded = np.expand_dims(labels_cpu, axis=-1)
    instance_ids_expanded = np.expand_dims(instance_ids_cpu, axis=-1)

    # Concatenating along the second dimension
    bboxes = np.concatenate(
        [bboxes_cpu, scores_expanded, labels_expanded, instance_ids_expanded], axis=1
    )
    bboxes = bboxes[
        np.logical_and(labels_cpu == DET_CAT_ID, scores_cpu > args.bbox_thr)
    ]
    bboxes = bboxes[nms(bboxes, args.nms_thr), :4]
    # predict pose
    pose_results = inference_topdown(pose_estimator, img, bboxes)
    data_samples = merge_data_samples(pose_results)

    if isinstance(img, str):
        img = mmcv.imread(img, channel_order="rgb")
    elif isinstance(img, np.ndarray):
        img = mmcv.bgr2rgb(img)

    visualizer.add_datasample(
        "result",
        img,
        data_sample=data_samples,
        draw_gt=False,
        draw_heatmap=args.draw_heatmap,
        draw_bbox=args.draw_bbox,
        show_kpt_idx=args.show_kpt_idx,
        skeleton_style=args.skeleton_style,
        show=args.show,
        wait_time=show_interval,
        kpt_thr=args.kpt_thr,
    )

    # if there is no instance detected, return None
    return data_samples.get("pred_instances", None)


def setup_video_writer(output_file, video_reader):
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    frame_shape = video_reader[0].shape
    return cv2.VideoWriter(output_file, fourcc, 30, (frame_shape[1], frame_shape[0]))


def save_predictions(pred_save_path, dataset_meta, pred_instances_list):
    with open(pred_save_path, "w") as f:
        json.dump(
            dict(meta_info=dataset_meta, instance_info=pred_instances_list),
            f,
            indent="\t",
        )
    print(f"Predictions have been saved at {pred_save_path}")


def main():
    """Demo for monkey detection, tracking, and pose estimation."""
    parser = ArgumentParser()
    parser.add_argument("det_config", help="Config file for detection")
    parser.add_argument("det_checkpoint", help="Checkpoint file for detection")
    parser.add_argument("pose_config", help="Config file for pose")
    parser.add_argument("pose_checkpoint", help="Checkpoint file for pose")
    parser.add_argument("--input", type=str, default="", help="Input Video file")
    parser.add_argument(
        "--output-root",
        type=str,
        default="",
        help="root of the output img file. "
        "Default not saving the visualization images.",
    )
    parser.add_argument(
        "--save-predictions",
        action="store_true",
        default=False,
        help="whether to save predicted results",
    )
    parser.add_argument("--device", default="cuda:0", help="Device used for inference")
    parser.add_argument(
        "--bbox-thr", type=float, default=0.01, help="Bounding box score threshold"
    )
    parser.add_argument(
        "--nms-thr", type=float, default=0.9, help="IoU threshold for bounding box NMS"
    )
    parser.add_argument(
        "--kpt-thr", type=float, default=0.3, help="Visualizing keypoint thresholds"
    )
    parser.add_argument(
        "--radius", type=int, default=5, help="Keypoint radius for visualization"
    )
    parser.add_argument(
        "--thickness", type=int, default=1, help="Link thickness for visualization"
    )
    parser.add_argument(
        "--draw-bbox", action="store_true", help="Draw bboxes of instances"
    )

    args = parser.parse_args()

    assert args.input and args.det_config and args.det_checkpoint

    if args.output_root:
        output_file = os.path.join(args.output_root, os.path.basename(args.input))
    else:
        output_file = args.output

    # Ensure the directory for the output file exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    if args.save_predictions:
        args.pred_save_path = f"{args.output_root}/results_{os.path.splitext(os.path.basename(args.input))[0]}.json"

    obj_score_thrs_high = 0.4
    obj_score_thrs_low = 0.1
    init_track_thr = 0.45
    num_frames_retain = 30

    tracker = dict(
        type="mmdet.models.trackers.ByteTracker",
        motion=dict(type="KalmanFilter"),
        obj_score_thrs=dict(high=obj_score_thrs_high, low=obj_score_thrs_low),
        init_track_thr=init_track_thr,
        weight_iou_with_det_scores=True,
        match_iou_thrs=dict(high=0.1, low=0.5, tentative=0.3),
        num_frames_retain=num_frames_retain,
    )
    tracker = MMENGINE_MODELS.build(cfg=tracker)

    # build detector
    detector = init_detector(args.det_config, args.det_checkpoint, device=args.device)
    detector.cfg = adapt_mmdet_pipeline(detector.cfg)

    # build pose estimator
    pose_estimator = init_pose_estimator(
        args.pose_config,
        args.pose_checkpoint,
        device=args.device,
        cfg_options=dict(model=dict(test_cfg=dict(output_heatmaps=args.draw_heatmap))),
    )

    # build visualizer
    pose_estimator.cfg.visualizer.radius = args.radius
    pose_estimator.cfg.visualizer.alpha = args.alpha
    pose_estimator.cfg.visualizer.line_width = args.thickness
    visualizer = VISUALIZERS.build(pose_estimator.cfg.visualizer)
    visualizer.set_dataset_meta(
        pose_estimator.dataset_meta, skeleton_style=args.skeleton_style
    )

    video_reader = mmcv.VideoReader(args.input)
    video_writer = (
        None if not output_file else setup_video_writer(output_file, video_reader)
    )

    print(f"Successfully opened video file: {args.input}")

    pred_instances_list = []
    for frame_idx, frame in enumerate(tqdm(video_reader)):
        pred_instances = process_one_image(
            args, frame, detector, pose_estimator, visualizer, tracker, 0.001, frame_idx
        )
        if args.save_predictions:
            pred_instances_list.append(
                dict(frame_id=frame_idx, instances=split_instance(pred_instances))
            )
        if video_writer:
            video_writer.write(mmcv.rgb2bgr(visualizer.get_image()))

    if video_writer:
        video_writer.release()

    if args.save_predictions:
        save_predictions(
            args.pred_save_path, pose_estimator.dataset_meta, pred_instances_list
        )


if __name__ == "__main__":
    main()
