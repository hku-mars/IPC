import rosbag
import numpy as np
import open3d as o3d
import sensor_msgs.point_cloud2 as pc2
import ctypes
import struct
import argparse
import time
from tqdm import tqdm, trange
import os

parser = argparse.ArgumentParser(description='Process bag sensor_msg/PointCloud2 to PCD file.')
parser.add_argument('--bag_name', type=str,
                    help='the bag file name')
parser.add_argument('--topic_name', type=str,
                    help='the point cloud2 topic name')
bagfile = parser.parse_args().bag_name
split_list = os.path.basename(bagfile)
bag_name = split_list.split(".",2)[0]
parser.add_argument('--outfile_name', type=str,default=bag_name+".pcd",
                    help='the point cloud2 topic name')
args = parser.parse_args()
topic_name = args.topic_name
bag = rosbag.Bag(bagfile)
total_len = 0
for topic, msg, t in bag.read_messages(topics=[topic_name]):
    rgb = np.array([[0, 0, 0]])
    gen = pc2.read_points(msg, skip_nans=True)
    int_data = list(gen)
    total_len = total_len + len(int_data)
xyz = np.empty((total_len, 3))
rgb = np.empty((total_len, 3))
cnt = 0

with tqdm(total=total_len, desc='Process', leave=True, ncols=100, unit='pt', unit_scale=True) as pbar:
    for topic, msg, t in bag.read_messages(topics=[topic_name]):
        gen = pc2.read_points(msg, skip_nans=True)
        int_data = list(gen)
        for x in int_data:
            test = x[3]
            # cast float32 to int so that bitwise operations are possible
            s = struct.pack('>f' ,test)
            i = struct.unpack('>l',s)[0]
            # you can get back the float value by the inverse operations
            pack = ctypes.c_uint32(i).value
            r = (pack & 0x00FF0000)>> 16
            g = (pack & 0x0000FF00)>> 8
            b = (pack & 0x000000FF)
            # prints r,g,b values in the 0-255 range
            # x,y,z can be retrieved from the x[0],x[1],x[2]
            xyz[cnt,:] =[x[0],x[1],x[2]]# np.append(xyz,[[x[0],x[1],x[2]]], axis = 0)
            cnt = cnt+1
            pbar.update(1)


out_pcd = o3d.geometry.PointCloud()
out_pcd.points = o3d.utility.Vector3dVector(xyz)
o3d.io.write_point_cloud(args.outfile_name,out_pcd)
bag.close()