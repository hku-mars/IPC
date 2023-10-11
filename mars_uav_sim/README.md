# mars_uav_sim

# 1 安装依赖

安装[fmt](./fmt.zip)：解压后直接编译安装

```bash
cd fmt
mkdir build
cd build
cmake ..
make
sudo make install
```

# 2 运行

```bash
catkin_make
source devel/setup.zsh
roslaunch test_interface simulator.launch

```

`/depth_gene/depth` 节点中输出深度图

`/depth_gene/depth_pc` 节点中输出的是深度图中像素导出的点云。

Rviz中按`G`后点击目标，飞机会飞向点击位置。代码控制方法为：向`planning/pos_cmd`中发送控制期望值。
