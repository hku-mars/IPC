import csv

import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import butter, lfilter, freqz


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def reader(master_table):
    file = './time_sudden_avoid.csv'
    df = pd.read_csv(file, encoding='GB18030')  # 如果不写encoding='GB18030'，则可能会出现中文乱码
    labels = list(df.columns.values)  # 获取所有字段
    print(labels)

    mapping = df['mapping']
    replan = df[' replan']
    sfc = df[' sfc']
    mpc = df[' mpc']

    total = mapping + replan + sfc + mpc

    # Filtering and plotting
    fig, ax = plt.subplots(1, 1, figsize=(4, 10))

    plt.gcf().subplots_adjust(bottom=0.1, left=0.3)

    font_size = 25
    tick_size = 20
    lw = 3
    # ax.boxplot(total)

    # 设置属性
    boxprops = {'color': '#76c893', 'linestyle': '-', 'linewidth': lw}
    whiskerprops = {'color': 'black', 'linestyle': '--', 'linewidth': lw}
    flierprops = {'color': 'r', 'marker': 'o', 'markersize': 5, 'linewidth': lw}
    medianprops = {'color': '#76c893', 'linestyle': '-', 'linewidth': lw}
    meanprops = {'color': '#168aad', 'linestyle': '-', 'linewidth': lw,'markersize': 15}
    widths = 0.5
    notch = True


    ax.boxplot(total,
               boxprops=boxprops,
               whiskerprops=whiskerprops,
               flierprops=flierprops,
               medianprops=medianprops,
               meanprops=meanprops,
               showmeans=True,
               showfliers=True,
               sym='o',
               widths=widths,
               notch=notch)

    plt.xticks([],fontsize=tick_size)
    plt.yticks([1,2,3,4,5,6,7,8,9,10],fontsize=tick_size)
    plt.xlabel('IPC', fontsize=font_size)
    plt.ylabel('Computation Time [ms]', fontsize=font_size)
    # plt.text(0.8, -1.3, 'Sudden Crossing\nObstacle Appear', fontsize=font_size, weight='bold', color='#f94144', ha='center')
    plt.grid()
    # plt.show()
    plt.savefig('sudden_obs.svg',dpi=1000)


plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

if __name__ == "__main__":
    reader('test')
