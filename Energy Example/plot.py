import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
import pyreadr
rc('font', family='serif')
rc('text', usetex=True)
rc('xtick',labelsize=14)
rc('ytick',labelsize=14)


def plot_and_save(data, fname, ylabel):
    data = np.reshape(data, [-1, 24])
    data_mean = np.mean(data, axis=0)
    data_std = np.std(data, axis=0)
    mask = data_std > data_mean
    data_std[mask] = 0.99 * data_mean[mask]
    data_csv = np.vstack((data_mean, data_std))
    np.savetxt('./csv/' + fname + '.csv', data_csv, delimiter=',')
    fig, ax = plt.subplots(1)
    ax.plot(np.arange(1, 25), data_mean, linewidth=2, color=[0, 0.4470, 0.7410])
    ax.plot(np.arange(1, 25), data_mean + data_std, linestyle='--', linewidth=2, color=[0, 0.4470, 0.7410])
    ax.plot(np.arange(1, 25), data_mean - data_std, linestyle='--', linewidth=2, color=[0, 0.4470, 0.7410])
    for i in range(len(data)):
        ax.plot(np.arange(1, 25), data[i], alpha=0.1, color=[0, 0.4470, 0.7410])
    ax.set_xlabel('Hour', fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_xlim(1, 24)
    ax.set_xticks([1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
    ax.grid(True)
    fig.savefig('./pdf/' + fname + '.pdf', format='pdf', dpi=1000, bbox_inches = 'tight', pad_inches = 0.02)
    plt.close(fig)
    return

def data_process(fname, start_date, end_date, ylabel):
    df = pyreadr.read_r(fname)[None]
    start_date = pd.to_datetime(start_date, format='%Y-%m-%d %H:%M:%S')
    end_date = pd.to_datetime(end_date, format='%Y-%m-%d %H:%M:%S')
    data_size = int((end_date - start_date).days * 24)
    df = df[(df['utc'] >= start_date) & (df['utc'] < end_date)]
    data = df.groupby(pd.Grouper(key='utc', freq='H')).sum()['sum'].values
    fname = fname[:-4].split('/')[-1]
    if data.size == data_size:
        plot_and_save(data, fname, ylabel)
    return


price_files = sorted(glob.glob(os.path.join('./price/', '*.csv')))
price = []
for fname in price_files:
    df = pd.read_csv(fname)
    price.extend(list(df['RRP'].values))
price = np.array(price)
price = np.reshape(price, [-1, 2])
price = np.mean(price, axis=1)
plot_and_save(price, 'price', r'Price (\$)')


battery_file = './battery/anonymous_public_battery_data.rds'
df = pyreadr.read_r(battery_file)[None]
data_csv = df.values.astype(float)
for data in data_csv:
    fname = './csv/battery_unit_' + str(int(data[0])) + '.csv'
    if data[1] > data[2]:
        np.savetxt(fname, data[1:], delimiter=',')

sort_files = lambda x: int(x.split('_')[-1].split('.')[0])
start_date = '2021-07-01 14:00:00'
end_date = '2021-10-01 14:00:00'


pv_files = sorted(glob.glob(os.path.join('./pv/', '*.rds')), key=sort_files)
ylabel_pv = 'Photovoltaic Power Generation (kW)'
for fname in pv_files:
    data_process(fname, start_date, end_date, ylabel_pv)


load_files = sorted(glob.glob(os.path.join('./load/', '*.rds')), key=sort_files)
ylabel_load = 'Power Consumption (kW)'
for fname in load_files:
    data_process(fname, start_date, end_date, ylabel_load)