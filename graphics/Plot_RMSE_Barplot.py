#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size':16})

rmse_file = '/scratch/leeck/ABC-DA_1.5da_leeck_364/Plots/cycle-averaged_RMSE.txt'
n_runs = 5

u_data = []
v_data = []
w_data = []
rp_data = []
bp_data = []

with open(rmse_file, 'r') as f:
  for l in range(n_runs):
    line = f.readline()
    u_data.append(line.split()[-1])
  for l in range(n_runs):
    line = f.readline()
    v_data.append(line.split()[-1])
  for l in range(n_runs):
    line = f.readline()
    w_data.append(line.split()[-1])
  for l in range(n_runs):
    line = f.readline()
    rp_data.append(line.split()[-1])
  for l in range(n_runs):
    line = f.readline()
    bp_data.append(line.split()[-1])

index = ['u','v','w','r_prime','b_prime']
#ExpNames    = ['100% Be', '80% Be, 20% Bc', '50% Be, 50% Be', '100% Bc']
#ExpNames     = ['100% Be-EBV', '100% Be-EnSRF', '100% Bc']
ExpNames     = ['100% Be-30mem', '100% Be-20mem', '100% Be-10mem', '100% Bc']

rmse = {}
for i in range(n_runs-1):
  ctrl_u = float(u_data[-1])
  ctrl_v = float(v_data[-1])
  ctrl_w = float(w_data[-1])
  ctrl_rp = float(rp_data[-1])
  ctrl_bp = float(bp_data[-1])
  rmse[ExpNames[i]] = [float(u_data[i])/ctrl_u,float(v_data[i])/ctrl_v,float(w_data[i])/ctrl_w,float(rp_data[i])/ctrl_rp,float(bp_data[i])/ctrl_bp]

dataFrame = pd.DataFrame(data=rmse)
cols = dataFrame.columns.tolist()

# Optional reshuffling of colors
#tmp=[]; tmp.append(cols[1]); tmp.append(cols[3]); tmp.append(cols[2]); tmp.append(cols[0])
#tmp=[]; tmp.append(cols[1]); tmp.append(cols[2]); tmp.append(cols[0])
tmp=[]; tmp.append(cols[3]); tmp.append(cols[2]); tmp.append(cols[1]); tmp.append(cols[0])

dataFrame = dataFrame[tmp]
dataFrame.index = index
dataFrame = dataFrame[::-1]

#dataFrame.plot.barh(title='Cycle-averaged RMSE', color=['red','magenta','purple','blue'])
#dataFrame.plot.barh(title='Cycle-averaged RMSE', color=['green','orange','blue'])
dataFrame.plot.barh(title='Cycle-averaged RMSE', color=['red', 'cyan', 'brown', 'blue'])

plt.xlim(0.6,1.8)
plt.axvline(x=1, color='grey', linestyle='--', linewidth=2)
plt.xlabel('Ratio with respect to FreeBG')
plt.tight_layout()
plt.savefig('cycle-averaged_RMSE.png',dpi=300)
plt.show(block=True)

