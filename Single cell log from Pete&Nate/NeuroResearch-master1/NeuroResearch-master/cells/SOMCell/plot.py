from neuron import h
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

h.load_file('stdrun.hoc')
h.load_file('cell_plotter.hoc')

def i_inj_plot(t):
    if t <= 400 or t >= 900:
        return 0
    else:
        return 0.3

trace = h.membrane_voltage.as_numpy()
time = h.time.as_numpy()
time = time[10000:40000]
trace = trace[10000:40000]

inj = np.array([i_inj_plot(t) for t in time])

current_inj = h.current_inj.as_numpy()
fi_target = h.target.as_numpy()
fi_reported = h.reported.as_numpy()


fig = plt.figure(figsize = (4,3.5), constrained_layout = False)
gs = mpl.gridspec.GridSpec(3, 1, height_ratios=[4,4,1])

ax = plt.subplot(gs[0])
plt.title("Simulated SOM+ Cell")
plt.plot(current_inj, fi_reported, label = 'Simulated', lw = 0.9)
plt.plot(current_inj, fi_target, label = 'Target', lw = 0.9)
plt.legend(loc='lower right', fontsize='xx-small', frameon=False)
plt.ylabel('frequency (Hz)')
plt.xlabel('current (nA)')


ax = plt.subplot(gs[1])


plt.plot(time, trace, color = 'red', lw = 0.9)
plt.ylabel('voltage (mV)')
ax.set_xticks([])
ax.set_yticks([-70, -20, 40])
ax.spines['bottom'].set_visible(False)


ax = plt.subplot(gs[2])
ax.spines['top'].set_visible(False)
plt.plot(time, inj, color = 'black', lw = 0.9)
plt.ylabel('current (nA)')
plt.xlabel('time (ms)')
ax.set_yticks([0,0.3])
ax.yaxis.set_label_position('right')
ax.yaxis.set_ticks_position('right') 
plt.tight_layout()


plt.show()