import glob, os, re, sys
import random

import h5py
import numpy as np
from numpy import genfromtxt
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

#from bmtk.analyzer.cell_vars import plot_report


def original_connection_totals(synapses_file):
    data = genfromtxt(synapses_file,delimiter=' ')
    cell_types = ['EC','CA3e','CA3o','CA3b','DGg','DGh','DGb']
    total_cell_types = len(cell_types)
    cell_nums = [30,63,8,8,384,32,32]
    cell_start = [sum(cell_nums[:i]) for i,v in enumerate(cell_nums)]

    e_matrix = np.zeros((total_cell_types,total_cell_types))
    for i, row in enumerate(data):
        e_matrix[int(row[1]),int(row[3])]+=1

    return e_matrix, cell_types

def original_percent_connectivity(synapses_file):
    conn_totals, cell_types = original_connection_totals(synapses_file)
    
    cell_nums = [30,63,8,8,384,32,32]
    total_cell_types = len(cell_nums)

    max_connect = np.ones((total_cell_types,total_cell_types),dtype=np.float)

    for a, i in enumerate(cell_nums):
        for b, j in enumerate(cell_nums):
            max_connect[a,b] = i*j
    ret = conn_totals/max_connect
    ret = ret*100
    ret = np.around(ret, decimals=1)

    return ret, cell_types

def original_connection_divergence_average(synapses_file,convergence=False):
    conn_totals, cell_types = original_connection_totals(synapses_file)
    
    cell_nums = [30,63,8,8,384,32,32]
    total_cell_types = len(cell_nums)

    e_matrix = np.ones((total_cell_types,total_cell_types),dtype=np.float)

    for a, i in enumerate(cell_nums):
        for b, j in enumerate(cell_nums):
            c = b if convergence else a
            e_matrix[a,b] = conn_totals[a,b]/cell_nums[c]

    ret = np.around(e_matrix, decimals=1)

    return ret, cell_types

def plot_connection_info(data, labels, title):

    fig, ax = plt.subplots()
    im = ax.imshow(data)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)


    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(labels)):
        for j in range(len(labels)):
            text = ax.text(j, i, data[i, j],
                        ha="center", va="center", color="w")
    ax.set_ylabel('Source')
    ax.set_xlabel('Target')
    ax.set_title(title)
    fig.tight_layout()
    #plt.show()
    return fig, ax

def plot_raster(seed, files_dir='legacy', file_regex='_SpikeTime[0-6].txt'):
    """
    Author: Tyler Banks
    Loop through all spike files, create a list of when each cell spikes, plot.
    """
    #seed will be 667 for example
    files = glob.glob(os.path.join(files_dir,'bmtk_'+seed+file_regex))
    data = None
    
    cell_types = ['EC','CA3e','CA3o','CA3b','DGg','DGh','DGb']#,'CA3noise','DGnoise']
    cell_nums = [30,63,8,8,384,32,32]#,63,384]
    d = [[] for _ in range(sum(cell_nums))]
    
    color_picker=['red','orange','yellow','green','blue','purple','black']#,'orange','blue']
    colors = []
    offset=0
    files.sort() 
    for n, file in enumerate(files):
        data = genfromtxt(file,delimiter=',')
        for i, row in enumerate(data):
            d[int(row[0])+offset].append(row[1])
        offset+=cell_nums[n]-1
        
    for i, n in enumerate(cell_nums):
        for _ in range(n):
            colors.append(color_picker[i])
        
    fig, axs = plt.subplots(1,1)
    axs.eventplot(d,colors=colors)
    axs.set_title('bmtk_'+seed+file_regex)
    axs.set_ylabel('Cell Number')
    axs.set_xlabel('Time (ms)')
    axs.legend(cell_types[::-1])
    
    leg = axs.get_legend()
    for i,c in enumerate(color_picker):
        leg.legendHandles[-i-1].set_color(c)
    
    plt.savefig('raster.png')
    
    plt.show()
    
    return

def print_avg_frequencies(spikesfile='479_output/spikes.csv', population='hippocampus'):
    cell_names = ['EC','CA3e', 'CA3o', 'CA3b', 'DGg', 'DGh', 'DGb']
    cell_nums = [30,63,8,8,384,32,32]
    cell_start = np.cumsum(cell_nums)-cell_nums
    
    df = pd.read_csv(spikesfile,delimiter=' ')
    df = df[df['population'] == population]
    df = df[['timestamps','node_ids']]
    
    spikeraster = np.array(df)
    
    bins = np.append(cell_start, cell_start[-1]+cell_nums[-1])
    
    pltarr = plt.hist(spikeraster[:,1],bins=bins)
    freq = pltarr[0]/cell_nums
    
    print(dict(zip(cell_names,list(freq))))
    

def plot_positions(files_dir='Outputs', file_regex='008sepsLoc*.txt'):
    """
    Author: Tyler Banks
    """
    files = glob.glob(os.path.join(files_dir,file_regex))
    data = None
    
    cell_types = ['EC','CA3e','CA3o','CA3b','DGg','DGh','DGb','Feng']
    cell_nums = [30,63,8,8,384,32,32,27001]
    d = [[] for _ in range(sum(cell_nums))]
    
    color_picker=['red','blue','orange','yellow','green','black','purple','pink']
    colors = []
    offset=0
    
    fig = plt.figure()
    ax = Axes3D(fig)
    handles = []

    for n, file in enumerate(files):
        data = genfromtxt(file,delimiter=' ')
        handle = ax.scatter(data[:,0],data[:,1],data[:,2],color=color_picker[n],label=cell_types[n])
        handles.append(handle)
    
    plt.title('Hummos Hippocampus')
    plt.legend(handles=handles)
    plt.show()
    
    return

def plot_lfp(output_dir, lfp_file='ecp.h5',channel=0):

    ecp_path = os.path.join(output_dir,lfp_file)
    
    f = h5py.File(ecp_path)
    lfp = list(f['ecp']['data'])
    lfp_arr = np.asarray(lfp)

    plt.figure()
    plt.plot(lfp_arr[:,channel])
    plt.show()


def plot_firing_rates(seed, cellGroup, files_dir='legacy'):
    '''
    Author: Pete Canfield
    Generate some statistics about the firing rates of each cell in the network.

    Arguments:
        1) The seed that you want to analyze.
        2) What group of cells do you want to look at in the network
        3) The directory with the spiking information
    '''
    cellGroups = ['EC','CA3e', 'CA3o', 'CA3b', 'DGg', 'DGh', 'DGb']
    cell_nums = [30,63,8,8,384,32,32]
    groupIndex = cellGroups.index(cellGroup)


    #Get the files from the directory.
    file = open(os.path.join(files_dir,'bmtk_'+str(seed)+'_SpikeTime'+str(groupIndex)+'.txt'))


    #First check to see if all the cells in that cell fire.
    data = genfromtxt(file,delimiter=',', dtype=None)
    fireIds = [data[i][0] for i in range(data.size)]
    fireTimes = [data[i][1] for i in range(data.size)]


    dontFire = []
    fireFreqs = []

    fireIds.sort()


    startInd = 0
    ind = 0

    for id in range(cell_nums[groupIndex]):
        while ind < len(fireIds) and fireIds[ind] == id:
            ind += 1
            
        count = ind - startInd
        if count == 0:
            dontFire.append(id)
                                                
        fireFreqs.append(count/10.0) #divided by 10 becuase the simulation runtime is 10 seconds.
        startInd = ind

    print ("Max firing freqeuncy: " + str(max(fireFreqs)))

    print("Cells in group " + cellGroup + " that do not fire: ")
    print(dontFire)

    fig, (ax1,ax2) = plt.subplots(1,2)

    print("The average spiking frequency of " + cellGroup + " cells is: " + str(sum(fireFreqs)/len(fireFreqs)))
    
    ax2.hist(fireFreqs,int(1 + 3.322 * np.log(len(fireFreqs))))
    ax2.plot()
    ax2.set_title(cellGroup + " frequency distribution")
    ax2.set_xlabel("Frequency")
    ax2.set_ylabel("Number of " + cellGroup + " cells with reported average frequency")

    ax1.hist(fireTimes,1000)
    ax1.plot()
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Number of " + cellGroup + " cells active")
    ax1.set_title("Firing times of " + cellGroup + " cells (10 ms bins)")
    
    fig.tight_layout()
    plt.show()


def plot_instantaneous_frequencies(seed, cellGroup, cell_num = 2, files_dir='legacy'):
    '''
    Author: Pete Canfield
    Generate a n x n matrix which graphs the instantaneous frequencies of n^2 randomly selected
    cells in the specified cell group.

    Arguments:
        1) The seed that you want to analyze.
        2) What group of cells do you want to look at in the network
        3) The number of cells you want to sample (cell_num ^ 2)
        4) The directory with the spiking information

    '''

    cellGroups = ['EC','CA3e', 'CA3o', 'CA3b', 'DGg', 'DGh', 'DGb']
    cell_nums = [30,63,8,8,384,32,32]
    groupIndex = cellGroups.index(cellGroup)

    cellSquared = cell_num ** 2

    if cellSquared > cell_nums[groupIndex]:
        print("There arent enough cells in " + cellGroup + " to pick " + str(cellSquared) + " cells.")
        print("Exiting")
        return

    toGenerate = random.sample(range(cell_nums[groupIndex]),cellSquared)

    #Get the files from the directory.
    file = open(os.path.join(files_dir,'bmtk_'+str(seed)+'_SpikeTime'+str(groupIndex)+'.txt'))


    #First check to see if all the cells in that cell fire.
    data = genfromtxt(file,delimiter=',', dtype=None)

    #Create a dictionary where each key is the cell id we are intereted in 
    #each value is a pair of the time information and the corresponding instantaneous
    #frequency information.
    freqs = {}
    for i in toGenerate: 
        freqs[i] = ([0],[0])


    for elem in data:
        for id in toGenerate:
            if elem[0] == id:
                freqs[id][0].append(elem[1])
                freqs[id][1].append(1/((elem[1] - freqs[id][1][-1])/1000.0))

    fig, axs = plt.subplots(cell_num,cell_num, sharex=True,sharey=True)

    item = 0
    for i in range(cell_num):
        for j in range(cell_num):
            axs[i,j].plot(freqs[toGenerate[item]][0],freqs[toGenerate[item]][1])
            #axs[i,j].set_title("Cell " + str(toGenerate[item]))
            item += 1

    fig.suptitle("Instantaneous freqeuncy of randomly selected " + cellGroup + " cells")
    
    fig.text(0.5, 0.03, 'Time (ms)', ha='center')
    fig.text(0.03, 0.5, 'Frequncy (Hz)', va='center', rotation='vertical')
    
    #orginize the data by firing ids, then by firing times.
    data = np.sort(data,order=['f0','f1'])
    
    #Go through each ID id in the sorted data and extract a sublist for each cell.
    #sort the sublist by firing time and then extract instantaneous firing info from it

    #Find the maximum insatntaneous freqeuncy
    maxFreq = 0
    index = 0

    cellID = None
    firingStartTime = None
    firingEndTime = None

    while index < data.size - 1:
        nextID = data[index + 1][0]
        if nextID == data[index][0]:
            freq = 1/((data[index + 1][1] - data[index][1])/1000)
            if freq > maxFreq:
                maxFreq = freq
                firingStartTime = data[index][1]
                firingEndTime = data[index + 1][1]
                cellID = data[index][0]
            
        index += 1


    print("Cell with ID: " + str(cellID) + " fired at: " + str(firingStartTime) + " and " + str(firingEndTime) + " ms leading to a maxium frequency of : " + str(maxFreq) + " Hz.")

    plt.show()

    #Now lets generate the 
    

if __name__ == '__main__':

    #Example plot report https://github.com/AllenInstitute/bmtk/blob/develop/bmtk/analyzer/cell_vars.py
    #plot_report(config_file=None, report_file=None, report_name=None, variables=None, node_ids=None):
    #plot_report(config_file='100_output/100_simulation_config.json',node_ids=0)

    ### Hummos model ###   
    #data,labels = original_connection_totals('./Outputs/008sepsConnections.txt')
    #plot_connection_info(data, labels, 'Hummos Model Connection Totals')

    #data,labels = original_percent_connectivity('./008sepsConnections.txt')
    #plot_connection_info(data, labels, 'Hummos Model Connection Percentage')

    #data, labels = original_connection_divergence_average('./Outputs/008sepsConnections.txt',convergence=False)
    #plot_connection_info(data, labels, 'Hummos Model Average Divergence')

    #data, labels = original_connection_divergence_average('./Outputs/008sepsConnections.txt',convergence=True)
    #plot_connection_info(data, labels, 'Hummos Model Average Convergence')


    seed = sys.argv[1]
    #output_dir = seed + "_output"
    plot_raster(seed)
    

    #group = sys.argv[2]
    #plot_firing_rates(seed, group)

    #nCells = int(sys.argv[3])

    #(seed,group,nCells)
    #plot_lfp(output_dir)
    #print_average_frequencies('./legacy/'+seed+'_output/spikes.csv')
    
    #plot_instantaneous_frequencies(seed,group,nCells)

