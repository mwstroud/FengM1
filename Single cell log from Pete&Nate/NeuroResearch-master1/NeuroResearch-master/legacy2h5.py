import os, sys
import numpy as np
import h5py

def convert(input_dir, output_dir, seed = '0041'):

    cell_offsets = np.cumsum(np.array([30,63,8,8,384,32,32])) - 30
    

    f_out = h5py.File(output_dir,'a')
    f_out.create_group('spikes/hippocampus')
    hipp = f_out['spikes/hippocampus']

    fireIds = np.empty((0,),dtype=np.uint64)
    fireTimes = np.empty((0,),dtype=np.dtype('<f8'))
    for i in range(7):
        file = open(os.path.join(input_dir,seed+"SpikeTime"+str(i)+".txt"))

        data = np.genfromtxt(file,delimiter=',', dtype=None)

        fireIds = np.append(fireIds,np.array([data[j][0] + cell_offsets[i] for j in range(data.size)]))
        fireTimes = np.append(fireTimes,np.array([data[j][1] for j in range(data.size)]))

    hipp.create_dataset('node_ids',dtype=np.uint64,shape=(fireIds.size,))
    hipp.create_dataset('timestamps',dtype=np.dtype('<f8'),shape=(fireTimes.size,))

    hipp['node_ids'].write_direct(fireIds)
    hipp['timestamps'].write_direct(fireTimes)





if __name__ == '__main__':
    convert(sys.argv[1],sys.argv[2])
