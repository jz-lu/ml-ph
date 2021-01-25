from multiprocessing import Pool

NUM_AVAILABLE_CORES = 8

def prod(m, n):
    print("Multiplying m = %d with n = %d"%(m, n))
    return m * n

print('Constructing multiprocessing worker pool...')
pool = Pool(NUM_AVAILABLE_CORES) # Parallelize as much as possible.
print('Pool loaded.')

# Run asynchronous isolated calculations for each displacement over the pool.
print('Starting parallel computation...')
arr = [pool.apply_async(prod, args=(i, i+1)).get() for i in range(16)]

print(arr)
print(len(arr))