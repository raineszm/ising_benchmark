import multiprocessing as mp

import numpy as np

from .simulation import make_simulation


def worker(input, output, N, n_evolve, n_average):
    sim = make_simulation(N)
    for t in iter(input.get, 'STOP'):
        r1 = np.random.randint(0, N, size=(n_evolve + n_average, 2))
        r2 = np.random.rand(n_evolve + n_average)
        (M, U) = sim.ensemble_av(1 / t, n_evolve, n_average, r1, r2)
        output.put((t, M, U))


def psimulate(fname, ts, N, n_evolve=1000, n_average=1000):
    t_queue = mp.Queue()
    out_queue = mp.Queue()

    args = (t_queue, out_queue, N, n_evolve, n_average)

    threads = [mp.Process(target=worker, args=args) for _ in range(4)]

    for t in threads:
        t.start()

    for t in ts:
        t_queue.put(t)

    for _ in threads:
        t_queue.put('STOP')

    with open(fname, 'w') as f:
        for _ in range(len(ts)):
            print('{},{},{}'.format(*out_queue.get()), file=f)

    for t in threads:
        t.join()
