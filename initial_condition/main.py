import os
import sys
import numpy as np


def init_time(dest):
    # iterator and time
    step = np.array(0, dtype=np.uint64)
    time = np.array(0, dtype=np.float64)
    np.save(f"{dest}/step.npy", step)
    np.save(f"{dest}/time.npy", time)
    return


def init_domain(lengths, glsizes, dest):
    np.save(f"{dest}/glsizes.npy", np.array(glsizes, dtype=np.uint64))
    np.save(f"{dest}/lengths.npy", np.array(lengths, dtype=np.float64))
    return


def init_fluid(lengths, glsizes, dest):
    shape0 = (glsizes[2], glsizes[1], glsizes[0] + 1)
    shape1 = (glsizes[2], glsizes[1], glsizes[0] + 2)
    ux = np.zeros(shape0, dtype=np.float64)
    uy = np.zeros(shape1, dtype=np.float64)
    uz = np.zeros(shape1, dtype=np.float64)
    p  = np.zeros(shape1, dtype=np.float64)
    t  = -0.5 + np.random.random_sample(shape1)
    np.save(f"{dest}/ux.npy", ux)
    np.save(f"{dest}/uy.npy", uy)
    np.save(f"{dest}/uz.npy", uz)
    np.save(f"{dest}/p.npy", p)
    np.save(f"{dest}/t.npy", t)
    return


def init_particles(lengths, glsizes, dest):
    nitems = 1
    rs = [0.25]
    ds = [1.]
    xs = [0.4]
    ys = [0.5]
    zs = [0.5]
    uxs = [0.]
    uys = [0.]
    uzs = [0.]
    vxs = [0.]
    vys = [0.]
    vzs = [0.]
    np.save(f"{dest}/p_nitems.npy", np.uint64(nitems))
    np.save(f"{dest}/p_rs.npy", np.float64(rs))
    np.save(f"{dest}/p_ds.npy", np.float64(ds))
    np.save(f"{dest}/p_xs.npy", np.float64(xs))
    np.save(f"{dest}/p_ys.npy", np.float64(ys))
    np.save(f"{dest}/p_zs.npy", np.float64(zs))
    np.save(f"{dest}/p_uxs.npy", np.float64(uxs))
    np.save(f"{dest}/p_uys.npy", np.float64(uys))
    np.save(f"{dest}/p_uzs.npy", np.float64(uzs))
    np.save(f"{dest}/p_vxs.npy", np.float64(vxs))
    np.save(f"{dest}/p_vys.npy", np.float64(vys))
    np.save(f"{dest}/p_vzs.npy", np.float64(vzs))
    return


def main():
    argv = sys.argv
    assert 2 == len(argv)
    lengths = [
            float(os.environ["lx"]),
            float(os.environ["ly"]),
            float(os.environ["lz"]),
    ]
    glsizes = [
            int(os.environ["glisize"]),
            int(os.environ["gljsize"]),
            int(os.environ["glksize"]),
    ]
    dest = argv[1]
    init_time(dest)
    init_domain(lengths, glsizes, dest)
    init_fluid(lengths, glsizes, dest)
    init_particles(lengths, glsizes, dest)
    return


main()
