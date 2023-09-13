import os
import sys
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot
from matplotlib import patches


def load_data(dname):
    glsizes = np.load(f"{dname}/glsizes.npy")
    lengths = np.load(f"{dname}/lengths.npy")
    ly, lx = lengths
    ny = glsizes[0]
    nx = glsizes[1]
    dx = lx / nx
    dy = ly / ny
    xc = np.linspace(0.5 * dx, lx - 0.5 * dx, nx)
    yc = np.linspace(0.5 * dy, ly - 0.5 * dy, ny)
    temp = np.load("{}/t.npy".format(dname))[:, 1:-1].T
    prs = np.load(f"{dname}/p_rs.npy")
    pys = np.load(f"{dname}/p_xs.npy")
    pxs = np.load(f"{dname}/p_ys.npy")
    return lx, ly, xc, yc, temp, prs, pxs, pys


def plot(ax, lx, ly, xs, ys, temp, prs, pxs, pys):
    ax.clear()
    ax.contourf(xs, ys, temp, vmin=-0.5, vmax=+0.5, levels=51, cmap="bwr")
    fc = "#AAAAAA"
    ec = "#000000"
    for pr, px, py in zip(prs, pxs, pys):
        c = patches.Circle(xy=(px,      py), radius=pr, fc=fc, ec=ec)
        ax.add_patch(c)
        c = patches.Circle(xy=(px - lx, py), radius=pr, fc=fc, ec=ec)
        ax.add_patch(c)
        c = patches.Circle(xy=(px + lx, py), radius=pr, fc=fc, ec=ec)
        ax.add_patch(c)
    kwrds = {
            "title": "",
            "aspect": "equal",
            "xlim": [0.0, lx],
            "ylim": [0.0, ly],
            "xlabel": "",
            "ylabel": "",
            "xticks": [],
            "yticks": [],
    }
    ax.set(**kwrds)


def main(ax, root):
    lx, ly, xs, ys, temp, prs, pxs, pys = load_data(root)
    plot(ax, lx, ly, xs, ys, temp, prs, pxs, pys)


if __name__ == "__main__":
    argv = sys.argv
    assert 2 == len(argv)
    root = argv[1]
    fig = pyplot.figure(figsize=(8, 4), frameon=False)
    ax = fig.add_axes(rect=[0.0, 0.0, 1.0, 1.0])
    if "step" in root:
        main(ax, root)
        pyplot.show(block=True)
    else:
        roots = [f"{root}/{dname}" for dname in os.listdir(root) if "step" in dname]
        roots = sorted(roots)
        for cnt, root in enumerate(tqdm(roots)):
            main(ax, root)
            pyplot.show(block=False)
            pyplot.pause(1.e-1)
            # pyplot.savefig(f"images/image{cnt:03d}.jpg", dpi=160)
    pyplot.close()
