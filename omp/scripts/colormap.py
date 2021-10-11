import numpy as np
from matplotlib import pyplot as plt, animation
import argparse

parser = argparse.ArgumentParser(description='Generate heatmap gif')
parser.add_argument('--frames', dest='frames', type=int, help='Amount of frames (please, use the same number as when running heat)')
parser.add_argument('--input', dest='input', type=str, help='Path to frames file')

args = parser.parse_args()

f = open(args.input, "r")
def read_frame():
    frame = []
    s = f.readline()
    while s != "\n":
        vals = [float(x) for x in s.split()]
        if len(vals) != 0:
            frame.append(vals)
        s = f.readline()
    return frame

def animate(i):
    plt.clf()

    plt.title(f"Temperature at t = {i} unit time")
    plt.xlabel("x")
    plt.ylabel("y")

    # This is to plot u_k (u at time-step k)
    plt.pcolormesh(read_frame(), cmap=plt.cm.jet, vmin=0, vmax=100)
    plt.colorbar()

    return plt

anim = animation.FuncAnimation(plt.figure(), animate, interval=50, frames=args.frames - 1)
anim.save('heat.gif')
f.close()