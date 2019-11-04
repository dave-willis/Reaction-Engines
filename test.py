#import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

n = len(np.arange(5, 11))
phi = len(np.arange(0.3, 0.9, 0.1))
psi = len(np.arange(0.8, 2, 0.1))
Lambda = len(np.arange(0.45,0.55,0.05))
AR = len(np.arange(1.0, 2.0, 0.2))
dho = len(np.arange(1, 1.2, 0.02))
a1 = len(np.arange(0, 20, 2))

num = AR*phi**2*psi**2*dho*a1*n

print(num)

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
t = np.arange(0.0, 1.0, 0.001)
a0 = 5
f0 = 3
s = a0 * np.sin(2 * np.pi * f0 * t)
l, = plt.plot(t, s, lw=2)
ax.margins(x=0)
ax.set_ybound(np.amin(a0*np.sin(2*np.pi*f0*t)), np.amax(a0*np.sin(2*np.pi*f0*t)))
ax.set_xbound(np.amin(t), np.amax(t))

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

sfreq = Slider(axfreq, 'Freq', 0.1, 30.0, valinit=f0)
samp = Slider(axamp, 'Amp', 0.1, 10.0, valinit=a0)


def update(val):
    amp = samp.val
    freq = sfreq.val
    l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    ax.set_ybound(np.amin(amp*np.sin(2*np.pi*freq*t)), np.amax(amp*np.sin(2*np.pi*freq*t)))
    ax.set_xbound(np.amin(t), np.amax(t))
    fig.canvas.draw_idle()


sfreq.on_changed(update)
samp.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sfreq.reset()
    samp.reset()
button.on_clicked(reset)


plt.show()
