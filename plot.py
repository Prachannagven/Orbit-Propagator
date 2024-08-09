import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numpy as np


f = open("Values.txt")

x,y,z = [],[],[]

for line in f.readlines():
    a = list(map(float,line.split(" ")))
    x.append(a[0])
    y.append(a[1])
    z.append(a[2])

x = x[::10]
y = y[::10]
z = z[::10]

print(len(x))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
line, = ax.plot([], [], [], lw=0.2)
sat, = ax.plot3D(x[0],y[0],z[0],'green')
ax.set_xlim(min(x)-10,max(x)+10)
ax.set_ylim(min(y),max(y))
ax.set_zlim(min(z),max(z))

def update(i):
    line.set_data(x[0:i], y[0:i])
    line.set_3d_properties(z[0:i])
    sat.set_data((x[i],),(y[i],))
    sat.set_3d_properties((z[i],))
    return line

anim = animation.FuncAnimation(fig, update, frames=len(x), interval=1, )
line.set_data(x[:],y[:])
line.set_3d_properties(z[:])
plt.show()