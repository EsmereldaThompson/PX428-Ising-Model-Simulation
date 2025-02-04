import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

# Animation all adapted from a video by Numeryst

# AIM
# to create a plotter that creates an animated heatmap to show the evolution of the spins as the metropolis algorithm iterates
# if plots every iteration is going to be very slow, so will need to reduce the data slightly so the change is still smooth without being super slow
print("Running successfully")

N_dim = 2
N_spins = 20

power = 3
N_iter = np.power(2,power) * np.power(N_spins,N_dim)

path = "/home/ewt/UNI Y3/PX442 - Laboratory for Mathematics & Physics Students/Simulation Lab/output_spins_file.txt" # had numerous issues if the file path was not fully specified for some reason?
spins_file = open(path, "r") # loads the output file from c

data = np.empty((0,N_spins,N_spins), dtype = int) # creating an empty array with appropriate dimensions
spins_current_iteration = np.empty((0,N_spins), dtype = int)

for line in spins_file:
    if line.rstrip(): # formatted the c file to have a blank line between iteration so this flags when a new iteration has started
        line = line.rstrip()
        #line.replace("\n"," ") # removing the \n from the end of the line
        spins_current_line = np.array([line.split(",")]) # storing all the information on the line as separate items in the array
        spins_current_iteration = np.append(spins_current_iteration,spins_current_line, axis = 0)
    else:
        data = np.append(data,[spins_current_iteration], axis = 0)
        spins_current_iteration = np.empty((0,N_spins), dtype = int)

data = data.astype(int)
shape = np.shape(data)
frame_num = shape[0]
print("File extracted successfully")

figure, axis = plt.subplots() # creating the figure
axis.set_title(f'Particle spins over {N_iter} iterations of the Metropolis algorithm')
axis.set_xticklabels([])
axis.set_yticklabels([])

heatmap = axis.imshow(data[0], cmap='Greys', extent = (-N_spins,N_spins,-N_spins,N_spins)) # Black and white heatmap
# essentially creates an empty 2d raster to animate in the next steps

def update(frame):
    heatmap.set_data(data[frame]) # should update to the next iteration
    return [heatmap]

animation = FuncAnimation(figure, update, frames = frame_num, blit = True, interval = 50, repeat = False)
animation.save("Metropolis_algorithm_animation.gif")
plt.show()

print("All Done :)")
#print(spins)
#spins_file.close()