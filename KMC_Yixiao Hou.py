# Import the Number Python model and use it as "np"
import numpy as np

# Convert variables into numpy version for the scientific calculation used later
n_a=np.int(1000)  # Initial number of A molecules for KMC model
n_b=np.int(0)  # Initial number of B molecules for KMC model
n_c=np.int(0)  # Initial number of C molecules for KMC model
t=np.float(0)  # Time elapsed
k1=np.float(1)  # Constant for reaction A→B in KMC model
k2=np.float(1)  # Constant for reaction B→C in KMC model
V=np.float(1000)  # Constant volume for analytical equations.
ca0=np.float(n_a/V)  # Initial concentration of A molecule
na=np.float(n_a)  # Instantaneous molecules for A
cb0=np.float(n_b/V)  # Initial concentration of B molecule
nb=np.float(n_b)  # Instantaneous molecules for B
cc0=np.float(n_c/V)  # Initial concentration of C molecule
nc=np.float(n_c)  # Instantaneous molecules for C
k1_=np.float(k1)  # Rate constant for reaction A→B
k2_=np.float(k2)  # Rate constant for reaction B→C

# Define variable sets for data logging and plotting
data=[]
data_t=[t]
data_a=[n_a]
data_b=[n_b]
data_c=[n_c]
data_na=[na]
data_nb=[nb]
data_nc=[nc]

import math  # Import the math Python model for exponential calculation

while n_a != 0 or n_b != 0:  # While there is still A and B left, there will be reaction occurring.
    p1 = np.random.random_sample()  # Get random number from 0 to 1.
    p2 = np.random.random_sample()
    r1 = k1 * n_a  # Definition for the reaction mechanism.
    r2 = k2 * n_b

    if k1_ = k2_ :
        na=V*ca0*math.exp(-1*k1_*t)  # Calculating A molecule by analytical solution.
        nb=V*(cb0*math.exp(-1*k2_*t)+ca0*k1*t*math.exp(-k1_*t))  # Calculating B molecule by analytical solution.
        nc=V*(ca0+cb0+cc0)-na-nb  # Calculating C molecule by analytical solution.
    else:


    if p2 < (r1/(r1+r2)):  # Monte Carlo step
        n_a = n_a-1  # Reaction A →B happens once.
        n_b = n_b+1
    else:
        n_b = n_b-1  # Reaction B→C happens once.
        n_c = n_c+1

    data.append((t, n_a, n_b, n_c,na,nb,nc))  # Add the data generated in this circulation to the dataset for logging.
    data_t.append(t)  # Add the data generated in this circulation to the dataset for plotting.
    data_a.append(n_a)
    data_b.append(n_b)
    data_c.append(n_c)
    data_na.append(na)
    data_nb.append(nb)
    data_nc.append(nc)
    dt = -np.log(p1) / (r1 + r2)  # This is the time interval until the next reaction to occur.
    t=t+dt  # To refresh the time for next circulation.
    open('KMC_data.txt','w+').write('\n'.join('%f, %i, %i, %i, %f, %f, %f' %x for x in data))
    # Record data of this circulation to the .txt file.

# Plot with the data acquired above.
import matplotlib.pyplot as plt  # Use the model named matplolib.pyplot to plot
plt.plot(data_t,data_a, label="A molecules")  # Configure for the first line.
plt.plot(data_t,data_b, label="B molecules")
plt.plot(data_t,data_c, label="C molecules")
plt.plot(data_t,data_na, label="A concentration")
plt.plot(data_t,data_nb, label="B concentration")
plt.plot(data_t,data_nc, label="C concentration")
plt.xlabel('t')  # Labeling the axis.
plt.ylabel('Molecules')
plt.title('KMC model')  # Labeling the graph.
plt.legend()
plt.show()  # Show the plotted graph.

# Return all the initial constants
print("Initial number of A molecules for KMC model: ", 1000)
print("Initial number of B molecules for KMC model: ", 0)
print("Initial number of C molecules for KMC model: ", 0)
print("Constant for reaction A→B in KMC model: ", 1)
print("Constant for reaction B→C in KMC model: ", 1)
print("Constant volume for analytical equation: ", 1000)
print("Initial concentration of A molecule: ", 1)
print("Initial concentration of B molecule: ", 0)
print("Initial concentration of C molecule: ", 0)
print("Rate constant for reaction A→B: ",1)
print("Rate constant for reaction B→C: ",1)
