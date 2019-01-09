# Import the Number Python model and use it as "np"
import numpy as np

# Define variables and allow users to input them
n_a=input("Initial A molecules: ")
n_b=input("Initial B molecules: ")
n_c=input("Initial C molecules: ")
k1=input("Rate constant for A→B: ")
k2=input("Rate constant for B→C: ")

# Convert variables into numpy version for the scientific calculation used later
n_a=np.int(n_a)
n_b=np.int(n_b)
n_c=np.int(n_c)
t=np.float(0)
k1=np.float(k1)
k2=np.float(k2)

# Define variable sets for data logging and ploting
data=[]
data_t=[t]
data_a=[n_a]
data_b=[n_b]
data_c=[n_c]

while n_a != 0 or n_b != 0:  # While there is still A and B left, there will be reaction occurring.
    p1 = np.random.random_sample()  # Get random number from 0 to 1.
    p2 = np.random.random_sample()
    r1 = k1 * n_a  # Definition for the reaction mechanism.
    r2 = k2 * n_b
    if p2 < (r1/(r1+r2)):  # Monte Carlo step
        n_a = n_a-1  # Reaction A →B happens once.
        n_b = n_b+1
    else:
        n_b = n_b-1  # Reaction B→C happens once.
        n_c = n_c+1
    data.append((t, n_a, n_b, n_c))  # Add the data generated in this circulation to the dataset for logging.
    data_t.append(t)  # Add the data generated in this circulation to the dataset for plotting.
    data_a.append(n_a)
    data_b.append(n_b)
    data_c.append(n_c)
    dt = -np.log(p1) / (r1 + r2)  # This is the time interval until the next reaction to occur.
    t=t+dt  # To refresh the time for next circulation.
    open('KMC_data.txt','w+').write('\n'.join('%f, %i, %i, %i' %x for x in data))
    # Record data of this circulation to the .txt file.

# Plot with the data acquired above.
import matplotlib.pyplot as plt  # Use the model named matplolib.pyplot to plot
plt.plot(data_t,data_a, label="A molecules")  # Configure for the first line.
plt.plot(data_t,data_b, label="B molecules")
plt.plot(data_t,data_c, label="C molecules")
plt.xlabel('t')  # Labeling the axis.
plt.ylabel('Molecules')
plt.title('KMC model')  # Labeling the graph.
plt.show()  # Show the plotted graph. 