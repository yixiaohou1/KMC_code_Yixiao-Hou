# Import the Number Python model and use it as "np"
import numpy as np

# Set initial value lists
NA = [1000, 400, 4000, 500, 1000, 1000, 1000, 1000]
NB = [0, 0, 0, 500, 0, 0, 0, 0]
NC = [0]*8
T = [0]*8
K1 = [1, 1, 1, 1, 10, 1, 10, 1]
K2 = [1, 1, 1, 1, 1, 10, 10, 1]
V = [1000]*7 + [1]


# Set variables into numpy version for the scientific calculation used later
n_a = np.int()  # Initial number of A molecules for KMC model
n_b = np.int()  # Initial number of B molecules for KMC model
n_c = np.int()  # Initial number of C molecules for KMC model
t = np.float()  # Time elapsed
k1 = np.float()  # Constant for reaction A→B in KMC model
k2 = np.float()  # Constant for reaction B→C in KMC model
v = np.float()  # Constant volume for analytical equations.
ca0=np.float()  # Initial concentration of A molecule
na=np.float()  # Instantaneous molecules for A
cb0=np.float()  # Initial concentration of B molecule
nb=np.float()  # Instantaneous molecules for B
cc0=np.float()  # Initial concentration of C molecule
nc=np.float()  # Instantaneous molecules for C
k1_=np.float()  # Rate constant for reaction A→B
k2_=np.float()  # Rate constant for reaction B→C
n = 1  # Circulation counting parameter

# Define variable sets for data logging and plotting
data = []
data_t = [t]
data_a = [n_a]
data_b = [n_b]
data_c = [n_c]
data_na = [na]
data_nb = [nb]
data_nc = [nc]

import math  # Import the math Python model for exponential calculation

while n <= 8:  # Repeat the calculation and plotting procedure for each dataset

    # Picking values from initial data lists
    n_a = NA[n-1]
    n_b = NB[n-1]
    n_c = NC[n-1]
    t = T[n-1]
    k1 = K1[n-1]
    k2 = K2[n-1]
    v = V[n-1]

    # Assign values for variables relative to analysis calculation
    ca0 = n_a / v
    na = n_a
    cb0 = n_b / v
    nb = n_b
    cc0 = n_c / v
    nc = n_c
    k1_ = k1
    k2_ = k2

    while n_a != 0 or n_b != 0:  # While there is still A and B left, there will be reaction occurring.
        p1 = np.random.random_sample()  # Get random number from 0 to 1.
        p2 = np.random.random_sample()
        r1 = k1 * n_a  # Definition for the reaction mechanism.
        r2 = k2 * n_b

        if k1_ == k2_ :
            na = v * ca0 * math.exp(-1 * k1_ * t)  # Calculating A molecule by analytical solution.
            nb = v * (cb0 * math.exp(-1 * k2_ * t) + ca0 * k1 * t * math.exp(-k1_ * t))  # Calculating B molecule by analytical solution.
            nc = v * (ca0 + cb0 + cc0) - na - nb  # Calculating C molecule by analytical solution.
        else:
            na = v * ca0 * math.exp(-1 * k1_ * t)
            nb = v * (cb0 * math.exp(-1 * k2_ * t) + ca0 * k1 / (k2 - k1) * (math.exp(-k1_ * t) - math.exp(-k2_ * t)))
            nc = v * (ca0 + cb0 + cc0) - na - nb

        if p2 < (r1 / (r1 + r2)):  # Monte Carlo step
            n_a = n_a - 1  # Reaction A →B happens once.
            n_b = n_b + 1
        else:
            n_b = n_b - 1  # Reaction B→C happens once.
            n_c = n_c + 1

        data.append((t, n_a, n_b, n_c, na, nb, nc))  # Add the data generated in this circulation to the dataset for logging.
        data_t.append(t)  # Add the data generated in this circulation to the dataset for plotting.
        data_a.append(n_a)
        data_b.append(n_b)
        data_c.append(n_c)
        data_na.append(na)
        data_nb.append(nb)
        data_nc.append(nc)
        dt = -np.log(p1) / (r1 + r2)  # This is the time interval until the next reaction to occur.
        t = t + dt  # To refresh the time for next circulation.
        open('KMC_data_%i.txt' % n, 'w+').write('\n'.join('%f, %i, %i, %i, %f, %f, %f' % x for x in data))  # Record data of this circulation to the .txt file.

    # Plot with the data acquired above.
    import matplotlib.pyplot as plt  # Use the model named matplolib.pyplot to plot

    plt.plot(data_t, data_a, label="A molecules")  # Configure for the first line.
    plt.plot(data_t, data_b, label="B molecules")
    plt.plot(data_t, data_c, label="C molecules")
    plt.plot(data_t, data_na, label="A concentration")
    plt.plot(data_t, data_nb, label="B concentration")
    plt.plot(data_t, data_nc, label="C concentration")
    plt.xlabel('t')  # Labeling the axis.
    plt.ylabel('Molecules')
    plt.title('KMC model data %i'% n)  # Labeling the graph.
    plt.legend()
    plt.show()  # Show the plotted graph.

# Return all the initial constants
    print("data: ", n)
    print("Initial number of A molecules for KMC model: ", NA[n-1])
    print("Initial number of B molecules for KMC model: ", NB[n-1])
    print("Initial number of C molecules for KMC model: ", NC[n-1])
    print("Constant for reaction A→B in KMC model: ", K1[n-1])
    print("Constant for reaction B→C in KMC model: ", K2[n-1])
    print("Constant volume for analytical equation: ", V[n-1])
    print("Initial concentration of A molecule: ", NA[n-1]/v)
    print("Initial concentration of B molecule: ", NB[n-1]/v)
    print("Initial concentration of C molecule: ", NC[n-1]/v)
    print("Rate constant for reaction A→B: ", k1_)
    print("Rate constant for reaction B→C: ", k2_, "\n\n")

    # Clear data lists for next circulation
    data.clear()
    data.clear()
    data_t.clear()
    data_a.clear()
    data_b.clear()
    data_c.clear()
    data_na.clear()
    data_nb.clear()
    data_nc.clear()

    # Add one to circulation counting parameter before next circulation begins
    n=n+1
