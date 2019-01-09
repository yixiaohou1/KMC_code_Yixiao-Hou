# Import the Number Python model and use it as "np"
import numpy as np

# Set initial value lists
NA = [1000, 400, 4000, 500, 1000, 1000, 1000, 1000]  # Initial molecules of A
NB = [0, 0, 0, 500, 0, 0, 0, 0]  # Initial molecules of B
NC = [0]*8  # Initial molecules of C
T = [0]*8   # Time initial
K1 = [1, 1, 1, 1, 10, 1, 10, 1]  # Rate constants for A->B
K2 = [1, 1, 1, 1, 1, 10, 10, 1]  # Rate constants for B->C

n = 1  # Loops counting parameter

# Define variable sets for data logging and plotting
data = []  # Array for data to be logged into a .txt file latter
data_t = []  # Array for plotting x-axis of time
data_a = []  # Array for plotting y-axis of A molecules
data_b = []  # Array for plotting y-axis of B molecules
data_c = []  # Array for plotting y-axis of C molecules

while n <= 7:  # Repeat the calculation and plotting procedure for each dataset

    # Picking values from initial data lists
    n_a = NA[n-1]  # A molecules for KMC model
    n_b = NB[n-1]  # B molecules for KMC model
    n_c = NC[n-1]  # C molecules for KMC model
    t = T[n-1]  # Initial time (to be set zero)
    k1 = K1[n-1]  # Rate constant for A->B in KMC model
    k2 = K2[n-1]  # Rate constant for B->C in KMC model

    # Assign values for variables relative to analysis calculation
    k1_ = k1  # Rate constant for A->B in analytical model
    k2_ = k2  # Rate constant for B->C in analytical model

    es = open('KMC_data_%i.txt' % n, 'w+')  # Open a .txt file. If not existing, then create it. If existing, then erase it.
    log = open('KMC_data_%i.txt' % n, 'a+')  # Open a .txt file. If not existing, then create it. If existing, then append it.
    es.write("Data set: %i, \n" % n)  # Erase or create the file and append the data set number to it.
    log.write("Initial number of A molecules for KMC model: %i, \n" % NA[n - 1])  # Append the initial A molecules to the file.
    log.write("Initial number of B molecules for KMC model: %i, \n" % NB[n - 1])  # Append the initial B molecules to the file.
    log.write("Initial number of C molecules for KMC model: %i, \n" % NC[n - 1])  # Append the initial C molecules to the file.
    log.write("Constant for reaction A->B: %f, \n" % K1[n - 1])  # Append the rate constant of A->B to the file.
    log.write("Constant for reaction B->C: %f, \n\n" % K2[n - 1])  # Append the rate constant of B->C to the file.
    log.write("Time      NA   NB  NC Reaction \n")  # # Append the meaning of the data following to the file.

    # KMC model
    while n_a != 0 or n_b != 0:  # While there is still A and B left, there will be reaction occurring.
        p1 = np.random.random_sample()  # Get random number from 0 to 1.
        p2 = np.random.random_sample()  # Get random number from 0 to 1.
        r1 = k1 * n_a  # Definition for the reaction mechanism.
        r2 = k2 * n_b

        if p2 < (r1 / (r1 + r2)):  # Monte Carlo step
            n_a = n_a - 1  # Reaction A->B happens once.
            n_b = n_b + 1
            rxn = "A->B"  # Mark out the reaction that happened once for logging.
        else:
            n_b = n_b - 1  # Reaction B->C happens once.
            n_c = n_c + 1
            rxn = "B->C"  # Mark out the reaction that happened once for logging.

        data.append((t, n_a, n_b, n_c,rxn))  # Add the data generated in this loop to the dataset for logging.
        data_t.append(t)  # Add the data generated in this loop to the dataset for plotting.
        data_a.append(n_a)  # Add the data generated in this loop to the dataset for plotting.
        data_b.append(n_b)  # Add the data generated in this loop to the dataset for plotting.
        data_c.append(n_c)  # Add the data generated in this loop to the dataset for plotting.
        dt = -np.log(p1) / (r1 + r2)  # This is the time interval until the next reaction to occur.
        t = t + dt  # To refresh the time for next circulation.

    log.write('\n'.join('%f, %i, %i, %i, %s' % x for x in data))  # Log data of this loop to the .txt file.

    # Analytical model
    t_ = np.linspace(0, t, 10000)  # Set the array of time for analytical calculation. The maximum of time is same with the ending time of KMC model, for this amount of time would be enough.

    if k1_ == k2_:  # Judging which analytical solution to use.
        na = NA[n - 1] * np.exp(-1 * k1_ * t_)  # Calculating A molecule by analytical solution.
        nb = NB[n - 1] * np.exp(-1 * k2_ * t_) + NA[n - 1] * k1_ * t_ * np.exp(
            -k1_ * t_)  # Calculating B molecule by analytical solution.
        nc = (NA[n - 1] + NB[n - 1] + NC[n - 1]) - na - nb  # Calculating C molecule by analytical solution.
    else:
        na = NA[n - 1] * np.exp(-1 * k1_ * t_)  # Calculating A molecule by analytical solution.
        nb = NB[n - 1] * np.exp(-1 * k2_ * t_) + NA[n - 1] * k1_ / (k2_ - k1_) * (np.exp(-k1_ * t_) - np.exp(-k2_ * t_))  # Calculating B molecule by analytical solution.
        nc = (NA[n - 1] + NB[n - 1] + NC[n - 1]) - na - nb  # Calculating C molecule by analytical solution.

    # Plot with the data acquired above.
    import matplotlib.pyplot as plt  # Use the model named matplolib.pyplot to plot

    plt.plot(data_t, data_a, label="A molecules for KMC model")  # Configure for the line of A molecules of KMC model.
    plt.plot(data_t, data_b, label="B molecules for KMC model")  # Configure for the line of B molecules of KMC model.
    plt.plot(data_t, data_c, label="C molecules for KMC model")  # Configure for the line of C molecules of KMC model.
    plt.plot(t_, na, label="A molecules for analytical solution")  # Configure for the line of A molecules of analytical model.
    plt.plot(t_, nb, label="B molecules for analytical solution")  # Configure for the line of B molecules of analytical model.
    plt.plot(t_, nc, label="C molecules for analytical solution")  # Configure for the line of C molecules of analytical model.
    plt.xlim(0,t)  # Set the limit for time axis
    plt.ylim(0,max(NA[n-1],NB[n-1],NC[n-1]))  # Set the limit for molecules axis. The Largest amount of A, B and C initial molecules should work.
    plt.xlabel('Time(s)')  # Labeling the x-axis.
    plt.ylabel('Molecules')  # Labeling the y-axis.
    plt.title('KMC model data %i: NA=%i, NB=%i, K1=%i, K2=%i'% (n,NA[n-1],NB[n-1],K1[n-1],K2[n-1]))  # Labeling the graph with the dataset.
    plt.legend()
    plt.show()  # Show the plotted graph.

    # Return all the initial constants
    print("data: ", n)
    print("Initial number of A molecules for KMC model: ", NA[n-1])
    print("Initial number of B molecules for KMC model: ", NB[n-1])
    print("Initial number of C molecules for KMC model: ", NC[n-1])
    print("Constant for reaction A->B in KMC model: ", K1[n-1])
    print("Constant for reaction B->C in KMC model: ", K2[n-1],"\n\n")

    # Clear data lists for next circulation
    data.clear()
    data_t.clear()
    data_a.clear()
    data_b.clear()
    data_c.clear()

    # Add one to loop counting parameter before next loop begins
    n = n+1
