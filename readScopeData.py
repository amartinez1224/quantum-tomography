import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numba import jit
import json
from tqdm import tqdm


@jit(nopython=True, parallel=True)
def gauss(x, amp, mu, sig):
    '''Numba gaussian function'''
    return amp*np.exp(-(x-mu)**2/(2.*sig**2))


@jit(nopython=True, parallel=True)
def seno(x, a1, a2, a3, a4):
    '''Numba sine function'''
    return a1*np.sin(a2*x+a3)+a4


def gausshist(n, bins):
    '''Fits a gaussian to a histogram'''
    # Get bins centers
    x = (np.abs(bins[1:]-bins[:-1])/2)+bins[:-1]
    # Perform gaussian fit over histogram
    arguments, variance = curve_fit(
        gauss, x, n, p0=[0.5, np.mean(x), np.std(x)])
    return arguments


def gen(data, density=400, v1=-12, v2=12):
    '''Generates histograms and fits from data'''
    # Generate voltage interval
    xd = np.linspace(v1, v2, density)
    # Generate matrix to store sinograms, each row is a sinogram
    m = np.empty((np.size(data[:, 0]), density))
    # Track gaussian data
    param = []
    files = tqdm(range(np.size(data[:, 0])))

    # For each data file
    for i in files:
        # Generate histogram
        n, bins = np.histogram(data[i, :], bins=10, density=True)
        # Fit histogram
        arguments = gausshist(n, bins)
        # Evaluate and append data to matrix
        m[i, :] = gauss(xd, arguments[0], arguments[1], arguments[2])

        param.append(arguments)
        files.set_description("Generating sinograms")
    return m, xd, np.array(param)


def loadData(fileNames):
    '''Loads data recorded by scope'''
    fileNames = tqdm(fileNames)
    data = []
    # Matrix of data loaded
    for file in fileNames:
        # must match data format
        data.append(np.genfromtxt(file, delimiter=","))
        fileNames.set_description("Loading %s" % file)
    return np.array(data)


def saveData(data, filename):
    '''Saves data in a json file'''
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)


def isVoltage(s):
    '''Checks if s is a valid voltage '''
    try:
        number = float(s)
        if number > -120 and number < 120:
            return True
        return False
    except ValueError:
        return False


def isAngle(s):
    '''Checks if s is a valid angle '''
    try:
        number = float(s)
        if number >= -10 and number < 10:
            return True
        return False
    except ValueError:
        return False


def main():
    print("\nWelcome to the scope processing tool!\n")

    print("\nNote: you will need to have the scope data saved in individual .csv files in a subdirectory named 'scope'. Each file should have a single column containing the values for the data measured for a single LO (piezo) angle. The angles for each file should be equidistant and the file name needs to be 'tekN.csv' where 'N' is a number from 0 to the total-1 number of measurements.\n")

    v1, v2 = -8, 8
    print("To begin, please enter the voltage range used in the scope during measurement.")
    while (voltage := input("Enter the voltage range (default -8,8): ")):
        if voltage == "":
            break
        elif voltage.count(',') == 1:
            v1I, v2I = voltage.split(',')
            print(v1I, v2I)
            if isVoltage(v1I) and isVoltage(v2I) and (v1I := float(v1I)) < (v2I := float(v2I)):
                v1 = v1I
                v2 = v2I
                break
            else:
                print("Invalid voltage range")
        else:
            print("Invalid voltage range")

    print(f'Voltage range set to {v1} to {v2} V\n')

    angle1, angle2 = 0, 2
    print("Please enter the angle range that was used for the measurements.")
    while (angle := input("Enter the angle range in PI multiples (default 0,2): ")):
        if angle == "":
            break
        elif angle.count(',') == 1:
            angle1I, angle2I = voltage.split(',')
            print(angle1I, angle2I)
            if isAngle(angle1I) and isAngle(angle2I) and (angle1I := float(angle1I)) < (angle2I := float(angle2I)):
                angle1 = angle1I
                angle2 = angle2I
                break
            else:
                print("Invalid angle range")
        else:
            print("Invalid angle range")

    print(f'Angle range set to {angle1} to {angle2} PI\n')

    angle1, angle2 = angle1*np.pi, angle2*np.pi

    measurements = 25
    while (measurementsI := input("Enter the number of measurements (default 25): ")):
        if measurementsI == "":
            break
        elif measurementsI.isdigit():
            measurements = int(measurementsI)
            break
        else:
            print("Invalid number of measurements")
    print(f'Number of measurements set to {measurements}\n')

    fileNames = ["scope/tek"+str(i)+".csv" for i in range(measurements)]

    data = loadData(fileNames)
    m, x, param = gen(data, v1=v1, v2=v2)

    phi = np.linspace(angle1, angle2, np.size(m[:, 0]))

    data = {
        "voltage range": f"{v1},{v2}",
        "angle range": f"{angle1},{angle2}",
        "phi": phi.tolist(),
        "x": x.tolist(),
        "pr": m.tolist()
    }
    saveData(data, "projections_from_scope.json")

    print("\nThe data has been generated and saved to the current directory\n")

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x, phi)
    ax.set_ylabel(r"$Piezo(\mu m)$")
    ax.set_zlabel(r"$P_{\theta}$")
    ax.plot_surface(X, Y, m, rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
    plt.show()


if __name__ == "__main__":
    main()
