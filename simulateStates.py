import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numba import jit
import json


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


def gen(mus, sig, density=100, v1=-12, v2=12):
    '''Generates sinograms using mu and sigma'''
    m = np.zeros((np.size(mus), density))
    y = np.linspace(v1, v2, density)
    for i in range(np.size(mus)):
        r = np.random.normal(mus[i], sig[i], 10000)
        n, bins = np.histogram(r, bins=density, density=True)
        arguments = gausshist(n, bins)
        m[i, :] = gauss(y, arguments[0], arguments[1], arguments[2])
    return m, y


def isAngle(s):
    '''Checks if s is a valid angle '''
    try:
        number = float(s)
        if number >= 0 and number < 360:
            return True
        return False
    except ValueError:
        return False


def isVoltage(s):
    '''Checks if s is a valid voltage '''
    try:
        number = float(s)
        if number > -120 and number < 120:
            return True
        return False
    except ValueError:
        return False


def isAlpha(s):
    '''Checks if s is a valid alpha '''
    try:
        number = float(s)
        if number >= 0 and number <= 10:
            return True
        return False
    except ValueError:
        return False


def saveData(data, filename):
    '''Saves data in a json file'''
    with open(filename, 'w') as f:
        json.dump(data, f)


def main():
    print("\nWelcome to the quadrature distribution data generator!\n")
    v1, v2 = -8, 8
    print("To begin, please enter the voltage range of the simulated scope")
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
    print("Now, please select the state to simulate")
    while (mode := input("Choose state (1: Vacuum, 2: Coherent, 3: Thermal, 4: Squeezed): ")) not in ["1", "2", "3", "4"]:
        print("Invalid mode")
    print(f'State {mode} selected\n')
    alpha = 3
    if mode == "2" or mode == "4":
        print("Now, please enter the alpha value")
        while (alphaI := input("Choose alpha (0-10) (default 3): ")):
            if alphaI == "":
                break
            elif isAlpha(alphaI):
                alpha = float(alphaI)
                break
            else:
                print("Invalid alpha")
        print(f'Alpha set to {alpha}\n')

    phi = np.linspace(0, 2*np.pi, 100)
    if mode == "1":
        name = "vacuum"
        mus = np.zeros_like(phi)
        sig = np.zeros_like(phi)+(1/np.sqrt(2))
        pr, x = gen(mus, sig, v1=v1, v2=v2)

    elif mode == "2":
        name = "coherent"
        sig = np.zeros_like(phi)+(1/np.sqrt(2))
        mus = np.sin(phi)*alpha
        pr, x = gen(mus, sig, v1=v1, v2=v2)

    elif mode == "3":
        name = "thermal"
        mus = np.zeros_like(phi)
        sig = np.zeros_like(phi)+3
        pr, x = gen(mus, sig, v1=v1, v2=v2)

    elif mode == "4":
        name = "squeezed"
        eta = 2*np.pi/4
        sig = ((np.sin(2*(phi+(eta)))+1)*0.5 *
               (1/np.sqrt(2)))+(1/np.sqrt(2)) - 0.25
        mus = (np.sin(phi)+np.cos(phi))*alpha
        pr, x = gen(mus, sig, v1=v1, v2=v2)

    data = {"state": name,
            "voltage range": f"{v1},{v2}",
            "phi": phi.tolist(),
            "x": x.tolist(),
            "pr": pr.tolist()
            }
    saveData(data, "simulated_"+name+".json")

    print("\nThe data has been generated and saved to the current directory\n")

    view = input("If you want to see the generated data, enter 'y': ")
    if view == "y":
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i in range(np.size(phi)-1, 0, -1):
            ax.plot(x, np.zeros_like(x)+phi[i], pr[i, :], c='mediumvioletred')
        ax.set_ylabel(r"$\phi$")
        ax.set_xlabel(r"$x(v)$")
        ax.set_zlabel(r"$pr(x,\phi)$")
        ax.view_init(elev=0, azim=0)
        plt.show()
        plt.close()


if __name__ == "__main__":
    main()
