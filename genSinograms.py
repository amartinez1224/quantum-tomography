import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
import numba
from numba import jit
from tqdm import tqdm
from scipy.special import factorial
import matplotlib.colors as colors
import matplotlib.cm as cm

# Gaussian function
@jit(nopython=True, parallel=True)
def gauss(x, amp, mu, sig):
    return amp*np.exp(-(x-mu)**2/(2.*sig**2))

#Sine function
@jit(nopython=True, parallel=True)
def seno(x, a1, a2, a3, a4):
    return a1*np.sin(a2*x+a3)+a4

# Gaussian fits from histogram
def gausshist(n, bins):
    # Get bins centers
    x = (np.abs(bins[1:]-bins[:-1])/2)+bins[:-1]
    # Perform gaussian fit over histogram
    arguments, variance = curve_fit(gauss, x, n, p0=[0.5,np.mean(x),np.std(x)])
    return arguments

# Loads data recorded by scope
def loadData(fileNames):
    fileNames = tqdm(fileNames)
    data = []
    # Matrix of data loaded
    for file in fileNames:
        data.append(np.genfromtxt(file,skip_header=15, delimiter=",",usecols=5)) # optional arguments must match data format
        fileNames.set_description("Loading %s" % file)
    return np.array(data)

# Generates histograms and fits from data
def alma(data, density=400, v1=-12, v2=12):

    # Generate voltage interval
    xd= np.linspace(v1, v2, density)
    # Generate matrix to store sinograms, each row is a sinogram
    m=np.empty((np.size(data[:,0]),density))
    # Track gaussian data
    param=[]
    files = tqdm(range(np.size(data[:,0])))

    # For each data file
    for i in files:
        # Generate histogram
        n,bins=np.histogram(data[i,:],bins=10,density=True)
        # Fit histogram
        arguments = gausshist(n, bins)
        # Evaluate and append data to matrix
        m[i,:]=gauss(xd,arguments[0],arguments[1],arguments[2])

        plt.hist(data[i,:],bins=10,density=True)
        plt.plot(xd,m[i,:])
        plt.title('id:%d amp:%.3f mu:%.3f sig:%.3f'%(i,arguments[0],arguments[1],arguments[2]))
        plt.savefig('hist/%d.pdf'%i)
        plt.close()

        param.append(arguments)
        files.set_description("Generating sinograms")
    return m,xd, np.array(param)

# Generates sinograms using mu
def generar(mus,density=50,v1=-12,v2=12):
    m=np.zeros((np.size(mus),density))
    y=np.linspace(v1, v2, density)
    for i in range(np.size(mus)):
        m[i,:]= gauss(y, 1, mus[i], 0.707107)
    return m,y

# Generates sinograms using mu and sigma
def gen(mus,sig,density=100,v1=-12,v2=12):
    m=np.zeros((np.size(mus),density))
    y=np.linspace(v1, v2, density)
    for i in range(np.size(mus)):
        r = np.random.normal(mus[i],sig[i],10000)
        n,bins=np.histogram(r,bins=density,density=True)
        arguments = gausshist(n, bins)
        m[i,:]=gauss(y,arguments[0],arguments[1],arguments[2])

        '''
        plt.scatter(np.linspace(0,1,np.size(r)),r,s=1,alpha=0.5)
        plt.ylim(-10,10)
        plt.xlabel(r"$t(au)$")
        plt.ylabel(r"$i_1-i_2(au)$")
        plt.xticks([])
        plt.savefig("samples.png",dpi=300)
        plt.close()
        plt.plot(y[np.size(y)//2:],m[i,np.size(y)//2:],linestyle='--')
        plt.xlabel(r"$i_1-i_2(au)$")
        plt.savefig("hist.png",dpi=300)
        plt.close()
        break
        '''

    return m,y

def fockMatrixCoh(alpha,n,m):
    rho = np.empty((n+1,m+1))
    n = np.arange(n+1)
    m = np.arange(m+1)
    for ni in n:
        rho[ni]=np.abs( np.exp(-1*np.abs(alpha)*np.abs(alpha))*(alpha**ni)*(alpha**m)/np.sqrt(factorial(ni)*factorial(m)) )
    return rho, n, m

def fockMatrixTerm(beta,n,m):
    rho = np.zeros((n+1,m+1))
    n = np.arange(n+1)
    m = np.arange(m+1)
    for ni in n:
        rho[ni][ni]=np.exp(-1*ni*beta)*(1-beta)
    return rho, n, m

def fockMatrixSque(r,t,n,m):
    rho = np.zeros((n+1,m+1))
    n = np.arange(n+1)
    m = np.arange(m+1)
    for ni in n[::2]:
        rho[ni][::2]= (((-1)**ni)*np.sqrt(factorial(2*ni))*np.exp(1j*ni*t)*(np.tanh(r)**ni)/(np.sqrt(np.cosh(r))*(2**ni)*factorial(ni))) * (((-1)**m[::2])*np.sqrt(factorial(2*m[::2]))*np.exp(-1j*m[::2]*t)*(np.tanh(r)**m[::2])/(np.sqrt(np.cosh(r))*(2**m[::2])*factorial(m[::2])))
    return rho, n, m

def generateRhoMatrix():
    rho, n, m = fockMatrixCoh(0,10,10)
    top = rho.ravel()
    X, Y = np.meshgrid(n, m)
    x,y = X.ravel(),Y.ravel()
    bottom = np.zeros_like(top)
    width = depth = 1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    dz = top
    offset = dz + np.abs(dz.min())
    fracs = offset.astype(float)/offset.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.Oranges(norm(fracs.tolist()))
    ax.bar3d(x, y, bottom, width, depth, top, shade=True,color=color_values)
    ax.set_ylabel("n")
    ax.set_xlabel("m")
    ax.set_zlabel(r"$\rho_{nm}$")
    ax.set_title(r"$\rho_{nm}$ for the vacuum state")
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    plt.savefig("rhoVacuum.png",dpi=300)
    plt.close()

    rho, n, m = fockMatrixCoh(3,20,20)
    top = rho.ravel()
    X, Y = np.meshgrid(n, m)
    x,y = X.ravel(),Y.ravel()
    bottom = np.zeros_like(top)
    width = depth = 1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    dz = top
    offset = dz + np.abs(dz.min())
    fracs = offset.astype(float)/offset.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.Oranges(norm(fracs.tolist()))
    ax.bar3d(x, y, bottom, width, depth, top, shade=True,color=color_values)
    ax.set_ylabel("n")
    ax.set_xlabel("m")
    ax.set_zlabel(r"$\rho_{nm}$")
    ax.set_title(r"$\rho_{nm}$ for a coherent state")
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    plt.savefig("rhoCoherent.png",dpi=300)
    plt.close()

    rho, n, m = fockMatrixTerm(0.2,15,15)
    top = rho.ravel()
    X, Y = np.meshgrid(n, m)
    x,y = X.ravel(),Y.ravel()
    bottom = np.zeros_like(top)
    width = depth = 1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    dz = top
    offset = dz + np.abs(dz.min())
    fracs = offset.astype(float)/offset.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.Oranges(norm(fracs.tolist()))
    ax.bar3d(x, y, bottom, width, depth, top, shade=True,color=color_values)
    ax.set_ylabel("n")
    ax.set_xlabel("m")
    ax.set_zlabel(r"$\rho_{nm}$")
    ax.set_title(r"$\rho_{nm}$ for a thermal state")
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    plt.savefig("rhoTerm.png",dpi=300)
    plt.close()

    rho, n, m = fockMatrixSque(2,0,12,12)
    top = rho.ravel()
    X, Y = np.meshgrid(n, m)
    x,y = X.ravel(),Y.ravel()
    bottom = np.zeros_like(top)
    width = depth = 1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    dz = top
    offset = dz + np.abs(dz.min())
    fracs = offset.astype(float)/offset.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.Oranges(norm(fracs.tolist()))
    ax.bar3d(x, y, bottom, width, depth, top, shade=True,color=color_values)
    ax.set_ylabel("n")
    ax.set_xlabel("m")
    ax.set_zlabel(r"$\rho_{nm}$")
    ax.set_title(r"$\rho_{nm}$ for an squeezed vacuum state")
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    plt.savefig("rhoSqueezedVac.png",dpi=300)
    plt.close()

if __name__ == "__main__":

    '''
    # Load data
    fileNames = [ "../tek00"+str(i)+"ALL.csv" if i>=10 else "../tek000"+str(i)+"ALL.csv" for i in range(2,24) ]
    #fileNames = [ "../tek00"+str(i)+"ALL.csv" if i>=10 else "../tek000"+str(i)+"ALL.csv" for i in range(7,21) ]
    data = loadData(fileNames)
    m,x,param=alma(data)

    #n = np.linspace(0, np.size(m[:,0])-1, np.size(m[:,0]))
    n = np.linspace(0,2*np.pi,np.size(m[:,0]))

    print(param)

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x, n)
    ax.set_ylabel(r"$Piezo(\mu m)$")
    ax.set_zlabel(r"$P_{\theta}$")
    ax.plot_surface(X, Y, m, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    plt.show()

    n = np.linspace(0,np.pi,100)
    mus = np.sin(n)
    m,x = generar(mus,v1=-5,v2=5)

    '''

    #generateRhoMatrix()
    # Generate data
    alpha=0
    n = np.linspace(0,2*np.pi,70)
    mus = (np.sin(n)+np.cos(n))*alpha
    #mus = np.zeros_like(n)
    sig = np.zeros_like(n)+(1/np.sqrt(2))
    sig = ((np.sin(2*(n+(1*np.pi/4)))+1)*0.35)+(1/np.sqrt(2))-0.25
    #sig = np.zeros_like(n)+3
    m,x = gen(mus,sig,v1=-8,v2=8)



    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    '''
    X, Y = np.meshgrid(x, n)
    ax.set_ylabel(r"$Piezo(\mu m)$")
    ax.set_zlabel(r"$P_{\theta}$")
    ax.plot_surface(X, Y, m, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    '''
    for i in range(np.size(n)-1,0,-1):
        ax.plot(x, np.zeros_like(x)+n[i], m[i,:],c='mediumvioletred')
    ax.set_ylabel(r"$\phi$")
    ax.set_xlabel(r"$x(v)$")
    ax.set_zlabel(r"$pr(x,\phi)$")
    ax.view_init(elev=0, azim=0)
    plt.savefig("simulatedSqqCoh.png",dpi=300)
    plt.show()
    plt.close()



    # Save data
    with open("rhoSq.dat", "wb") as f:
        f.write(m.tobytes())

    with open("phiSq.dat", "wb") as f:
        f.write(n.tobytes())

    with open("xSq.dat", "wb") as f:
        f.write(x.tobytes())
