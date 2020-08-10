import numpy as np
from threading import Thread, Lock
from queue import Queue
from scipy import interpolate
from scipy.special import gamma, factorial, hyp1f1, eval_hermite
from scipy.optimize import curve_fit

W = None
rhonm = None
progressVar = None
count = 0


def gauss2d(t, amp, muX, muY, sigX, sigY, theta):
    x,y = t
    a = (np.cos(theta)**2)/(2*sigX**2) + (np.sin(theta)**2)/(2*sigY**2)
    b = -(np.sin(2*theta))/(4*sigX**2) + (np.sin(2*theta))/(4*sigY**2)
    c = (np.sin(theta)**2)/(2*sigX**2) + (np.cos(theta)**2)/(2*sigY**2)
    f = amp*np.exp( - (a*((x-muX)**2) + 2*b*(x-muX)*(y-muY) + c*((y-muY)**2)))
    return f.ravel()

def loadData(mf,nf,xf):
    m = None
    n = None
    x = None
    with open(xf, "rb") as f:
        x = np.frombuffer(f.read(),'float64')
    with open(nf, "rb") as f:
        n = np.frombuffer(f.read(),'float64')
    with open(mf, "rb") as f:
        m = (np.frombuffer(f.read(),'float64')).reshape((np.size(n),np.size(x)))
    return m,n,x

def wigner(iq,ip,q,p,m,angles,volt,kc):
    int=0
    global W
    for angle in range(np.size(angles)):
        convolution = np.sum( m[angle,:] * Kcomp(q,p,angles[angle],volt,kc) )
        int = int + convolution
    W[iq,ip] = int*np.abs(angles[1]-angles[0])*np.abs(volt[1]-volt[0])/(2*np.pi*np.pi)

def K(arg,kc):
    return ((kc**2)/2.)*(1-( (kc**2)*(arg**2)/4. )+( (kc**4)*(arg**4)/72. )-( (kc**6)*(arg**6)/2880. )+( (kc**8)*(arg**8)/201600. ))

def Kor(arg,kc):
    return (np.cos(kc*arg) + kc*arg*np.sin(kc*arg) - 1)/(arg**2)

def Kcomp(q,p,angle,volt,kc):
    turn = 0.01
    arg = ( q*np.cos(angle) ) + ( p*np.sin(angle) ) - volt
    arg[np.abs(arg*kc)<turn] = K(arg[np.abs(arg*kc)<turn],kc)
    arg[np.abs(arg*kc)>=turn] = Kor(arg[np.abs(arg*kc)>=turn],kc)
    return arg

que = Queue()
queRho = Queue()
lock = Lock()

def upadateBar():
    global progressVar,count, norm
    with lock:
        count+=1
        progressVar.set(count)

def worker():
    while True:
        item = que.get()
        wigner(*item)
        upadateBar()
        que.task_done()

def workerRho():
    while True:
        item = queRho.get()
        quadratureToFock(*item)
        queRho.task_done()

for i in range(4):
     t = Thread(target=worker)
     t.daemon = True
     t.start()
     t = Thread(target=workerRho)
     t.daemon = True
     t.start()

def quadratureToRho(w,q,p):
    listPosP = np.empty(np.size(p)*np.size(p))
    listNegP = np.empty(np.size(p)*np.size(p))
    resultP = np.empty(np.size(p)*np.size(p))
    listPosQ = np.empty(np.size(q)*np.size(q))
    listNegQ = np.empty(np.size(q)*np.size(q))
    resultQ = np.empty(np.size(q)*np.size(q))
    k=0
    for i in range(np.size(p)):
        for j in range(np.size(p)):
            listPosP[k] = p[i]+p[j]
            listNegP[k] = p[i]-p[j]
            int = w[:,i]*np.exp(2j*q*p[j])
            resultP[k] = np.abs(np.sum((int[1:]+int[:-1])*np.abs(q[1:]-q[:-1])/2))
            listPosQ[k] = q[i]+q[j]
            listNegQ[k] = q[i]-q[j]
            int = w[i,:]*np.exp(2j*p*q[j])
            resultQ[k] = np.abs(np.sum((int[1:]+int[:-1])*np.abs(p[1:]-p[:-1])/2))
            k+=1
    return resultP, listPosP, listNegP, resultQ, listPosQ, listNegQ

def quadratureToFock(n,m,rho,x,xp):
    integral = []
    global rhoNM
    for i in range(np.size(x)):
        int = rho[i,:]*np.exp(-0.5*((x[i]*x[i])+(xp*xp)))*eval_hermite(n,x[i])*eval_hermite(m,xp)
        integral.append(np.sum(np.abs(xp[1:]-xp[:-1])*(int[1:]+int[:-1])/2))
    integral = np.array(integral)
    integral = np.sum(np.abs(x[1:]-x[:-1])*(integral[1:]+integral[:-1])/2)
    rhoNM[n][m] = integral/(np.sqrt(np.pi*(2**m)*(2**n)*factorial(n)*factorial(m)))

# Interpolate rho qq usig cubic splines
def rhoInterpolate(rho,q,qp,qmax,qmin,density=100):
    # Get data in limits
    index = np.logical_and(np.logical_and(q>qmin,q<qmax),np.logical_and(qp>qmin,qp<qmax))
    rho,q,qp=rho[index],q[index],qp[index]
    # Perform interpolation
    f = interpolate.interp2d(q, qp, rho, kind='quintic')
    # Genrate new space
    x,y = np.linspace(qmin,qmax,density),np.linspace(qmin,qmax,density)
    # Obtain values of rho in new space
    rho = f(x,y)
    # Return new space and rho
    return x,y,rho

def rhoFitting(rho,q,qp,qmax,qmin,density=100):
    # Get data in limits
    index = np.logical_and(np.logical_and(q>qmin,q<qmax),np.logical_and(qp>qmin,qp<qmax))
    rho,q,qp=rho[index],q[index],qp[index]
    # Perform fitting
    arguments, variance = curve_fit(gauss2d, (q,qp), rho, p0=[1, q[np.argmax(rho)], qp[np.argmax(rho)], 0.7, 0.7, 0], bounds=([0,-100,-100,0.05,0.05,0],[10,100,100,100,100,2*np.pi]))
    # Genrate new space
    x,y = np.linspace(qmin,qmax,density),np.linspace(qmin,qmax,density)
    X,Y = np.meshgrid(x,y)
    # Obtain values of rho in new space
    rho = gauss2d((X,Y),arguments[0],arguments[1],arguments[2],arguments[3],arguments[4],arguments[5]).reshape(density, density)
    # Return new space and rho
    return x,y,rho

def rhoFock(rho,x,xp,n=20,m=20):
    global rhoNM
    rhoNM = np.empty((n,m))
    n = np.arange(n)
    m = np.arange(m)
    for i in n:
        for j in m:
            queRho.put((i,j,rho,x,xp))
    queRho.join()
    return rhoNM, n, m

def tomo(m,angles,volt,prog,q1=-1,q2=1,p1=-1,p2=1,density=100,kc=2):
    global W, progressVar, count
    Q = np.linspace(q1,q2,density)
    P = np.linspace(p1,p2,density)
    W = np.zeros((np.size(Q),np.size(P)))
    progressVar = prog
    count = 0
    for q in range(density):
        for p in range(density):
            que.put((q,p,Q[q],P[p],m,angles,volt,kc))
    que.join()
    return Q,P,W

def comb(n,k):
    return factorial(n)/(factorial(k)*factorial(n-k))

def term1(n,m,d,phi,j1,j2):
    l1, l2 = np.meshgrid(np.arange(n-j1+((d-j2)/2)),np.arange(j1+(j2/2)))
    x = -1/8.
    up = (np.cos(phi)**(2*(n-j1-l1)-j2+d))+(np.sin(phi)**(2*(j1-l2)+j2))
    down = factorial(l1)*factorial(l2)*factorial(2*(n-j1-l1)-j2+d)*factorial(2*(j1-l2)+j2)
    return np.sum((up/down)*term2(n,m,d,l1,l2))

def term2(n,m,d,l1,l2):
    kron = m-l1-l2
    kron = (kron != 0).astype(int)
    return gamma(n+0.5*(d+(d%2))-l1-l2+1)*kron

def R(n,m,d,phi):
    sum = 0
    for j1 in np.arange(n+1):
        for j2 in np.arange(d+1):
            co = comb(n,j1)*comb(d,j2)
            if co != 0:
                sum += (-1**j2)*co*factorial(2*(n-j1)+d-j2)*factorial(2*j1+j2)*term1(n,m,d,phi,j1,j2)
    return sum*( ((2**(3/2)*1j)**(d%2)) * (2**(n+(d/2)+1)) * ((-1j)**(2*n+d)) )

def rhoElement(n,d,data,angles,volt):
    sum = 0
    x=volt/np.max(volt)
    for m in range(int(n+0.5*d)):
        for f in range(np.size(angles)):
            sum+=R(n,m,d,angles[f])*np.mean(np.exp(-2*x*x)*(x%2)*hyp1f1(m-n-0.5*(d-((d+1)%2)),0.5+(d%2),2*(x**2)))
    return sum/np.size(angles)

def rho(n_i,m_i,matriz,angles,volt):
    rho = np.empty((n_i,m_i),dtype="complex64")
    for n in range(n_i):
        for m in range(m_i):
            print(n,m)
            d = m-n
            rho[n,m] = rhoElement(n,d,matriz,angles,volt)
    return rho

n = np.linspace(-5,5,100)
def gauss(x, amp, mu, sig):
    return amp*np.exp(-(x-mu)**2/(2.*sig**2))
