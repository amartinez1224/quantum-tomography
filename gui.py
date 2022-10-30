import tomo

import numpy as np
from threading import Thread

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
import matplotlib.animation as animation
from matplotlib import style
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.ticker as ticker

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from PIL import ImageTk, Image

import os

Q, P, W, m, angles, volt = np.empty(0), np.empty(
    0), np.empty(0), np.empty(0), np.empty(0), np.empty(0)
q1C, q2C, p1C, p2C, densityC, kcC = -5, 5, -5, 5, 50, 5
color, angle1, angle2 = 'cividis', 45, 45
t = None
buttonQQ, buttonNM = None, None
fileName, change = "SampleStates/simulated_vacuum.json", False


class counter:
    def __init__(self):
        self.count = 0
        self.start = False
        self.end = False

    def get(self):
        return self.count

    def set(self, x):
        self.count = x

    def started(self):
        self.start = True

    def ended(self):
        self.end = True


class popupWindow(object):
    def __init__(self, master):
        top = self.top = tk.Toplevel(master)
        self.l = tk.Label(top, text="N value")
        self.l.pack()
        self.e = tk.Entry(top, validate='key', validatecommand=(master.register(
            self.validate), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W'))
        self.e.delete(0, tk.END)
        self.e.insert(0, "20")
        self.e.pack()
        self.l1 = tk.Label(top, text="M value")
        self.l1.pack()
        self.f = tk.Entry(top, validate='key', validatecommand=(master.register(
            self.validate), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W'))
        self.f.delete(0, tk.END)
        self.f.insert(0, "20")
        self.f.pack()
        self.b = tk.Button(top, text='Continue', command=self.cleanup)
        self.b.pack()

    def cleanup(self):
        self.n = self.e.get()
        self.m = self.f.get()
        self.top.destroy()
    def validate(self, action, index, value_if_allowed,
                 prior_value, text, validation_type, trigger_type, widget_name):
        if value_if_allowed:
            try:
                int(value_if_allowed)
                return True
            except ValueError:
                return False
        else:
            return False


contador = counter()


def task():
    global contador
    if contador.start:
        progress.set(contador.get()*100/(densityC*densityC))
        if contador.end:
            graficar()
            contador = counter()
            t = None
    win.after(200, task)


def styl(*args):
    global Q, P, W, color, angle1, angle2
    if color != col.get() or angle1 != a1.get() or angle2 != a2.get():
        color = col.get()
        angle1 = a1.get()
        angle2 = a2.get()
        if np.size(Q) < 1 or np.size(P) < 1 or np.size(W) < 1:
            return
        graficar()


def check():
    global q1C, q2C, p1C, p2C, densityC, kcC, t
    q1C, q2C, p1C, p2C, densityC, kcC = qmin.get(), qmax.get(
    ), pmin.get(), pmax.get(), dense.get(), cut.get()
    # b.step(50)
    # win.update()
    for widget in frame1.winfo_children():
        widget.destroy()
    t = Thread(target=data)
    t.start()
    bar()
    '''
    b=bar()
    while(t.is_alive()):
        b.wait_variable(progress)
        b.update()
    graficar()
    '''


def data():
    global Q, P, W, m, angles, volt, contador, fileName, change
    if np.size(m) < 1 or np.size(angles) < 1 or np.size(volt) < 1 or change:
        try:
            m, angles, volt = tomo.loadData(fileName)
        except:
            tk.messagebox.showinfo(
                'Error', f'The file {txtFile.get()} is not formatted as expected!')
            return
        change = False
    generateButton["state"] = 'disabled'
    contador.started()
    Q, P, W = tomo.tomo(m, angles, volt, contador, q1=q1C,
                        q2=q2C, p1=p1C, p2=p2C, density=densityC, kc=kcC)
    contador.ended()
    generateButton["state"] = 'normal'


def changeData():
    file = tk.filedialog.askopenfilename(initialdir=os.path.abspath(
        ""), title="Load data", filetypes=(("Json files", "*.json"), ("all files", "*.*")))
    if file == "":
        tk.messagebox.showinfo('Error', 'The file does not exist!')
    global fileName, change
    fileName = file
    change = True
    txtFile.set(fileName.split("/")[-1])


def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def densityMatrixQQ():
    global W, Q, P
    rhopP, listPosP, listNegP, rhopQ, listPosQ, listNegQ = tomo.quadratureToRho(
        W, Q, P)
    '''
    fig = plt.figure(num='Density matrix')
    ax = fig.add_subplot(121, projection='3d')
    ax.scatter(listPosP, listNegP,rhopP,c=rhopP,cmap='gray')
    ax.set_title("P")
    ax.view_init(elev=angle1, azim=angle2)
    ax = fig.add_subplot(122, projection='3d')
    ax.scatter(listPosQ, listNegQ,rhopQ,c=rhopQ,cmap='gray')
    ax.set_title("Q")
    ax.view_init(elev=angle1, azim=angle2)
    plt.show()
    '''
    Xp, Yp, Zp = tomo.rhoFitting(
        rhopP, listPosP, listNegP, np.max(P), np.min(P))
    Xq, Yq, Zq = tomo.rhoFitting(
        rhopQ, listPosQ, listNegQ, np.max(Q), np.min(Q))
    #Xp, Yp, Zp = tomo.rhoInterpolate(rhopP, listPosP, listNegP,np.max(P),np.min(P))
    #Xq, Yq, Zq = tomo.rhoInterpolate(rhopQ, listPosQ, listNegQ,np.max(Q),np.min(Q))

    plotDensityMatrixQQ(Xp, Yp, Zp, Xq, Yq, Zq)


state = 0


def changeState(dir):
    global state
    if dir:
        state += 1
    else:
        state -= 1
    if state == 0:
        generateButton["state"] = "normal"
    else:
        generateButton["state"] = "disabled"


def argumentsRhoNM():
    w = popupWindow(win)
    changeState(True)
    buttonNM["state"] = "disabled"
    win.wait_window(w.top)
    ni, mi = int(w.n), int(w.m)
    densityMatrixNM(ni, mi)
    buttonNM["state"] = "normal"
    changeState(False)


def argumentsMatrixQQ():
    changeState(True)
    buttonQQ["state"] = "disabled"
    densityMatrixQQ()
    buttonQQ["state"] = "normal"
    changeState(False)


def densityMatrixNM(ni, mi):
    global W, Q, P
    rhopP, listPosP, listNegP, rhopQ, listPosQ, listNegQ = tomo.quadratureToRho(
        W, Q, P)
    Xp, Yp, Zp = tomo.rhoFitting(
        rhopP, listPosP, listNegP, np.max(P), np.min(P))
    Xq, Yq, Zq = tomo.rhoFitting(
        rhopQ, listPosQ, listNegQ, np.max(Q), np.min(Q))
    #Xp, Yp, Zp = tomo.rhoInterpolate(rhopP, listPosP, listNegP,np.max(P),np.min(P))
    #Xq, Yq, Zq = tomo.rhoInterpolate(rhopQ, listPosQ, listNegQ,np.max(Q),np.min(Q))
    rhon, nr, mr = tomo.rhoFock(Zq, Xq, Yq, ni, mi)
    plotDensityMatrixNM(rhon, nr, mr)


def plotDensityMatrixQQ(Xp, Yp, Zp, Xq, Yq, Zq):
    fig, axes = plt.subplots(figsize=(
        10, 4), nrows=1, ncols=2, num="Density matrix in quadrature representation")
    l = np.array([1/100000., 1/10000., 1/1000., 1/100., 1/10., 1])
    ax = axes[0]
    h = ax.contour(Xp, Yp, Zp, levels=l, norm=colors.LogNorm(
        vmin=1/1000000., vmax=1), cmap='Blues')
    ax.set_xlabel(r'$p$')
    ax.set_ylabel(r'$p^\prime$')
    ax.set_aspect('equal')
    ax = axes[1]
    ax.contour(Xq, Yq, Zq, levels=l, norm=colors.LogNorm(
        vmin=1/1000000., vmax=1), cmap='Blues')
    ax.set_xlabel(r'$q$')
    ax.set_ylabel(r'$q^\prime$')
    ax.set_aspect('equal')
    fig.colorbar(h, ax=axes.ravel().tolist(), format=ticker.FuncFormatter(fmt))
    # plt.savefig("rhoQPVac.png",dpi=300)
    plt.show()


def plotDensityMatrixNM(rhon, nr, mr):
    '''
    rhon, nr, mr = tomo.rhoFock(rhopP, listPosP, listNegP,np.max(P),np.min(P))
    X, Y = np.meshgrid(nr, mr)
    x,y = X.ravel(),Y.ravel()
    top = rhon.ravel()
    top = np.abs(top)
    bottom = np.zeros_like(top)
    width = depth = 1
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(121, projection='3d')
    dz = top
    offset = dz + np.abs(dz.min())
    fracs = offset.astype(float)/offset.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.Greens(norm(fracs.tolist()))
    ax.bar3d(x, y, bottom, width, depth, top, shade=True,color=color_values)
    ax.set_ylabel("n")
    ax.set_xlabel("m")
    ax.set_zlabel(r"$|\rho_{nm}|$")
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    '''

    fig = plt.figure(num="Density matrix in Fock representation")
    X, Y = np.meshgrid(mr, nr)
    x, y = X.ravel(), Y.ravel()
    top = rhon.ravel()
    top = np.abs(top)
    bottom = np.zeros_like(top)
    width = depth = 1
    ax = fig.add_subplot(111, projection='3d')
    dz = top
    offset = dz + np.abs(dz.min())
    fracs = offset.astype(float)/offset.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.Greens(norm(fracs.tolist()))
    ax.bar3d(x, y, bottom, width, depth, top, shade=True, color=color_values)
    ax.set_ylabel("n")
    ax.set_xlabel("m")
    ax.set_zlabel(r"$|\rho_{nm}|$")
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    # plt.savefig("rhoNMVac.png",dpi=300)
    plt.show()


def graficar():
    for widget in frame1.winfo_children():
        widget.destroy()
    b = bar(indetermine=True)
    global Q, P, W, color, angle1, angle2, buttonQQ, buttonNM
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(P, Q)
    ax.set_ylabel(r"$q$")
    ax.set_xlabel(r"$p$")
    ax.set_zlabel(r"$W(q,p)$")
    ax.view_init(elev=angle1, azim=angle2)
    h = ax.plot_surface(X, Y, W, rstride=1, cstride=1,
                        cmap=color, edgecolor='none')
    #ax.contour(X, Y, W)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    # plt.savefig("number",dpi=300)
    plt.close()
    Plot = FigureCanvasTkAgg(fig, master=frame1)
    Plot.draw()
    Plot.get_tk_widget().pack(side='top', fill='both', expand=1)
    b.destroy()
    toolbar = NavigationToolbar2Tk(Plot, frame1)
    toolbar.children['!button2'].pack_forget()
    toolbar.children['!button3'].pack_forget()
    buttonQQ = tk.Button(
        master=toolbar, text="Density Matrix quadrature", command=lambda: argumentsMatrixQQ())
    buttonNM = tk.Button(
        master=toolbar, text="Density Matrix Fock", command=lambda: argumentsRhoNM())
    buttonNM.pack(side="right")
    buttonQQ.pack(side="right")
    toolbar.update()


def bar(indetermine=False):
    progress.set(0)
    s = ttk.Style()
    s.theme_use('clam')
    TROUGH_COLOR = 'white'
    BAR_COLOR = '#308fac'
    s.configure("bar.Horizontal.TProgressbar", troughcolor=TROUGH_COLOR,
                bordercolor=TROUGH_COLOR, background=BAR_COLOR, lightcolor=BAR_COLOR, darkcolor=BAR_COLOR)
    if indetermine:
        loadingBar = ttk.Progressbar(
            frame1, mode="indeterminate", style="bar.Horizontal.TProgressbar")
        loadingBar.start()
    else:
        loadingBar = ttk.Progressbar(
            frame1, style="bar.Horizontal.TProgressbar", variable=progress)
        # loadingBar.config(maximum=densityC*densityC)
    loadingBar.place(relx=.5, rely=.5, anchor="center", relwidth=0.5)

    return loadingBar


if __name__ == "__main__":

    win = tk.Tk()
    win.geometry('1125x900')
    win.title('Quantum Tomography')
    win.config(cursor="arrow")

    frame1 = tk.Frame(win, background='#308fac')
    frame1.place(relwidth=0.8, relheight=1)

    frame2 = tk.Frame(win)
    frame2.place(relx=0.8, relwidth=0.2, relheight=1)

    tk.Label(master=frame2).grid(row=0, column=0)
    tk.Label(master=frame2).grid(row=0, column=1)

    tk.Label(master=frame2, text="Q min:").grid(row=1, column=0)
    tk.Label(master=frame2, text="Q max:").grid(row=2, column=0)
    tk.Label(master=frame2, text="P min:").grid(row=3, column=0)
    tk.Label(master=frame2, text="P max:").grid(row=4, column=0)
    tk.Label(master=frame2, text="Density:").grid(row=5, column=0)
    tk.Label(master=frame2, text="Kc:").grid(row=6, column=0)
    generateButton = tk.Button(frame2, text='Tomography', command=check)
    generateButton.grid(row=7, column=1)

    tk.Label(master=frame2).grid(row=8, column=0)
    tk.Label(master=frame2).grid(row=8, column=1)

    tk.Label(master=frame2, text="Color:").grid(row=9, column=0)
    tk.Label(master=frame2, text="Angle 1:").grid(row=10, column=0)
    tk.Label(master=frame2, text="Angle 2:").grid(row=11, column=0)
    tk.Label(master=frame2).grid(row=12, column=0)
    tk.Label(master=frame2).grid(row=12, column=1)

    tk.Label(master=frame2, text="Data:").grid(row=15, column=0)

    txtFile = tk.StringVar()
    bFile = tk.Button(frame2, textvariable=txtFile, command=changeData)
    txtFile.set(fileName.split("/")[-1])
    bFile.grid(row=15, column=1)

    tk.Label(master=frame2).grid(row=16, column=0)
    tk.Label(master=frame2).grid(row=16, column=1)

    qmin = tk.DoubleVar()
    qmin.set(q1C)
    tk.Spinbox(frame2, from_=-1000000, to=1000000,
               textvariable=qmin, width=8).grid(row=1, column=1)
    qmax = tk.DoubleVar()
    qmax.set(q2C)
    tk.Spinbox(frame2, from_=-1000000, to=1000000,
               textvariable=qmax, width=8).grid(row=2, column=1)
    pmin = tk.DoubleVar()
    pmin.set(p1C)
    tk.Spinbox(frame2, from_=-1000000, to=1000000,
               textvariable=pmin, width=8).grid(row=3, column=1)
    pmax = tk.DoubleVar()
    pmax.set(p2C)
    tk.Spinbox(frame2, from_=-1000000, to=1000000,
               textvariable=pmax, width=8).grid(row=4, column=1)
    dense = tk.IntVar()
    dense.set(densityC)
    tk.Spinbox(frame2, from_=-1000000, to=1000000,
               textvariable=dense, width=8).grid(row=5, column=1)
    cut = tk.DoubleVar()
    cut.set(kcC)
    tk.Spinbox(frame2, from_=-100, to=100, textvariable=cut,
               width=8).grid(row=6, column=1)

    col = tk.StringVar()
    op = ['viridis', 'plasma', 'inferno', 'magma', 'cividis',
          'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds']
    col.set(color)
    menu = tk.OptionMenu(frame2, col, *op, command=styl)
    menu.config(width=8)
    menu.grid(row=9, column=1)

    a1 = tk.DoubleVar()
    a1.set(angle1)
    e = tk.Entry(frame2, textvariable=a1, width=8)
    e.grid(row=10, column=1)
    e.bind("<FocusOut>", styl)
    e.bind("<Return>", styl)

    a2 = tk.DoubleVar()
    a2.set(angle2)
    e = tk.Entry(frame2, textvariable=a2, width=8)
    e.grid(row=11, column=1)
    e.bind("<FocusOut>", styl)
    e.bind("<Return>", styl)

    frame2.grid_rowconfigure(8, weight=1)
    frame2.grid_rowconfigure(12, weight=1)
    frame2.grid_columnconfigure(0, weight=1)
    frame2.grid_columnconfigure(1, weight=1)

    progress = tk.IntVar()
    tk.Label(master=frame2, text="Progress %").grid(row=17, column=0)
    pr = tk.Entry(frame2, textvariable=progress, width=8)
    pr.grid(row=17, column=1)
    pr.config(state=tk.DISABLED)

    win.after(1000, task)
    win.mainloop()
