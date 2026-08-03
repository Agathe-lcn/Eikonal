# dataset_generator.py
# Version simplifiée prête à l'emploi

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

GRIDS=[100,200,500,1000,2000,5000,10000]

h0=1/100

line_angle=-20
rect_length=0.5
rect_height=0.15
rect_shift=0.1
contact_length= rect_length - rect_shift

circle_radius=0.28
ring1_width=0.10
ring2_width=0.10

cmap=ListedColormap(["white","#bfe3ff","#ffbfbf"])

def line_sources(n,angle_deg,length):
    a=np.deg2rad(angle_deg)
    c=np.array([n/2,n/2])
    d=np.array([np.cos(a),np.sin(a)])
    L=length*n
    pts=[]
    for s in np.linspace(-L/2,L/2,int(8*L)+1):
        p=np.round(c+s*d).astype(int)
        x,y=p
        if 0<=x<n and 0<=y<n: pts.append((x,y))
    return sorted(set(pts))

def circle_sources(n,r):
    R=r*n
    c=np.array([n/2,n/2])
    pts=[]
    for t in np.linspace(0,2*np.pi,max(100,int(12*R))):
        p=np.round(c+[R*np.cos(t),R*np.sin(t)]).astype(int)
        x,y=p
        if 0<=x<n and 0<=y<n: pts.append((x,y))
    return sorted(set(pts))

def line_materials(n):
    X,Y=np.meshgrid(np.arange(n),np.arange(n))
    a=np.deg2rad(line_angle)
    cx=cy=n/2
    xt=(X-cx)*np.cos(a)+(Y-cy)*np.sin(a)
    yn=-(X-cx)*np.sin(a)+(Y-cy)*np.cos(a)
    L=rect_length*n
    H=rect_height*n
    shift=rect_shift*n/2
    m=np.zeros((n,n),int)
    r1=(yn>=0)&(yn<=H)&(xt>=-L/2+shift)&(xt<=L/2+shift)
    r2=(yn<=0)&(yn>=-H)&(xt>=-L/2-shift)&(xt<=L/2-shift)
    m[r1]=1
    m[r2]=2
    return m

def circle_materials(n):
    X,Y=np.meshgrid(np.arange(n),np.arange(n))
    R=np.sqrt((X-n/2)**2+(Y-n/2)**2)
    r=circle_radius*n
    w1=ring1_width*n
    w2=ring2_width*n
    m=np.zeros((n,n),int)
    m[(R>=r-w1)&(R<=r)]=1
    m[(R>=r)&(R<=r+w2)]=2
    return m

def save_dataset(path,n,h,sources,scenario,params):
    with open(path,"w") as f:
        f.write(f"n = {n}\n")
        f.write(f"m = {n}\n")
        f.write(f"h = {h:.10f}\n\n")
        f.write(f"# Scenario : {scenario}\n")
        f.write(f"# Parametres : {params}\n")
        f.write(f"# Nombre de sources : {len(sources)}\n\n")
        f.write("sources:\n")
        for x,y in sources:
            f.write(f"{x} {y}\n")

def plot(mat,sources,out):
    plt.figure(figsize=(6,6))
    plt.imshow(mat,origin="lower",cmap=cmap,vmin=0,vmax=2,interpolation="nearest")
    plt.contour(mat,levels=[0.5,1.5],colors="black",linewidths=0.5)
    if sources:
        xs,ys=zip(*sources)
        plt.scatter(xs,ys,s=10,c="red",edgecolors="black",linewidths=0.2)
    plt.axis("equal"); plt.tight_layout(); plt.savefig(out,dpi=300); plt.close()

for d in ["datasets/line","datasets/circle","images/line","images/circle"]:
    os.makedirs(d,exist_ok=True)

for n in GRIDS:
    h=h0*100/n
    s=line_sources(n,line_angle,contact_length)
    save_dataset(f"datasets/line/line_{n}.txt",n,h,s,"ligne",
                 f"angle={line_angle}, longueur={rect_length}, hauteur={rect_height}, decalage={rect_shift}")
    plot(line_materials(n),s,f"images/line/line_{n}.png")
    s=circle_sources(n,circle_radius)
    save_dataset(f"datasets/circle/circle_{n}.txt",n,h,s,"anneaux",
                 f"rayon={circle_radius}, ep1={ring1_width}, ep2={ring2_width}")
    plot(circle_materials(n),s,f"images/circle/circle_{n}.png")
print("OK")
