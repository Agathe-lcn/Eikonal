import os
from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

GRIDS = [100, 200, 500, 1000, 2000, 5000, 10000]

h0 = 1 / 100

line_angle = -20
rect_length = 0.5
rect_height = 0.15
rect_shift = 0.1
contact_length = rect_length - rect_shift

circle_radius = 0.28
ring1_width = 0.10
ring2_width = 0.10

cmap = ListedColormap(["white", "#bfe3ff", "#ffbfbf"])

# Mettre a False pour desactiver completement toutes les sorties image.
GENERATE_IMAGES = True

# Modes d'execution disponibles pour la generation des fichiers de configuration.
# nproc = 1 correspond au mode sequentiel.
# max_depth = -1 desactive la limitation de profondeur.
EXECUTION_MODES = [
    {
        "name": "seq",
        "nproc": 1,
        "overlap": None,
        "max_depth": -1,
    },
    {
        "name": "mpi_np2_ov5_md10",
        "nproc": 2,
        "overlap": 5,
        "max_depth": 10,
    },
]

DATASETS_ROOT = Path("datasets")
IMAGES_ROOT = Path("images")
IMAGE_CACHE_ROOT = IMAGES_ROOT / "_cache"

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

def validate_execution_mode(mode):
    nproc = mode["nproc"]
    overlap = mode["overlap"]
    max_depth = mode["max_depth"]

    if nproc < 1:
        raise ValueError(f"Mode invalide {mode['name']}: nproc doit etre >= 1")
    if nproc == 1 and overlap not in (None, 0):
        raise ValueError(f"Mode invalide {mode['name']}: overlap doit etre absent en sequentiel")
    if nproc > 1 and (overlap is None or overlap < 1):
        raise ValueError(f"Mode invalide {mode['name']}: overlap doit etre >= 1 en MPI")
    if max_depth < -1:
        raise ValueError(f"Mode invalide {mode['name']}: max_depth doit etre >= 0 ou -1")


def mode_directory_name(mode):
    if mode.get("name"):
        return mode["name"]

    if mode["nproc"] == 1:
        suffix = "full" if mode["max_depth"] < 0 else f"md{mode['max_depth']}"
        return f"seq_{suffix}"

    depth_suffix = "full" if mode["max_depth"] < 0 else f"md{mode['max_depth']}"
    return f"mpi_np{mode['nproc']}_ov{mode['overlap']}_{depth_suffix}"


def mode_description(mode):
    depth_label = "desactive" if mode["max_depth"] < 0 else str(mode["max_depth"])
    if mode["nproc"] == 1:
        return f"mode=sequentiel, nproc=1, max_depth={depth_label}"
    return f"mode=mpi, nproc={mode['nproc']}, overlap={mode['overlap']}, max_depth={depth_label}"


def save_dataset(path, n, h, sources, scenario, params, mode):
    with open(path, "w") as f:
        f.write(f"n = {n}\n")
        f.write(f"m = {n}\n")
        f.write(f"h = {h:.10f}\n\n")
        if mode["nproc"] > 1:
            f.write(f"overlap = {mode['overlap']}\n")
        if mode["max_depth"] >= 0:
            f.write(f"max_depth = {mode['max_depth']}\n")
        if mode["nproc"] > 1 or mode["max_depth"] >= 0:
            f.write("\n")
        f.write(f"# Scenario : {scenario}\n")
        f.write(f"# Parametres : {params}\n")
        f.write(f"# Execution : {mode_description(mode)}\n")
        f.write(f"# Nombre de sources : {len(sources)}\n\n")
        f.write("sources:\n")
        for x, y in sources:
            f.write(f"{x} {y}\n")
        f.write("\nwalls:\n")

def plot(mat, sources, out):
    plt.figure(figsize=(6, 6))
    plt.imshow(mat, origin="lower", cmap=cmap, vmin=0, vmax=2, interpolation="nearest")
    plt.contour(mat, levels=[0.5, 1.5], colors="black", linewidths=0.5)
    if sources:
        xs, ys = zip(*sources)
        plt.scatter(xs, ys, s=10, c="red", edgecolors="black", linewidths=0.2)
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()

def ensure_mode_directories(mode_dir):
    for root in (DATASETS_ROOT, IMAGES_ROOT):
        for scenario_name in ("line", "circle"):
            os.makedirs(root / mode_dir / scenario_name, exist_ok=True)


def ensure_cache_directories():
    for scenario_name in ("line", "circle"):
        os.makedirs(IMAGE_CACHE_ROOT / scenario_name, exist_ok=True)


def scenario_payloads(n):
    return {
        "line": {
            "sources": line_sources(n, line_angle, contact_length),
            "materials": line_materials(n),
            "scenario": "ligne",
            "params": f"angle={line_angle}, longueur={rect_length}, hauteur={rect_height}, decalage={rect_shift}",
            "basename": f"line_{n}",
        },
        "circle": {
            "sources": circle_sources(n, circle_radius),
            "materials": circle_materials(n),
            "scenario": "anneaux",
            "params": f"rayon={circle_radius}, ep1={ring1_width}, ep2={ring2_width}",
            "basename": f"circle_{n}",
        },
    }


def build_image_cache(n):
    ensure_cache_directories()
    payloads = scenario_payloads(n)

    for scenario_name, payload in payloads.items():
        cache_path = IMAGE_CACHE_ROOT / scenario_name / f"{payload['basename']}.png"
        if not cache_path.exists():
            plot(payload["materials"], payload["sources"], cache_path)

    return payloads


def generate_for_mode(mode, write_images):
    validate_execution_mode(mode)
    mode_dir = mode_directory_name(mode)
    ensure_mode_directories(mode_dir)

    for n in GRIDS:
        h = h0 * 100 / n
        payloads = build_image_cache(n)

        line_path = DATASETS_ROOT / mode_dir / "line" / f"line_{n}.txt"
        line_img_path = IMAGES_ROOT / mode_dir / "line" / f"line_{n}.png"
        line_payload = payloads["line"]
        save_dataset(
            line_path,
            n,
            h,
            line_payload["sources"],
            line_payload["scenario"],
            line_payload["params"],
            mode,
        )
        if write_images:
            shutil.copyfile(IMAGE_CACHE_ROOT / "line" / f"{line_payload['basename']}.png", line_img_path)

        circle_path = DATASETS_ROOT / mode_dir / "circle" / f"circle_{n}.txt"
        circle_img_path = IMAGES_ROOT / mode_dir / "circle" / f"circle_{n}.png"
        circle_payload = payloads["circle"]
        save_dataset(
            circle_path,
            n,
            h,
            circle_payload["sources"],
            circle_payload["scenario"],
            circle_payload["params"],
            mode,
        )
        if write_images:
            shutil.copyfile(IMAGE_CACHE_ROOT / "circle" / f"{circle_payload['basename']}.png", circle_img_path)


def main():
    for index, execution_mode in enumerate(EXECUTION_MODES):
        generate_for_mode(execution_mode, write_images=(GENERATE_IMAGES and index == 0))
        print(f"generated {mode_directory_name(execution_mode)}")

    print("OK")


if __name__ == "__main__":
    main()
