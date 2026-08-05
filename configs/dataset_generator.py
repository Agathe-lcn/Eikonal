import os
from pathlib import Path
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

GRIDS = [100, 200, 500] #[100, 200, 500, 1000, 2000, 5000, 10000]

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
GENERATE_IMAGES = False

# Mettre a True pour generer aussi les images de resultat des solveurs dans chaque dossier de sortie.
GENERATE_RESULT_IMAGES = True

# Mettre a True pour executer les solveurs sur des fichiers de configuration existants.
RUN_SOLVERS = True

# Liste de motifs relatifs au dossier configs/ pour choisir un ou plusieurs fichiers .txt a executer.
# Exemple: ["datasets/seq/line/line_100.txt", "datasets/mpi_np2_ov5_md10/circle/*.txt"]
# Une liste vide signifie: tous les .txt sous datasets/.
RUN_CONFIG_PATTERNS = []

# Mettre a True pour supprimer le dossier de sortie d'un fichier avant de relancer la simulation.
CLEAN_OUTPUT_DIRECTORIES = False

# Modes d'execution disponibles pour la generation des fichiers de configuration.
# nproc = 1 correspond au mode sequentiel.
# max_depth = -1 desactive la limitation de profondeur.
EXECUTION_MODES = [
    {
        "name": "seq",
        "nproc": 1,
        "overlap": None,
        "max_depth": 10,
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
SCRIPT_ROOT = Path(__file__).resolve().parent
SEQUENTIAL_EXECUTABLE = (SCRIPT_ROOT.parent / "Eikonal2D" / "bin" / "generate_grid.exe").resolve()
MPI_EXECUTABLE = (SCRIPT_ROOT.parent / "Eikonal2D_MPI" / "bin" / "generate_grid_mpi.exe").resolve()
MPIEXEC_EXECUTABLE = Path(r"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe")
PYTHON_EXECUTABLE = Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "Python" / "Python312" / "python.exe"
SEQUENTIAL_VISUALIZATION_SCRIPT = (SCRIPT_ROOT.parent / "Eikonal2D" / "visualization" / "visu_fim.py").resolve()
MPI_GLOBAL_VISUALIZATION_SCRIPT = (SCRIPT_ROOT.parent / "Eikonal2D_MPI" / "visualization" / "visu_fim_global.py").resolve()
MPI_LOCAL_VISUALIZATION_SCRIPT = (SCRIPT_ROOT.parent / "Eikonal2D_MPI" / "visualization" / "visu_fim_local.py").resolve()

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


def dataset_mode_map():
    mode_map = {}
    for mode in EXECUTION_MODES:
        validate_execution_mode(mode)
        mode_map[mode_directory_name(mode)] = mode
    return mode_map


def resolve_config_files():
    patterns = RUN_CONFIG_PATTERNS or ["datasets/**/*.txt"]
    config_paths = []

    for pattern in patterns:
        config_paths.extend(SCRIPT_ROOT.glob(pattern))

    config_files = sorted({path.resolve() for path in config_paths if path.is_file() and path.suffix == ".txt"})
    return [path for path in config_files if DATASETS_ROOT.name in path.parts]


def execution_mode_for_config(config_path):
    dataset_root = (SCRIPT_ROOT / DATASETS_ROOT).resolve()
    relative_path = config_path.resolve().relative_to(dataset_root)
    parts = relative_path.parts
    if len(parts) < 2:
        raise ValueError(f"Configuration hors structure attendue: {config_path}")

    mode = dataset_mode_map().get(parts[0])
    if mode is None:
        raise ValueError(f"Impossible d'associer {config_path} a un mode declare dans EXECUTION_MODES")
    return mode


def output_directory_for_config(config_path):
    return config_path.with_suffix("")


def ensure_executable_exists(executable_path):
    if not executable_path.exists():
        raise FileNotFoundError(f"Executable introuvable: {executable_path}")


def ensure_python_exists():
    if not PYTHON_EXECUTABLE.exists():
        raise FileNotFoundError(f"Python introuvable: {PYTHON_EXECUTABLE}")


def command_for_config(config_path, mode):
    if mode["nproc"] == 1:
        ensure_executable_exists(SEQUENTIAL_EXECUTABLE)
        return [str(SEQUENTIAL_EXECUTABLE), str(config_path)]

    ensure_executable_exists(MPI_EXECUTABLE)
    ensure_executable_exists(MPIEXEC_EXECUTABLE)
    return [str(MPIEXEC_EXECUTABLE), "-n", str(mode["nproc"]), str(MPI_EXECUTABLE), str(config_path)]


def visualization_commands_for_mode(config_path, mode):
    ensure_python_exists()

    if mode["nproc"] == 1:
        ensure_executable_exists(SEQUENTIAL_VISUALIZATION_SCRIPT)
        return [[str(PYTHON_EXECUTABLE), str(SEQUENTIAL_VISUALIZATION_SCRIPT), str(config_path)]]

    ensure_executable_exists(MPI_GLOBAL_VISUALIZATION_SCRIPT)
    ensure_executable_exists(MPI_LOCAL_VISUALIZATION_SCRIPT)
    return [
        [str(PYTHON_EXECUTABLE), str(MPI_GLOBAL_VISUALIZATION_SCRIPT)],
        [str(PYTHON_EXECUTABLE), str(MPI_LOCAL_VISUALIZATION_SCRIPT)],
    ]


def generate_result_visualizations(config_path, output_dir, mode):
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"

    for command in visualization_commands_for_mode(config_path, mode):
        subprocess.run(command, cwd=output_dir, env=env, check=True)


def execute_config(config_path):
    mode = execution_mode_for_config(config_path)
    output_dir = output_directory_for_config(config_path)

    if CLEAN_OUTPUT_DIRECTORIES and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(config_path, output_dir / "config.txt")

    command = command_for_config(config_path, mode)
    print(f"run {config_path.relative_to(SCRIPT_ROOT)} -> {output_dir.relative_to(SCRIPT_ROOT)}")
    subprocess.run(command, cwd=output_dir, check=True)
    if GENERATE_RESULT_IMAGES:
        generate_result_visualizations(config_path, output_dir, mode)


def execute_selected_configs():
    config_files = resolve_config_files()
    if not config_files:
        print("no config selected for execution")
        return

    for config_path in config_files:
        execute_config(config_path)


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

    if RUN_SOLVERS:
        execute_selected_configs()

    print("OK")


if __name__ == "__main__":
    main()
