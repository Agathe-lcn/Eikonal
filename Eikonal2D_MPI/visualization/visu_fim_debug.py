import glob
import os
import re
import matplotlib.pyplot as plt
import numpy as np


def parse_trace_files():
    traces = []
    for path in sorted(glob.glob("mpi_exchange_rank*.txt")):
        rank = int(re.search(r"rank(\d+)", path).group(1))
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                entry = {}
                for token in parts:
                    if "=" not in token:
                        continue
                    key, value = token.split("=", 1)

                    # Conversion automatique des types
                    if key in {"cycle", "comm", "rank", "updated"}:
                        entry[key] = int(value)
                    elif key in {"value", "old_T", "new_T"}:
                        entry[key] = float(value)
                    else:
                        entry[key] = value

                if "direction" in entry and "stage" in entry:
                    # Garantir que le rang est présent
                    if "rank" not in entry:
                        entry["rank"] = rank
                    traces.append(entry)
    return traces


def summarize_by_cycle(traces):
    by_cycle = {}
    for entry in traces:
        key = (entry["cycle"], entry["comm"], entry["rank"])
        by_cycle.setdefault(key, []).append(entry)
    return by_cycle


def main():
    traces = parse_trace_files()
    if not traces:
        print("Aucune trace mpi_exchange_rank*.txt trouvée.")
        return

    by_cycle = summarize_by_cycle(traces)
    print(
        f"{len(traces)} lignes de trace lues sur {len(by_cycle)} couples (cycle, comm, rank)."
    )

    for (cycle, comm, rank), entries in sorted(by_cycle.items())[:20]:
        print(
            f"cycle={cycle} comm={comm} rank={rank} -> {len(entries)} événements"
        )
        for e in entries[:8]:
            print(" ", e)
        if len(entries) > 8:
            print("  ...")

    # Graphique des mises à jour de halos
    updates = [e for e in traces if e.get("updated") == 1]
    if updates:
        fig, ax = plt.subplots(figsize=(10, 4))
        ranks = [e["rank"] for e in updates]
        cycles = [e["cycle"] for e in updates]
        ax.scatter(cycles, ranks, c="tab:red", s=20, alpha=0.8)
        ax.set_xlabel("cycle")
        ax.set_ylabel("rank")
        ax.set_title("Mises à jour de halos observées")
        fig.tight_layout()
        fig.savefig("mpi_exchange_updates.png", dpi=200, bbox_inches="tight")
        print("Figure enregistrée: mpi_exchange_updates.png")
    else:
        print("Aucune mise à jour de halo (updated == 1) détectée.")


if __name__ == "__main__":
    main()