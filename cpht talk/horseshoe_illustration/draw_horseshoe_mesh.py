from pathlib import Path

import matplotlib.pyplot as plt


def parse_inp(path: Path):
    nodes = {}
    elements = []
    section = None

    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue

        if line.startswith("*"):
            upper = line.upper()
            if upper.startswith("*NODE"):
                section = "node"
            elif upper.startswith("*ELEMENT") and ("CPS3" in upper or "CPS4" in upper):
                section = "element"
            else:
                section = None
            continue

        if section == "node":
            parts = [x.strip() for x in line.split(",") if x.strip()]
            if len(parts) >= 3:
                node_id = int(parts[0])
                nodes[node_id] = (float(parts[1]), float(parts[2]))
            continue

        if section == "element":
            parts = [x.strip() for x in line.split(",") if x.strip()]
            if len(parts) >= 4:
                conn = tuple(int(x) for x in parts[1:])
                elements.append(conn)

    return nodes, elements


def mesh_edges(elements):
    edges = set()
    for conn in elements:
        n = len(conn)
        for i in range(n):
            a = conn[i]
            b = conn[(i + 1) % n]
            if a == b:
                continue
            edges.add(tuple(sorted((a, b))))
    return edges


def main():
    repo_root = Path(__file__).resolve().parents[2]
    inp_path = repo_root / "projects/simple_geometries/meshes/snake/horseshoe.inp"
    nodes, elements = parse_inp(inp_path)
    edges = mesh_edges(elements)

    fig, ax = plt.subplots(figsize=(6, 6), dpi=220)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    for a, b in edges:
        p0 = nodes[a]
        p1 = nodes[b]
        ax.plot(
            [p0[0], p1[0]],
            [p0[1], p1[1]],
            color="white",
            lw=0.45,
            alpha=0.72,
            solid_capstyle="butt",
            zorder=3,
        )

    xs = [x for x, _ in nodes.values()]
    ys = [y for _, y in nodes.values()]
    pad = 0.04
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect("equal")
    ax.axis("off")

    out_path = Path(__file__).with_name("horseshoe_mesh_illustration.png")
    fig.savefig(out_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.05)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
