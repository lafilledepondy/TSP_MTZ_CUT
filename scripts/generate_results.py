import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT / "build"
EXE = BUILD_DIR / "TSP_Gurobi"
DATA_DIR = ROOT / "data"
RESULTS_TEX = ROOT / "results.tex"

MODES = [
    ("MTZ", "MTZ"),
    ("CUT", "CUT"),
    ("CUT_LP", "CUT_LP"),
]


def parse_result_line(line: str) -> dict:
    parts = line.strip().split()
    if not parts or parts[0] != "RESULT":
        return {}
    data = {}
    for part in parts[1:]:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        data[key] = value
    return data


def run_instance(instance_path: Path, mode: str) -> dict:
    cmd = [str(EXE), str(instance_path), mode, "--summary"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )

    for line in (proc.stdout + "\n" + proc.stderr).splitlines():
        if line.startswith("RESULT "):
            data = parse_result_line(line)
            if data:
                return data

    raise RuntimeError(f"No RESULT line found for {instance_path.name} ({mode}).")


def format_num(value: str, decimals: int = 2) -> str:
    if value == "NA":
        return "-"
    try:
        num = float(value)
    except ValueError:
        return value
    if abs(num - round(num)) < 1e-6:
        return str(int(round(num)))
    return f"{num:.{decimals}f}"


def build_table(rows: list) -> str:
    lines = []
    lines.append("\\begin{table}[H]")
    lines.append("\\centering")
    lines.append("\\small")
    lines.append("\\setlength{\\tabcolsep}{3pt}")
    lines.append("\\resizebox{\\textwidth}{!}{")
    lines.append("\\begin{tabular}{l ccccc ccccc ccccc}")
    lines.append("\\toprule")
    lines.append(" & \\multicolumn{5}{c}{MTZ} & \\multicolumn{5}{c}{CUT (sol $\\in \\mathbb{N}$)} & \\multicolumn{5}{c}{CUT (sol $\\in \\mathbb{Q}$)} \\\\")
    lines.append("\\cmidrule(lr){2-6}\\cmidrule(lr){7-11}\\cmidrule(lr){12-16}")
    lines.append("Instance & Obj & Bound & Nodes & Status & Time (s) & Obj & Bound & Cuts & Status & Time (s) & Obj & Bound & Nodes & Status & Time (s)\\\\")
    lines.append("\\midrule")
    lines.extend(rows)
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    lines.append("}")
    lines.append("\\caption{MTZ vs CUT results (time limit 180 secondes)}")
    lines.append("\\end{table}")
    return "\n".join(lines) + "\n"


def main() -> int:
    if not EXE.exists():
        print(f"Missing executable: {EXE}")
        print("Build first: mkdir -p build && cd build && cmake .. && make TSP_Gurobi")
        return 1

    instances = sorted(
        [p for p in DATA_DIR.iterdir() if p.suffix.lower() in {".tsp", ".atsp"}]
    )
    if not instances:
        print("No instances found in data/")
        return 1

    rows = []
    for inst in instances:
        results = {}
        for _, mode in MODES:
            results[mode] = run_instance(inst, mode)

        mtz = results["MTZ"]
        cut = results["CUT"]
        cut_lp = results["CUT_LP"]

        row = " & ".join(
            [
                inst.name,
                format_num(mtz.get("obj", "NA")),
                format_num(mtz.get("bound", "NA")),
                format_num(mtz.get("nodes", "NA"), decimals=0),
                mtz.get("status", "NA"),
                format_num(mtz.get("time", "NA")),
                format_num(cut.get("obj", "NA")),
                format_num(cut.get("bound", "NA")),
                format_num(cut.get("cuts", "NA"), decimals=0),
                cut.get("status", "NA"),
                format_num(cut.get("time", "NA")),
                format_num(cut_lp.get("obj", "NA")),
                format_num(cut_lp.get("bound", "NA")),
                format_num(cut_lp.get("nodes", "NA"), decimals=0),
                cut_lp.get("status", "NA"),
                format_num(cut_lp.get("time", "NA")),
            ]
        )
        rows.append(row + " \\\\")

    RESULTS_TEX.write_text(build_table(rows), encoding="utf-8")
    print(f"Wrote {RESULTS_TEX}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
