import subprocess
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

from rna3db.utils import PathLike, read_json, write_json

DSSR_EXECUTABLE_PATH = "<ADD DSSR EXECUTABLE HERE>"
MAX_DSSR_TRIES = 3


def run_dssr(cif_file: PathLike, out_file: PathLike) -> None:
    if Path(out_file).exists():
        return

    for _ in range(MAX_DSSR_TRIES):
        try:
            subprocess.check_call(
                [
                    DSSR_EXECUTABLE_PATH,
                    f"--input={str(cif_file)}",
                    f"--output={str(out_file)}"
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
            return
        except Exception:
            continue


def parse_dssr_output(dssr_out_file: PathLike) -> dict[str, str]:
    sec_structs = {}

    with open(dssr_out_file, "r") as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if line.startswith(">") and "#" in line:
            chain_id = line.split()[0].split("-")[-1]
            seq = lines[idx + 1].rstrip()
            ss = lines[idx + 2].rstrip()

            if len(seq) < 2:
                continue

            # Remove ligand pairings
            while seq[-2] == "&" and seq[-1].islower():
                seq = seq[:-2]
                ss = ss[:-2]

            sec_structs[chain_id] = ss

    return sec_structs


def collect_pdb_ids(data_json: dict) -> set[str]:
    pdb_ids = set()

    for component_id in data_json:
        for cluster_id in data_json[component_id]:
            for struct_id in data_json[component_id][cluster_id]:
                pdb_id = struct_id.split("_")[0]
                pdb_ids.add(pdb_id)

    return pdb_ids


def calculate_ss(
    input_path: PathLike,
    output_path: PathLike,
    cif_dir: PathLike,
    dssr_out_dir: PathLike,
    num_workers: int = 1
) -> None:
    data_json = read_json(input_path)

    cif_dir = Path(cif_dir)
    dssr_out_dir = Path(dssr_out_dir)

    dssr_out_dir.mkdir(parents=True, exist_ok=True)
    
    pdb_ids = list(collect_pdb_ids(data_json))

    with ThreadPoolExecutor(max_workers=num_workers) as ex:
        _ = list(
                tqdm(
                    ex.map(
                        run_dssr,
                        [cif_dir / f"{pdb_id}.cif" for pdb_id in pdb_ids],
                        [dssr_out_dir / f"{pdb_id}.dssr" for pdb_id in pdb_ids],
                    ), total=len(pdb_ids)
            )
        )

    for component_id in data_json:
        for cluster_id in data_json[component_id]:
            for struct_id in data_json[component_id][cluster_id]:
                pdb_id, chain_id = struct_id.split("_")

                data_json[component_id][cluster_id][struct_id]["ss"] = ""

                dssr_out = parse_dssr_output(dssr_out_dir / f"{pdb_id}.dssr")
                data_json[component_id][cluster_id][struct_id]["ss"] = ""

                if chain_id in dssr_out:
                    data_json[component_id][cluster_id][struct_id]["ss"] = \
                        dssr_out[chain_id]

    write_json(data_json, output_path)
