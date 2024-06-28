import tempfile
import subprocess
from pathlib import Path
import tqdm

from rna3db.utils import PathLike, read_json, write_json

DSSR_EXECUTABLE_PATH = "<ADD PATH HERE>"


def run_dssr(mmcif_path: PathLike) -> dict:
    sec_structs = {}

    with tempfile.NamedTemporaryFile() as f_tmp:
        # Run DSSR
        subprocess.check_call(
            [
                DSSR_EXECUTABLE_PATH,
                f"--input={str(mmcif_path)}",
                f"--output={f_tmp.name}"
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )

        # Parse DSSR output
        lines = f_tmp.readlines()
        for idx, line in enumerate(lines):
            if line.startswith(">") and "#" in line:
                chain_id = line.split()[0].split("-")[-1]
                seq = lines[idx + 1]
                ss = lines[idx + 2]

                # Remove ligand pairings
                while seq[-2] == "&" and seq[-1].islower():
                    seq = seq[:-2]
                    ss = ss[:-2]

                sec_structs[chain_id] = ss

    return sec_structs


def dssr(
    input_path: PathLike,
    output_path: PathLike,
    mmcif_dir: PathLike
) -> None:
    data_json = read_json(input_path)
    mmcif_dir = Path(mmcif_dir)
    ss_cache = {}

    for component_id in tqdm(data_json):
        for cluster_id in data_json[component_id]:
            for struct_id in data_json[component_id][cluster_id]:
                entry_id, chain_id = struct_id.split("_")

                if entry_id not in ss_cache:
                    ss_cache[entry_id] = \
                        run_dssr(mmcif_dir / f"{entry_id}.cif")

                if chain_id not in ss_cache[entry_id]:
                    ss_cache[entry_id][chain_id] = ""

                data_json[component_id][cluster_id][struct_id]["ss"] = \
                    ss_cache[entry_id][chain_id]

    write_json(data_json, output_path)
