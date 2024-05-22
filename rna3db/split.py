import random

from typing import Sequence, Optional

from rna3db.utils import PathLike, read_json, write_json, read_struct_list


def find_optimal_components(lengths_dict, capacity):
    component_name = list(lengths_dict.keys())
    lengths = list(lengths_dict.values())

    dp = [0] * (capacity + 1)
    trace = [[] for i in range(capacity + 1)]
    for i in range(len(lengths)):
        for j in range(capacity, lengths[i] - 1, -1):
            if dp[j] < dp[j - lengths[i]] + lengths[i]:
                dp[j] = dp[j - lengths[i]] + lengths[i]
                trace[j] = trace[j - lengths[i]] + [component_name[i]]

    return set(trace[capacity])


def find_component_id(cluster_json: dict, struct_id: str):
    for component_id in cluster_json:
        for cluster_id in cluster_json[component_id]:
            if struct_id in cluster_json[component_id][cluster_id]:
                return component_id

    return None


def split(
    input_path: PathLike,
    output_path: PathLike,
    splits: Sequence[float] = [0.7, 0.0, 0.3],
    split_names: Sequence[str] = ["train_set", "valid_set", "test_set"],
    shuffle: bool = False,
    forced_ids_last_path: Optional[PathLike] = None,
    force_zero_last: bool = False,
):
    """A function that splits a JSON of components into a train/test set.

    The split is done by assigning components into the training set until a
    specified training set split percentage (train_size) is met. This is done
    starting with the largest component.

    Args:
        input_path (PathLike): path to JSON containing components
        output_path (PathLike): path to output JSON
    """
    if sum(splits) != 1.0:
        raise ValueError("Sum of splits must equal 1.0.")

    # read json
    cluster_json = read_json(input_path)

    # count number of repr sequences
    lengths = {k: len(v) for k, v in cluster_json.items()}
    total_repr_clusters = sum(lengths.values())

    # shuffle if we want to add randomness
    if shuffle:
        L = list(zip(component_name, lengths))
        random.shuffle(L)
        component_name, lengths = zip(*L)
        component_name, lengths = list(component_name), list(lengths)

    output = {k: {} for k in split_names}

    if force_zero_last:
        output[split_names[-1]]["component_0"] = cluster_json["component_0"]
        lengths.pop("component_0")

    # Force given structures into last split (test)
    if forced_ids_last_path is not None:
        struct_ids = read_struct_list(forced_ids_last_path)

        for struct_id in struct_ids:
            component_id = find_component_id(cluster_json, struct_id)
            if component_id not in output[split_names[-1]]:
                output[split_names[-1]][component_id] = \
                    cluster_json[component_id]
                lengths.pop(component_id)

    capacities = [round(total_repr_clusters * ratio) for ratio in splits]

    if split_names[-1] in output:
        capacities[-1] = \
            capacities[-1] - \
            sum(
                [
                    len(output[split_names[-1]][component_id])
                    for component_id in output[split_names[-1]]
                ]
            )  # Maybe not needed?

    for name, capacity in zip(split_names, capacities):
        components = find_optimal_components(lengths, capacity)
        for k in sorted(components):
            lengths.pop(k)
            output[name][k] = cluster_json[k]

    assert len(lengths) == 0

    write_json(output, output_path)
