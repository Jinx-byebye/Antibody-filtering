"""Lossless confidence output for the frozen OpenDDE benchmark."""
from pathlib import Path
import numpy as np


def flatten(value, prefix="", result=None):
    """Serialize tensor trees without rounding or object/pickle arrays."""
    if result is None:
        result = {}
    if isinstance(value, dict):
        for key, item in value.items():
            flatten(item, f"{prefix}__{key}" if prefix else str(key), result)
    elif value is not None:
        if hasattr(value, "detach"):
            value = value.detach().cpu()
            # NumPy has no native bfloat16; float32 represents every value exactly.
            if str(value.dtype) == "torch.bfloat16":
                value = value.float()
            value = value.numpy()
        array = np.asarray(value)
        if array.dtype.hasobject:
            raise TypeError(f"Unsupported confidence field: {prefix}")
        result[prefix] = array
    return result


def install_lossless_dumper():
    import runner.dumper as module
    from opendde.utils.file_io import save_json

    class LosslessDumper(module.DataDumper):
        def dump_predictions(self, pred_dict, dump_dir, pdb_id, atom_array,
                             entity_poly_type, seed):
            self._benchmark_atoms = atom_array
            return super().dump_predictions(pred_dict, dump_dir, pdb_id, atom_array,
                                            entity_poly_type, seed)

        def _save_confidence(self, data, prediction_save_dir, sample_name, seed,
                             sorted_indices):
            if not self.need_atom_confidence or not data.get("full_data"):
                raise ValueError("Full confidence is mandatory for this benchmark.")
            ranks = sorted_indices if sorted_indices is not None else range(len(data["full_data"]))
            atoms = self._benchmark_atoms
            for index, rank in enumerate(ranks):
                arrays = flatten(data["full_data"][index])
                arrays["sample_index"] = np.asarray(index)
                arrays["seed"] = np.asarray(seed)
                arrays["coordinate"] = flatten(data["coordinate"][index], "value")["value"]
                for key in ("contact_probs", "per_sample_contact_probs"):
                    if key in data:
                        value = data[key]
                        if key == "per_sample_contact_probs":
                            value = value[index]
                        arrays[f"model__{key}"] = flatten(value, "value")["value"]
                for key in ("chain_id", "res_id", "res_name", "atom_name", "element"):
                    arrays[f"structure__{key}"] = np.asarray(getattr(atoms, key))
                output = Path(prediction_save_dir)
                np.savez_compressed(output / f"{sample_name}_full_data_sample_{rank}.npz", **arrays)
                save_json(data["summary_confidence"][index],
                          str(output / f"{sample_name}_summary_confidence_sample_{rank}.json"), indent=4)

    module.DataDumper = LosslessDumper
