#!/usr/bin/env python3
"""Load a saved Cell2location model checkpoint and run `export_posterior` with
`accelerator="cpu"`, writing the resulting AnnData to disk.

This is NOT a user-facing CLI -- it's invoked by `cell2location_deconv.py` as a
subprocess to avoid all-NaN posteriors caused by Pyro's parameter store retaining training
state. Needed to avoid a leftover GPU "memory" and instead run a fresh subprocess."""

import argparse


def export_posterior(
    checkpoint_dir: str,
    model_class: str,
    output_h5ad: str,
    num_samples: int,
    batch_size: int,
) -> None:
    """
    Load a model checkpoint, run `export_posterior`, and write the result to disk.

    Parameters
    ----------
    checkpoint_dir : str
        Directory containing the saved model checkpoint.
    model_class : str
        "RegressionModel" or "Cell2location"; selects which class to load.
    output_h5ad : str
        Destination path for the posterior-exported AnnData.
    num_samples : int
        Posterior samples to draw.
    batch_size : int
        Batch size for posterior sampling.

    Returns
    -------
    None
    """
    if model_class == "RegressionModel":
        from cell2location.models import RegressionModel as ModelCls
    else:
        from cell2location.models import Cell2location as ModelCls

    # accelerator="cpu" here (not just in sample_kwargs below) is required: .load() defaults
    # to "auto", which picks CUDA when available, leaving the model's own buffers (e.g.
    # detection_mean_y_e) on cuda:0 while sample_posterior's "cpu" accelerator creates fresh
    # tensors on CPU — a device mismatch that only surfaces on an actual CUDA GPU, not MPS/CPU.
    model = ModelCls.load(checkpoint_dir, accelerator="cpu")
    adata = model.export_posterior(
        model.adata,
        sample_kwargs={
            "num_samples": num_samples,
            "batch_size": batch_size,
            "accelerator": "cpu",
        },
    )
    adata.write_h5ad(output_h5ad, compression="lzf")
    print(f"[INFO] Wrote posterior-exported adata to {output_h5ad}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--checkpoint-dir",
        required=True,
        help="Directory containing the saved model checkpoint.",
    )
    parser.add_argument(
        "--model-class",
        required=True,
        choices=["RegressionModel", "Cell2location"],
        help="Which class to load the checkpoint as.",
    )
    parser.add_argument(
        "--output-h5ad",
        required=True,
        help="Destination path for the posterior-exported AnnData.",
    )
    parser.add_argument(
        "--posterior-num-samples",
        type=int,
        required=True,
        help="Posterior samples to draw.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        required=True,
        help="Batch size for posterior sampling.",
    )
    args = parser.parse_args()

    export_posterior(
        checkpoint_dir=args.checkpoint_dir,
        model_class=args.model_class,
        output_h5ad=args.output_h5ad,
        num_samples=args.posterior_num_samples,
        batch_size=args.batch_size,
    )


if __name__ == "__main__":
    main()
