#!/usr/bin/env python3
"""
Load a saved Cell2location model checkpoint and run `export_posterior` with
`accelerator="cpu"`; write the resulting AnnData to disk.

Not a user-facing CLI — invoked by `cell2location_deconvolute.py` as a subprocess
to avoid all-NaN posteriors caused by Pyro's parameter store retaining training state.
"""

import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint-dir", required=True)
    parser.add_argument("--model-class", required=True, choices=["RegressionModel", "Cell2location"])
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--posterior-num-samples", type=int, required=True)
    parser.add_argument("--batch-size", type=int, required=True)
    args = parser.parse_args()

    if args.model_class == "RegressionModel":
        from cell2location.models import RegressionModel as ModelCls
    else:
        from cell2location.models import Cell2location as ModelCls

    # accelerator="cpu" here (not just in sample_kwargs below) is required: .load() defaults
    # to "auto", which picks CUDA when available, leaving the model's own buffers (e.g.
    # detection_mean_y_e) on cuda:0 while sample_posterior's "cpu" accelerator creates fresh
    # tensors on CPU — a device mismatch that only surfaces on an actual CUDA GPU, not MPS/CPU.
    model = ModelCls.load(args.checkpoint_dir, accelerator="cpu")
    adata = model.export_posterior(
        model.adata,
        sample_kwargs={
            "num_samples": args.posterior_num_samples,
            "batch_size": args.batch_size,
            "accelerator": "cpu",
        },
    )
    adata.write_h5ad(args.output_h5ad, compression="gzip")
    print(f"[INFO] Wrote posterior-exported adata to {args.output_h5ad}")


if __name__ == "__main__":
    main()
