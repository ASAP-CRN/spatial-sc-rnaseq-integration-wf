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

    model = ModelCls.load(args.checkpoint_dir)
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
