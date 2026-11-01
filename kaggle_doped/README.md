# Running the doped-instance job on Kaggle

1. Create a Kaggle dataset from this folder (all files), e.g. named `doped-run7`.
2. New notebook, accelerator **TPU VM v3-8** (the host has ~330 GB RAM and 96 vCPUs; the TPU itself is unused), add the dataset.
3. One cell:
   ```
   !bash /kaggle/input/doped-run7/run_kaggle.sh /kaggle/input/doped-run7 500
   ```
   (500 = minutes of budget; a session is capped at 9 h.) Outputs land in `/kaggle/working/doped/`:
   `gate_result.txt` (must say PASS), `shots_scored.txt` (uniform controls first, then IBM's shots, exact log2 Pr each),
   `samples.txt` (ACC lines are the samples; every line has its z), `run.log`.
4. Download `/kaggle/working/doped/` and give it to the session.
The first two passes (the gate) tell the per-pass time; expect minutes per pass on 96 cores.
