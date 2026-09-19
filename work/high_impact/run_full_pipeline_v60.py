#!/usr/bin/env python3
"""Run the v5 reproducibility pipeline followed by v6 prespecified upgrades."""
from pathlib import Path
import subprocess, sys

ROOT=Path(__file__).resolve().parents[2]
commands=[
    [sys.executable,'work/high_impact/run_full_pipeline_v50.py'],
    [sys.executable,'work/high_impact/run_temporal_transportability.py'],
    [sys.executable,'work/high_impact/run_multidonor_msc_atlas.py'],
    [sys.executable,'work/high_impact/make_figure5_top_tier.py'],
]
for command in commands:
    print('+',' '.join(command),flush=True)
    subprocess.run(command,cwd=ROOT,check=True)
