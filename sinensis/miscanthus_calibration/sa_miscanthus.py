"""Global sensitivity analysis of WOFOST for Miscanthus at Veenkampen.

Adapted from Allard de Wit's SA_mp_example.py (September 2023).
Runs Sobol' sensitivity analysis in parallel via multiprocessing, then
saves the parameter sets and target results to disk so the analysis and
plotting can be done in a separate notebook.

Usage:
    conda activate py313_pcse
    python sa_miscanthus.py

Outputs (saved next to this script):
    sa_paramsets.npy        Saltelli parameter samples, shape (N*(2k+2), k)
    sa_target_results.npy   WOFOST output for each paramset, shape (N*(2k+2),)
    sa_problem.yaml         Copy of the problem definition for reproducibility

Paulian Oprescu, May 2026.
"""

import sys
from pathlib import Path
from functools import partial, lru_cache
import multiprocessing as mp

import yaml
import numpy as np
from tqdm import tqdm
from SALib.sample import saltelli

import pcse
from pcse.models import Wofost72_PP
from pcse.base import ParameterProvider
from pcse.fileinput import YAMLCropDataProvider, CABOFileReader, ExcelWeatherDataProvider
from pcse.util import WOFOST72SiteDataProvider


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
project = Path("/Users/paulianoprescu/Projects/wofost_miscanthus/miscanthus_calibration")
input_dir = project / "input_params"

CROP_FILE    = input_dir / "miscanthus" / "miscanthus.yaml"
SOIL_FILE    = input_dir / "soil" / "soil_10110.soil"
WEATHER_FILE = input_dir / "weather" / "veenkampen_weather.xlsx"

# Output directory (same folder as this script)
OUTPUT_DIR = Path(__file__).parent


# ---------------------------------------------------------------------------
# Cached loaders — each worker loads inputs once, then reuses them
# ---------------------------------------------------------------------------
@lru_cache
def get_modelparameters():
    cropdata = YAMLCropDataProvider(fpath=CROP_FILE.parent, force_reload=True)
    cropdata.set_active_crop("miscanthus", "miscanthus_sinensis")
    soildata = CABOFileReader(SOIL_FILE)
    sitedata = WOFOST72SiteDataProvider(WAV=70)
    return ParameterProvider(cropdata=cropdata, soildata=soildata, sitedata=sitedata)


@lru_cache
def get_weatherdata():
    return ExcelWeatherDataProvider(str(WEATHER_FILE))


@lru_cache
def get_agromanagement():
    agro_yaml = """
    - 2018-01-01:
        CropCalendar:
            crop_name: miscanthus
            variety_name: miscanthus_sinensis
            crop_start_date: 2018-01-01
            crop_start_type: emergence
            crop_end_date: 2018-12-30
            crop_end_type: harvest
            max_duration: 300
        TimedEvents: null
        StateEvents: null
    """
    return yaml.safe_load(agro_yaml)


# ---------------------------------------------------------------------------
# Single simulation
# ---------------------------------------------------------------------------
def run_wofost_simulation(paramset, problem, target_variable):
    params = get_modelparameters()
    weather = get_weatherdata()
    agromanagement = get_agromanagement()

    params.clear_override()
    for name, value in zip(problem["names"], paramset):
        params.set_override(name, value)

    wofost = Wofost72_PP(params, weather, agromanagement)
    wofost.run_till_terminate()
    r = wofost.get_summary_output()
    target_result = r[0][target_variable]
    if target_result is None:
        print("Target variable is not available in summary output!")
    return target_result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    target_variable = "TAGP"

    problem_yaml = """
        num_vars: 6
        names:
        - TSUM1
        - TSUM2
        - RGRLAI
        - Q10
        - SPAN
        - RDI
        bounds:
        - [1500, 2000]
        - [100, 500]
        - [0.01, 0.1]
        - [1.8, 2.2]
        - [50, 100]
        - [5, 100]
    """
    problem = yaml.safe_load(problem_yaml)

    # Saltelli sample size — power of 2 recommended for Sobol sequences
    calc_second_order = False   # set True only if you need S2 (much more expensive)
    nsamples = 16
    paramsets = saltelli.sample(problem, nsamples, calc_second_order=calc_second_order)
    print(f"Generated {len(paramsets)} parameter sets "
          f"(N={nsamples}, k={problem['num_vars']}, second_order={calc_second_order})")

    run_wofost_partial = partial(run_wofost_simulation,
                                 problem=problem,
                                 target_variable=target_variable)

    n_cpu = mp.cpu_count()
    print(f"Running on {n_cpu} CPU cores...")
    target_results = []
    with tqdm(total=len(paramsets)) as pbar:
        with mp.Pool(n_cpu) as pool:
            for result in pool.imap(run_wofost_partial, paramsets):
                target_results.append(result)
                pbar.update()

    target_results = np.array(target_results)

    # Save outputs for downstream analysis in the notebook
    np.save(OUTPUT_DIR / "sa_paramsets.npy", paramsets)
    np.save(OUTPUT_DIR / "sa_target_results.npy", target_results)
    with open(OUTPUT_DIR / "sa_problem.yaml", "w") as f:
        f.write(problem_yaml)

    print(f"\nSaved:")
    print(f"  {OUTPUT_DIR / 'sa_paramsets.npy'}      shape={paramsets.shape}")
    print(f"  {OUTPUT_DIR / 'sa_target_results.npy'}  shape={target_results.shape}")
    print(f"  {OUTPUT_DIR / 'sa_problem.yaml'}")
    print(f"\nTarget variable: {target_variable}")
    print(f"Mean {target_variable}: {target_results.mean():.1f}   "
          f"Std: {target_results.std():.1f}   "
          f"Min/Max: {target_results.min():.1f} / {target_results.max():.1f}")


if __name__ == "__main__":
    main()
