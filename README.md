# An astrocytic intracellular calcium signaling model with a focus on calcium refractory period

This repository contains the code and data for an astrocyte intracellular calcium signaling model that explains the **refractory period** observed in experimental recordings.

---

## Repository Structure

### `functions/`
This directory contains the **MATLAB implementation** of the model and the **parameter optimization** code using a genetic algorithm. The code is organized into several subfolders based on functionality.

- The calcium signaling model is defined by a system of **ordinary differential equations (ODEs)**.
- A custom implementation of the **4th-order Runge-Kutta method** is used for simulation.
- The simulation loop can be found in `functions/run.m`.

### `Data/`
This directory contains **preprocessed experimental data** used for model fitting. The structure is as follows:

- Four subfolders represent data from **four individual mice**.
- Each mouse folder contains subfolders for **recording locations** on the motor cortex.
- Each recording location folder includes subfolders for **individual calcium recordings**.

Inside each calcium recording folder:
- A `.mat` file stores:
  - `avg_roi_temp_down`: the real calcium signal extracted from imaging data
  - `velocity`: the animal’s running speed time series
  - `TT`: a table containing visual stimuli and reward event timestamps
- Additionally, precomputed **model input variables** are included to enable faster simulation runs.

---

## How to Run

### 🧪 **Run a Simulation Demo**

1. Download or clone this repository.
2. Open `main_simulation.m` in MATLAB.
3. Set the following parameters in the script:
   - Mouse name
   - Recording location
   - Recording number  
   (Refer to the structure in the `Data/` folder)
4. Run the simulation and compare the model output with the experimental calcium signal.

---

### 🧬 **Parameter Optimization via Genetic Algorithm**

1. Open `main_param_opt.m` in MATLAB.
2. Set the following parameters:
   - `save_path`: directory to store train-test splits and cross-validation data
   - Paths to folders for intermediate GA results
   - GA parameters: `pop_size`, `tournament_size`, `max_generations`, `mutation_rate`
3. Run the script to start optimization.

> 💡 **Performance Note:**  
> On an Ubuntu server with four **Intel(R) Xeon(R) Platinum 8268 CPU @ 2.90GHz** (i.e., 96 workers to be used in MATLAB `parpool`), one generation of the genetic algorithm takes approximately **9000 seconds**, using:
> - `pop_size = 40`
> - `tournament_size = 15`
> - 28 recordings for training

---

## Tested MATLAB Version

The code has been tested and validated on **MATLAB** versions **R2023a, R2024b, R2025a**

---

## Citation

If you use this code or data in your research, please cite the corresponding publication (to be added here when available).
