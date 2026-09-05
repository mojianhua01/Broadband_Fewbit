# Broadband Few-Bit

This repository contains the MATLAB code for the paper:

> **Channel Estimation in Broadband Millimeter-Wave MIMO Systems with Few-Bit ADCs**

The code evaluates broadband millimeter-wave MIMO channel-estimation methods under low-resolution ADC quantization. It includes the simulation scripts, channel data, pilot-sequence utilities, antenna-radiation examples, and the GAMP implementations used by the experiments.

## Requirements

- MATLAB with support for matrix and signal-processing functions used by the scripts
- The repository root and its subdirectories available on the MATLAB path

## Quick Start

1. Clone or download this repository.
2. Open MATLAB and change the current folder to the repository root.
3. Add the repository folders to the MATLAB path:

   ```matlab
   addpath(genpath(pwd))
   ```
4. Run the main comparison script:

   ```matlab
   cd Broadband_Fewbit_Main
   Broadband_Few_Bit_Comparison
   ```

The default configuration uses the included `H_UPA_16_4_16_2clusters.mat` channel file, one-bit ADC quantization, a shifted Zadoff-Chu training sequence, and fixed random seed `42`. Simulation parameters can be changed near the top of `Broadband_Few_Bit_Comparison.m`.

## Repository Structure

| Directory                  | Description                                                             |
| -------------------------- | ----------------------------------------------------------------------- |
| `Broadband_Fewbit_Main/` | Main broadband few-bit channel-estimation experiments and channel data  |
| `main/`                  | GAMP estimators, linear transforms, and input/output estimation modules |
| `Phil_Pilot/`            | Pilot-generation and training utilities                                 |
| `Gold Sequences/`        | Gold-sequence generation and correlation utilities                      |
| `Narrowband_Fewbit/`     | Narrowband few-bit experiments                                          |
| `Antenna_Radiation/`     | Antenna and array-radiation examples                                    |
| `Test_of_*/`             | Additional tests and experimental scripts                               |
| `WLAN/`                  | WLAN-related utilities and experiments                                  |

## Reproducibility

The main script initializes MATLAB's random-number generator with `rng(42)`. To run a different experiment, update the model, ADC, pilot, SNR, and training-length parameters in the script before execution. Several alternative configurations are provided as commented examples.

## Citation

If you use this code in academic work, please cite the associated paper. The complete publication information is:

> J. Mo, P. Schniter, and R. W. Heath, "Channel Estimation in Broadband Millimeter Wave MIMO Systems With Few-Bit ADCs," *IEEE Transactions on Signal Processing*, vol. 66, no. 5, pp. 1141-1154, Mar. 2018, doi: [10.1109/TSP.2017.2781644](https://doi.org/10.1109/TSP.2017.2781644).

**Authors:** Jianhua Mo, Philip Schniter, and Robert W. Heath

**Journal:** IEEE Transactions on Signal Processing

**Publication details:** Volume 66, Issue 5, pages 1141-1154, published March 1, 2018

**DOI:** [10.1109/TSP.2017.2781644](https://doi.org/10.1109/TSP.2017.2781644)

### BibTeX

```bibtex
@article{Mo2018BroadbandFewBit,
	author  = {Mo, Jianhua and Schniter, Philip and Heath, Robert W.},
	title   = {Channel Estimation in Broadband Millimeter Wave MIMO Systems With Few-Bit ADCs},
	journal = {IEEE Transactions on Signal Processing},
	volume  = {66},
	number  = {5},
	pages   = {1141--1154},
	year    = {2018},
	month   = mar,
	doi     = {10.1109/TSP.2017.2781644}
}
```
