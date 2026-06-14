# Protein-Oligomerization-Equilibria
# 🔬 Protein Oligomerization Simulator

An interactive web tool for simulating fluorescence correlation spectroscopy (FCS) titration curves across three oligomerization equilibria: monomer–dimer, monomer–trimer, and monomer–dimer–tetramer.

> *This tool is intended to facilitate the implementation of quantitative FCS-based oligomerization assays and to promote the use of physically accurate binding models over empirical approximations that may produce inaccurate or misleading parameter estimates.*

---

## Overview

Determining the oligomeric state of a protein and its concentration-dependent equilibrium is a central challenge in biophysics. FCS provides a powerful, solution-based readout of apparent diffusion times that reflect the average hydrodynamic size of a labeled species, but interpreting titration curves requires physically accurate models that account for labeling efficiency, the concentration of labeled tracer, and the stoichiometry of the assembly.

This tool computes and displays the predicted **normalized apparent diffusion time** τ_n = τ_app / τ_m as a function of total protein concentration, given user-specified dissociation constants (*K*_d) and experimental parameters (*C*_L, *f*). It is designed for research groups using FCS to characterize self-associating proteins and serves as a companion to the theoretical framework described in the references below.

---

## Original Code

This webapp is based on the `<filename>.py` file, maintaining all the original functionality while adding an interactive web interface.

---

## The Three Equilibrium Cases

**Monomer–Dimer** — monomer fraction solved in closed form:

```
2M ⇌ D        (Kd)
```

**Monomer–Trimer** — monomer fraction from roots of a cubic equation; *K*_d here is an effective dissociation constant derived from a per-step equilibrium constant *K*_d,E:

```
3M ⇌ T        (Kd_eff)
```

**Monomer–Dimer–Tetramer** — sequential two-step assembly via a dimeric intermediate; species fractions from roots of a quartic, apparent diffusion time from roots of a cubic:

```
2M ⇌ D ⇌ T₄        (Kd2 for M→D step, Kd1 for D→T₄ step)
```

---

## Parameters

| Parameter | Description |
|-----------|-------------|
| *K*_d | Dissociation constant for the monomer–oligomer equilibrium (nM) |
| *K*_d1 | Dissociation constant for the dimer-to-tetramer step (nM) |
| *K*_d2 | Dissociation constant for the monomer-to-dimer step (nM) |
| *f* | Labeling efficiency (fraction of molecules carrying a fluorescent label) |
| *C*_L | Concentration of labeled protein in terms of the highest oligomer (nM) |
| Conc. range | Total protein concentration range for the simulation (nM) |

The tool outputs:
- **Upper panel:** τ_app / τ_oligomer vs. protein concentration (log scale), with a secondary *x*-axis showing total monomer concentration
- **Lower panel:** Fractional species concentrations (α_monomer, α_dimer, α_trimer, or α_tetramer) vs. protein concentration
- **CSV download:** All plotted values for downstream analysis

---

## Installation

### Requirements

- Python ≥ 3.9
- pip

### Setup

```bash
# 1. Clone the repository
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>

# 2. (Recommended) Create a virtual environment
python -m venv venv
source venv/bin/activate        # macOS / Linux
venv\Scripts\activate           # Windows

# 3. Install dependencies
pip install streamlit numpy matplotlib pandas
```

### Run locally

```bash
streamlit run app.py
```

The app will open automatically in your browser at `http://localhost:8501`.

---

## Usage

1. Select an **equilibrium model** from the sidebar dropdown.
2. Enter the dissociation constant(s) and experimental parameters (*f*, *C*_L).
3. Set the **concentration range** to match the span of your titration experiment.
4. Click **▶ Run Simulation** to update the plots.
5. Use **⬇ Download CSV** to export the simulated curves for plotting or fitting in your analysis software of choice.

> **Tip:** The plots auto-refresh when you switch equilibrium models. For all other parameter changes, click *Run Simulation* to apply.

---

## References

1. **[Author(s), Title, Journal, Year, DOI — placeholder]**

2. Kanno, D. M., & Levitus, M. (2014). Protein oligomerization equilibria and kinetics investigated
   by fluorescence correlation spectroscopy: A mathematical treatment.
   *The Journal of Physical Chemistry B*, 118(43), 12404–12415. https://doi.org/10.1021/jp507741r

---

## Citation

If you use this tool in your work, please cite reference [1] above and link to this repository.

---

## License

[MIT / GPLv3 / other — add your license here]
