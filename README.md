# Protein-Oligomerization-Equilibria
# 🔬 Protein Oligomerization Simulator

An interactive web tool for simulating fluorescence correlation spectroscopy (FCS) titration curves for three oligomerization equilibria: monomer–dimer, monomer–trimer, and monomer–dimer–tetramer.

> *This tool is intended to facilitate the implementation of quantitative FCS-based oligomerization studies and to promote the use of physically accurate binding models over empirical approximations that may result in inaccurate or misleading parameter estimates.*

---

## Overview

FCS measures fluorescence intensity fluctuations as labeled molecules diffuse through a femtoliter-sized confocal observation volume. The temporal behavior of these fluctuations, captured by the autocorrelation function, encodes the translational diffusion coefficients of the fluorescent species in the sample. Since the diffusion coefficient of a globular protein scales approximately with the cube root of its molecular mass, FCS can, in principle, distinguish between oligomeric states and report on the degree of dissociation as a function of concentration.
 
A key complication is that the autocorrelation decay of a mixture of oligomeric species is experimentally indistinguishable from that of a single species. The measured quantity is therefore an apparent diffusion time (τₐₚₚ) that reflects the weighted average of all species present, where the weighting depends on both the concentrations and relative brightnesses of each species. Because larger oligomers carry more fluorescent labels and contribute more strongly to the autocorrelation signal, simple intuitive interpretations of τₐₚₚ can be seriously misleading without a proper mathematical framework.
 
This tool is based on a rigorous mathematical treatment that establishes analytical relationships between the experimentally measured τₐₚₚ values and the thermodynamic parameters governing protein oligomerization, including dissociation equilibrium constants for three different oligomerization cases: dimer-monomer, trimer-monomer and tetramer-dimer-monomer. It computes the predicted **normalized apparent diffusion time** τₙ = τₐₚₚ / τₘ as a function of total protein concentration, given user-specified dissociation constants (Kd) and experimental parameters (CL, *f*).

---

## Original Code

This webapp is based on the `proteinequilibriafn.py` file, maintaining all the original functionality while adding an interactive web interface.

---

## The Three Equilibrium Cases

**Monomer–Dimer**: monomer fraction solved in closed form:

```
D ⇌ 2M        (K_d)
```

**Monomer–Trimer**: monomer fraction from roots of a cubic equation; Kd here is an effective dissociation constant derived from a per-step equilibrium constant Kd,E:

```
T ⇌ 3M        (K_d,eff)
```

**Monomer–Dimer–Tetramer**: sequential two-step assembly via a dimeric intermediate; species fractions from roots of a quartic, apparent diffusion time from roots of a cubic:

```
T₄ ⇌ 2D ⇌ 4M        (K_d1 for T₄→D step, K_d2 for D→M step)
```

---

## Parameters

| Parameter | Description |
|-----------|-------------|
| Kd | Dissociation constant for the monomer–oligomer equilibrium (nM); for the trimer case this is an effective Kd derived from a per-step equilibrium constant Kd,E |
| Kd1 | Dissociation constant for the dimer-to-tetramer step (nM) |
| Kd2 | Dissociation constant for the monomer-to-dimer step (nM) |
| *f* | Labeling efficiency (fraction of molecules carrying a fluorescent label) |
| CL | Concentration of labeled protein in terms of the highest oligomer (nM) |
| Conc. range | Total protein concentration range for the simulation (nM) |

The tool outputs:
- **Upper panel:** τₐₚₚ / τₘ vs. protein concentration (log scale), with a secondary *x*-axis showing total monomer concentration
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
2. Enter the dissociation constant(s) and experimental parameters (*f*, CL).
3. Set the **concentration range** to match the span of your titration experiment.
4. Click ** Run Simulation** to update the plots.
5. Use ** Download CSV** to export the simulated curves for plotting or fitting in your analysis software of choice.

> **Tip:** The plots auto-refresh when you switch equilibrium models. For all other parameter changes, click *Run Simulation* to apply.

---

## References

1. **[Author(s), Title, Journal, Year, DOI: placeholder]**

2. Kanno, D. M., & Levitus, M. (2014). Protein oligomerization equilibria and kinetics investigated
   by fluorescence correlation spectroscopy: A mathematical treatment.
   *The Journal of Physical Chemistry B*, 118(43), 12404–12415. https://doi.org/10.1021/jp507741r

---

## Citation

If you use this tool in your work, please cite reference [1] above.

---

## License

[MIT / GPLv3 / other: add your license here]
