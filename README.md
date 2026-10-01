# StabilityFunction

**StabilityFunction** is a SageMath package for computing GIT-stability and GIT-semistable reduction of projective plane curves ($n=2$). Its main functionality addresses three problems:

1. **Determine GIT-stability:** Decide whether a projective plane curve over a perfect field is GIT-semistable and whether it is GIT-stable.
2. **Find a model with GIT-semistable reduction:** Given a GIT-stable plane curve $X$ over a complete discretely valued field $K$ with perfect residue field, compute a plane model over the valuation ring of $K$ whose special fiber is GIT-semistable, or determine that no such model exists.
3. **Find an extension admitting GIT-semistable reduction:** When the residue field is finite, compute a finite extension $L/K$ such that $X_L$ admits a plane model with GIT-semistable special fiber.

The package also provides a function combining the second and third tasks: it computes a suitable extension $L/K$ together with a plane model of $X_L$ whose special fiber is GIT-semistable.

Additional functionality supports computations with diagonalizable valuations, Bruhat–Tits buildings, and the geometry of plane curves. Example notebooks are available in the `notebooks/` directory.

## Project Structure

* **semistable_model/**: The main Python package containing the mathematical logic.
* **notebooks/**: Interactive Jupyter notebooks demonstrating the usage of the project.
* **documentation/**: Contains supplementary documentation.

```text
StabilityFunction
├── documentation
│   └── commutative_diagrams
│       ├── base_change_alg.tex
│       ├── base_change_plane_model.tex
│       ├── base_change_vec.tex
│       └── basis_matrix_product.tex
├── notebooks
│   ├── example_plane_models.ipynb
│   ├── example_semistable_model_with_irrational_cusps.ipynb
│   ├── example_semistable_model_with_rational_cusps.ipynb
│   └── example_valuations.ipynb
├── README.md
├── scripts
│   ├── random_quartic_reduction_experiments.py
│   └── test_tame_assumption.sage
├── semistable_model
│   ├── curves
│   │   ├── approximate_factors.py
│   │   ├── approximate_solutions.py
│   │   ├── component_graphs_of_plane_curves.py
│   │   ├── component_graphs.py
│   │   ├── cusp_resolution.py
│   │   ├── genus3_reduction_types.py
│   │   ├── integral_plane_curves.py
│   │   ├── plane_curves.py
│   │   ├── plane_curves_valued.py
│   │   └── stable_reduction_of_quartics.py
│   ├── finite_schemes.py
│   ├── geometry_utils
│   │   └── transformations.py
│   ├── stability
│   │   ├── admissible_functions.py
│   │   ├── bruhat_tits_building.py
│   │   ├── extension_search.py
│   │   ├── parametric_optimization.py
│   │   └── stability_function.py
│   └── valuations
│       └── linear_valuations.py
├── setup.py
└── spherical_stability_function.py
```

## Prerequisites

This project depends on **SageMath** and the **MCLF** library. You must install it before using this package.

  ```bash
  sage -pip install git+https://github.com/MCLF/mclf
  ```

## Installation

Choose the option that best fits your needs.

### Option A: Install as a Library (Usage Only)
**Best if:** You just want to use the `semistable_model` package in your own scripts and do not need to run the provided notebooks or edit the source code.

You can install the package directly from GitHub without cloning the full repository:
  ```bash
  sage -pip install git+https://github.com/kst3rn/StabilityFunction
  ```
(Note: You do not need to configure `nbstripout` for this method.)

### Option B: Clone for Development & Notebooks
**Best if:** You want to run the interactive notebooks or contribute to the code.

#### 1. Clone the Repository
  ```bash
  git clone https://github.com/kst3rn/StabilityFunction.git
  cd StabilityFunction
  ```

#### 2. Configure Git Filters (Important!)
This repository uses a `.gitattributes` file to automatically strip outputs from Jupyter Notebooks to keep the git history clean. **Before committing changes to notebooks, configure `nbstripout` locally so that notebook outputs are removed automatically from commits.**

First, install the `nbstripout` tool:
  ```bash
  pipx install nbstripout
  ```

To **activate the filter**, run the following command once inside the **repository root**:
  ```bash
  nbstripout --install
  ```
(This sets up the necessary filter definitions in your local `.git/config`)
  
#### 3. Install the Package in Editable Mode

Run the following command from the repository root to install the project in "editable" mode. This allows you to edit the source and see changes immediately without reinstalling.
  ```bash
  sage -pip install -e .
  ```

## Usage (Running the Notebooks)
To run the notebooks, you must start the Jupyter server from the repository root directory to ensure file paths resolve correctly.

  1. **Start Jupyter via Sage:**
  ```bash
  sage -n jupyter
  ```

  2. **Open a Notebook:** In the browser file tree, click on the `notebooks` folder and open the desired `.ipynb` file.
