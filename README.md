# StabilityFunction

...

## Project Structure

* **semistable_model/**: The main Python package containing the mathematical logic.
* **notebooks/**: Interactive Jupyter notebooks demonstrating the usage of the project.
* **documentation/**: Contains supplementary documentation.

## Prerequisites

This project depends on the **MCLF** library. You must install it before using this package.

  ```bash
  sage -pip install git+https://github.com/MCLF/mclf
  ```

## Installation

Choose the option that best fits your needs.

### Option A: Install as a Library (Usage Only)
**Best if:** You just want to use the `semistable_model` in your owm scripts and do not need to run the provided notebooks or edit the source code.

You can install the package directly from GitHub without cloning the full repository:
  ```bash
  sage -pip intall git+https://github.com/kst3rn/StabilityFunction
  ```
(Note: You do not need to configure `nbstripout` for this method.)

### Option B: Clone for Development & Notebooks
**Best if:** You want to run the interactive notebooks or contribute to the code.

### 1. Clone the Repository
  ```bash
  git clone https://github.com/kst3rn/StabilityFunction.git
  cd StabilityFunction
  ```

### 2. Configure Git Filters (Important!)
This repository uses a `.gitattributes` file to automatically strip outputs from Jupyter Notebooks to keep the git history clean. **You must configure this locally, or Git may report errors.**

First, install the `nbstripout` tool:
  ```bash
  pipx install nbstripout
  ```

To **activate the filter**, run the following command once inside the **repository root**:
  ```bash
  nbstripout --install
  ```
(This sets up the necessary filter definitions in your local `.git/config`)
  
### 3. Install the Package in Editable Mode

Run the following command from the repository root to install the project in "editable" mode. This allows you to edit the source and see changes immediately without reinstalling.
  ```bash
  sage -pip install -e .
  ```

## Usage (Running the Notebooks)
To run the notebooks, you must start the Jupyter server from the repository root directory to ensure file paths resolve correctly.

  1. **Navigate to the repository root:**
  ```bash
  cd StabilityFunction
  ```

  2. **Start Jupyter via Sage:**
  ```bash
  sage -n jupyter
  ```

  3. **Open a Notebook:** In the browser file tree, click on the `notebooks` folder and open the desired `ipynb` file.
