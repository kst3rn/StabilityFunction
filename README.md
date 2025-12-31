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

This project is designed to be installed as a local Python package within SageMath.

### 1. Clone the Repository
  ```bash
  git clone https://github.com/kst3rn/StabilityFunction.git
  cd StabilityFunction
  ```

### 2. Install the Package

Run the following command from the repository root to install the project in "editable" mode. This allows you to edit the source and see changes immediately without reinstalling.
  ```bash
  sage -pip install -e .
  ```

## Usage (Running the Notebooks)
To run the notebooks, you must start the Jupyter server from the repository root directory.

### 1. Navigate to the repository root:
  ```bash
  cd StabilityFunction
  ```

### 2. Start Jupyter:
  ```bash
  sage -n jupyter
  ```

### 3. In the browser file tree click on the notebooks folder and open the desired notebook.

## Development Setup (Handling Notebooks)

To keep the git history clean and the repository size small, this project uses **automated notebook stripping**. This removes output cells and execution counts from Jupyter Notebooks before they are committed.

**Requirement:**
Because the repository includes a `.gitattributes` file enforcing this filter, **you must configure `nbstripout` locally**, otherwise Git might report errors about a missing filter driver.

### 1. Install nbstripout
You can install it globally (recommended) or within your Sage/Python environment.

**Option A: System-wide (Recommended)**
Use `pipx` to install it isolated from your system packages:
  ```bash
  pipx install nbstripout
  ```

### 1. Activate the filter
Even if installed, you need to tell Git how to use it for this repository. Run this once inside the repository root:
  ```bash
  nbstripout --install
  ```
(This sets up the necessary filter definitions in your local `.git/config`)
