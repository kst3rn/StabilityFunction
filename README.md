# StabilityFunction

...

## Project Structure

* **semistable_model/**: The main Python package containing the mathematical logic.
* **notebooks/**: Interactive Jupyter notebooks demonstrating the usage of the project.
* **documentation/**: Contains supplementary documentation.

## Prerequisites

This project depends on the **MCLF** (MacLane Valuations and Berkovich Spaces) library. You must install it before using this package.

```bash
sage -pip install git+[https://github.com/MCLF/mclf](https://github.com/MCLF/mclf)
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

## For Contributors: Keeping Git Clean
Jupyter notebooks often clutter git history with execution counts and output data. To keep the repository clean, we highly recommend using `nbstripout`. This tool automatically removes outputs before you commit.

## Setup (Run once)
```bash
# 1. Install nbstripout
sage -pip install nbstripout

# 2. Configure the git filter for this repository
nbstripout --install
```
Once installed, you can saves your notebooks with outputs visible on your screen, but `git` will only commit the clean code.
