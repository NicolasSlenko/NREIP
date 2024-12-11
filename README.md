# MATLAB Setup Instructions

To run the MATLAB scripts provided in this repository, follow the instructions below to set up MATLAB and install any necessary package dependencies.

## Prerequisites

Ensure you have MATLAB installed on your system. You can download MATLAB from the [MathWorks website](https://www.mathworks.com/products/matlab.html).

## Installation

1. **Download and Install MATLAB:**
    - Visit the [MathWorks website](https://www.mathworks.com/products/matlab.html).
    - Follow the instructions to download and install MATLAB on your system.

2. **Verify Installation:**
    - Open MATLAB.
    - Type `ver` in the Command Window and press Enter. This will display the installed MATLAB version and toolboxes.

## Package Dependencies

Some scripts may require additional MATLAB toolboxes or packages. Follow the steps below to install any necessary dependencies:

1. **Check Required Toolboxes:**
    - Each script will list the required toolboxes at the beginning of the file or in the documentation.

2. **Install Toolboxes:**
    - Open MATLAB.
    - Go to the "Home" tab and click on "Add-Ons" -> "Get Add-Ons".
    - Search for the required toolbox and click "Install".

3. **Install Additional Packages:**
    - If additional packages are required, they will be listed in the script documentation.
    - Use the MATLAB Add-On Explorer to search for and install these packages.

## Running Scripts

Once MATLAB and the necessary dependencies are installed, you can run the scripts by:

1. Opening MATLAB.
2. Navigating to the directory containing the script.
3. Typing the script name (without the `.m` extension) in the Command Window and pressing Enter.

Example:
```matlab
cd('/path/to/your/scripts')
exampleScript
```

For any issues or further assistance, refer to the MATLAB documentation or contact the repository maintainer. 


# Reinforcement Learning Model Training with Blade Element Momentum Theory

This script trains a reinforcement learning (RL) model using the Blade Element Momentum (BEM) theory.

## Instructions for Running the Scripts

1. Ensure all files are in the same directory.
2. Open the `RL_Rough_Structure5.slx` Simulink model in MATLAB.
3. Run the `testRL.m` script to start the RL training process.
Note: All of these files are found in the Original MatLab Scripts folder 

## File Descriptions

- **`nacaCoordinates_new.m`**:
  - Generates the airfoil coordinates and plots the airfoil shape.
  - Outputs these coordinates to the `xfoil.m` script.

- **`xfoil.m`**:
  - Uses the airfoil coordinates to calculate lift and drag coefficients.
  
- **`bem.m`**:
  - Performs Blade Element Momentum (BEM) theory calculations.
  - Outputs the efficiency of the system.

- **`testRL.m`**:
  - Contains the configuration for the RL model.
  - Trains the model for 2000 episodes using a DDPG (Deep Deterministic Policy Gradient) agent. The number of episodes can be adjusted in this script.

- **`RL_Rough_Structure5.slx`**:
  - A Simulink model that contains the "Next State" and "Reward" functions:
    - **Next State Function**: Runs the physics model to generate the output.
    - **Reward Function**: Calculates the reward based on the RL agent's action.

## Future Improvements

- **Physics Model**:
  - The current physics model occasionally produces efficiencies greater than 1. A more accurate physics model or physics-informed neural networks could be explored as alternatives.

- **Reward Function**:
  - The reward function should be optimized to suit the new physics model, based on a chosen efficiency threshold.

