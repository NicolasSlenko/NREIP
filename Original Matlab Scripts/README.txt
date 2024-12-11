***READ ME***

The files in this script are used to train the reinforcement learning model using the blade element momentum
theorem.

To run these scripts, ensure all files are in the same path and run the testRL.m script with RL_Rough_Structure5.slx open.

The nacaCoordinates_new.m script generates the airfoil coordinates and plots the shape of the airfoil. These coordinates are passed
to the xfoil.m script to find lift and drag coefficients, which are used by bem.m to carry out the blade
element momentum theory calculations and output the efficiency. 

The testRL.m script contains the information for the RL model in conjunction with the Simulink file. The model is set to train for 2000
episodes using a DDPG agent, which can be adjusted in the testRL.m file. The Simulink file contains the Next State and Reward functions, which run
the physics model to generate an output and calculate the reward for the RL agent's action, resepectively. 

***Future Improvements***

- The current physics model occasionally outputs efficiencies greater than 1, so another physics model or physics-informed
  neural networks could be utilized
- The reward function will need to be optimized for this new physics model based on a chosen efficiency threshold 

