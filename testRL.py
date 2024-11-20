# Define observation and action specifications
obsInfo = rlNumericSpec([3, 1], LowerLimit=0, UpperLimit=1)
actInfo = rlNumericSpec([3, 1], LowerLimit=-0.5, UpperLimit=0.5)

# Create the environment
env = rlSimulinkEnv("RL_Rough_Structure5", "RL_Rough_Structure5/RL Agent", obsInfo, actInfo)

# Configure the agent options
opts = rlDDPGAgentOptions(DiscountFactor=0.6)

# Create the DDPG agent
agent = rlDDPGAgent(obsInfo, actInfo, opts)

# Set training options
opts = rlTrainingOptions(
    MaxEpisodes=2000,
    MaxStepsPerEpisode=10,
    ScoreAveragingWindowLength=100,
    StopTrainingCriteria="AverageReward",
    StopTrainingValue=10000
)

# Train the agent
train(agent, env, opts)
