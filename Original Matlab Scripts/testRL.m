obsInfo = rlNumericSpec([3, 1], LowerLimit = 0, UpperLimit = 1);
actInfo = rlNumericSpec([3, 1], LowerLimit = -0.5 , UpperLimit = 0.5);

env = rlSimulinkEnv ("RL_Rough_Structure5", "RL_Rough_Structure5/RL Agent", obsInfo, actInfo);

opts = rlDDPGAgentOptions ("DiscountFactor", 0.6);

agent = rlDDPGAgent(obsInfo,actInfo, opts);

opts = rlTrainingOptions("MaxEpisodes", 2000,"MaxStepsPerEpisode", 10, "ScoreAveragingWindowLength", 100, "StopTrainingCriteria", "AverageReward", "StopTrainingValue", 10000);

train(agent, env, opts);