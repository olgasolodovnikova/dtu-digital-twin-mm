# Power goal
power_goal = [0.8, 1.3, 2.0, 1.5 ,1.8, 2.3, 2.8, 2.5, 2.7, 1.5] # [MW]

for i in range(len(power_goal)):
    power_goal[i] = power_goal[i]*1e5 #to W
