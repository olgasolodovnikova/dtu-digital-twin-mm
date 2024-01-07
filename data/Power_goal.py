# Power goal
power_goal = [0.2, 0.3, 0.5, 1.0 ,2.0, 3.0, 3.5, 3.7, 4.0, 4.1] # [MW]

for i in range(len(power_goal)):
    power_goal[i] = power_goal[i]*1e5 #to W
