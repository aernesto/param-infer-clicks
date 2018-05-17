import pickle
import matplotlib.pyplot as plt

results = pickle.load(open("data/mse.pkl", "rb"))

# first file
file = results[0]['file']
trial_numbers = file[1]
# lin 2 lin fit
data = results[0]['stats']['lin']

# MSE plot
plt.subplot(1, 2, 1)
plt.plot(trial_numbers, [x[0] for x in data])
plt.xlabel('trial nb')
plt.ylabel('MSE')

# Width plot
plt.subplot(1, 2, 2)
plt.plot(trial_numbers, [x[1] for x in data])
plt.xlabel('trial nb')
plt.ylabel('width')

plt.show()

