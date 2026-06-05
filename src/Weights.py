import matplotlib.pyplot as plt
import numpy as np


L = 150.0

water_raw = [
    (2.5,   3.5667891682785307),
    (7.5,   2.6969926716060204),
    (12.5,  0.9080153525061379),
    (17.5,  0.5380893150797977),
    (27.5,  0.15000000000000002),
    (32.5,  0.17777777777777778),
    (42.5,  0.07500000000000001),
    (47.5,  0.1),
]
sediment_raw = [   # (bin_centre_cm, concentration)
    (0.0125,               634.7414656833154),
    (0.037500000000000006, 478.3568408220815),
    (0.0625,               970.1680149695655),
    (0.08750000000000001,  2343.5408338909124),
    (0.1125,               1420.2905042168202),
    (0.1375,               1237.2282190132125),
    (0.16250000000000003,  1026.581643277907),
    (0.1875,               299.62609220263096),
    (0.21250000000000002,  681.4776292062617),
    (0.2375,               458.20696863375105),
    (0.2625,               422.2530388157952),
    (0.3125,               640.1515738798312),
    (0.3875,               63.55031042845192),
    (0.41250000000000003,  187.55853823063677),
    (0.48750000000000004,  30.846954308574094),
    (0.5375000000000001,   30.516169700062356),
    (0.6125,               17.81769496827207),
]

water_z   = np.array([r[0]                for r in water_raw])
water_val = np.array([r[1]                for r in water_raw])
sed_z     = np.array([r[0] + 50.0 for r in sediment_raw])  # cm→m + offset
sed_val   = np.array([r[1]                for r in sediment_raw])

# Zero-padding: penalise mass past the last observed depth
pad_z   = np.arange(sed_z[-1] + 0.5, L + 0.5, 10)
pad_val = np.zeros(len(pad_z))

obs_z   = np.concatenate([water_z, sed_z, pad_z])
obs_val = np.concatenate([water_val, sed_val, pad_val])
# --- fixed solver inputs ---
t_final = 70.0
t_array = np.linspace(0, t_final, 2000)
z_grid  = np.linspace(0, 150, 2000)



sigma_obs    = np.where(obs_val > 0, 0.5 * obs_val, 1.0)
weights_real = 1.0 / sigma_obs**2
plt.plot(weights_real, obs_z, 'o-')
plt.grid()
plt.gca().invert_yaxis()
plt.xlabel("Weight")
plt.ylabel("Depth (m)")
plt.show()

weights_real = obs_val / np.max(obs_val)
plt.axhline(y=50, linestyle='--', color='gray', label='seabed')
plt.plot(weights_real, obs_z, 'o-', color='tab:blue')
plt.grid()
plt.gca().invert_yaxis()
plt.xlabel("Weight")
plt.ylabel("Depth (m)")
plt.legend()
plt.show()


other_weights = np.concatenate([
    water_val / np.max(water_val),
    sed_val   / np.max(sed_val),
    np.ones_like(pad_z)*0.5
])

other_weights = other_weights + 0.5
other_weights = other_weights/np.max(other_weights)

plt.axhline(y=50, linestyle='--', color='gray', label='seabed')
plt.plot(other_weights, obs_z, '--', color='black', alpha = 0.5)

plt.plot(other_weights, obs_z, 'o', color='tab:orange')
plt.grid()
plt.gca().invert_yaxis()
plt.xlabel("Weight")
plt.ylabel("Depth (m)")
plt.legend()
plt.show()

print(np.min(other_weights), np.max(other_weights))
