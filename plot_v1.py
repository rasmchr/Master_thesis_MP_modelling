import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = pd.read_csv("Binned_mean_mean_concentration.csv")

import matplotlib.pyplot as plt

data_water = data[data['Type'] == 'Water']
data_bottom = data[data["Type"] == 'Bottom']
data_sediment = data[data["Type"] == 'Sediment']

fig, (ax_water, ax_sed) = plt.subplots(
    2, 1,            # 2 rows, 1 column
    sharex=False,    # different x axes
    sharey=False,    # different y axes
    figsize=(8, 10),
    gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.075}  # minimal space
)

# --- Water column plot ---
ax_water.scatter(
    data_water["Mean MP particles pr. litre"],
    data_water["Bin center (m)"],
    color="blue"
)


lower_w = data_water["Min value pr. litre"]

upper_w = data_water["Max value pr. litre"]


xerr_water = np.vstack([lower_w, upper_w])



ax_water.errorbar(data_water["Mean MP particles pr. litre"],
                  data_water["Bin center (m)"],
                  xerr = xerr_water,
                  fmt = 'o',
                  color="#63ace5",
                  label = "Water column samples"
                  )

ax_water.axhspan(0, 49.5, color='#63ace5', alpha=0.2, zorder=0, label = 'Water column')
ax_water.axhspan(49.5, 50.9, color='tab:brown', alpha=0.8, zorder=0, label = 'Seabed')
ax_water.scatter(
    data_bottom["Mean MP particles pr. litre"],
    len(data_bottom["Bin center (m)"])*[50],
    color="gold"
    )

ax_water.errorbar(data_bottom["Mean MP particles pr. litre"],
                  50,
                  xerr = np.vstack([data_bottom["Min value pr. litre"], data_bottom["Max value pr. litre"]]),
                  fmt = 'o',
                  color="#ffd3b6",
                  label = "At/Near seabed samples"
                  )


#ax_water.set_title("Water column")
ax_water.set_ylabel("Depth (m)")
ax_water.set_ylim(0, 50.9)

ax_water.invert_yaxis()   # optional for depth plots



# --- Sediment plot ---
ax_sed.scatter(
    data_sediment["Mean MP particles pr. litre"],
    data_sediment["Bin center (m)"] + 50,
    color="green"
    )

lower_s = data_sediment["Min value pr. litre"]

upper_s = data_sediment["Max value pr. litre"]


xerr_s = np.vstack([lower_s, upper_s])

ax_sed.errorbar(data_sediment["Mean MP particles pr. litre"],
                  data_sediment["Bin center (m)"]  + 50,
                  xerr = xerr_s,
                  fmt = 'o',
                  color="#7bc043",
                  label = 'Sediment samples'
                  )


ax_sed.axhspan(49.975, 50.001, color='tab:brown', alpha=0.8, zorder=0, label = 'Seabed')
ax_sed.axhspan(50.001, 51, color='tab:brown', alpha=0.1, zorder=0, label = 'Sediment')
ax_sed.set_ylim(49.975, 50.65)
#ax_sed.set_title("Sediment")
ax_sed.set_ylabel("Depth (m)")
ax_sed.set_xlabel("MP concentration (particles/L)")
ax_sed.invert_yaxis()

ax_sed.grid()
ax_water.grid()
ax_sed.legend(loc = 'right')
ax_water.legend(loc = 'right')



plt.show()


for size in data_sediment:
    fig, (ax_water, ax_sed) = plt.subplots(
        2, 1,            # 2 rows, 1 column
        sharex=False,    # different x axes
        sharey=False,    # different y axes
        figsize=(8, 10),
        gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.075}  # minimal space
    )

    # --- Water column plot ---
    ax_water.scatter(
        data_water["Mean MP particles pr. litre"],
        data_water["Bin center (m)"],
        color="blue"
    )


    lower_w = data_water["Min value pr. litre"]

    upper_w = data_water["Max value pr. litre"]


    xerr_water = np.vstack([lower_w, upper_w])



    ax_water.errorbar(data_water["Mean MP particles pr. litre"],
                    data_water["Bin center (m)"],
                    xerr = xerr_water,
                    fmt = 'o',
                    color="#63ace5",
                    label = "Water column samples"
                    )

    ax_water.axhspan(0, 49.5, color='#63ace5', alpha=0.2, zorder=0, label = 'Water column')
    ax_water.axhspan(49.5, 50.9, color='tab:brown', alpha=0.8, zorder=0, label = 'Seabed')
    ax_water.scatter(
        data_bottom["Mean MP particles pr. litre"],
        len(data_bottom["Bin center (m)"])*[50],
        color="gold"
        )

    ax_water.errorbar(data_bottom["Mean MP particles pr. litre"],
                    50,
                    xerr = np.vstack([data_bottom["Min value pr. litre"], data_bottom["Max value pr. litre"]]),
                    fmt = 'o',
                    color="#ffd3b6",
                    label = "At/Near seabed samples"
                    )





    # --- Sediment plot ---
    ax_sed.scatter(
        data_sediment["Mean MP particles pr. litre"],
        data_sediment["Bin center (m)"] + 50,
        color="green"
        )

    lower_s = data_sediment["Min value pr. litre"]

    upper_s = data_sediment["Max value pr. litre"]


    xerr_s = np.vstack([lower_s, upper_s])

    ax_sed.errorbar(data_sediment["Mean MP particles pr. litre"],
                    data_sediment["Bin center (m)"]  + 50,
                    xerr = xerr_s,
                    fmt = 'o',
                    color="#7bc043",
                    label = 'Sediment samples'
                    )

#ax_water.set_title("Water column")
ax_water.set_ylabel("Depth (m)")
ax_water.set_ylim(0, 50.9)

ax_water.invert_yaxis()   # optional for depth plots

ax_sed.axhspan(49.975, 50.001, color='tab:brown', alpha=0.8, zorder=0, label = 'Seabed')
ax_sed.axhspan(50.001, 51, color='tab:brown', alpha=0.1, zorder=0, label = 'Sediment')
ax_sed.set_ylim(49.975, 50.65)
#ax_sed.set_title("Sediment")
ax_sed.set_ylabel("Depth (m)")
ax_sed.set_xlabel("MP concentration (particles/L)")
ax_sed.invert_yaxis()

ax_sed.grid()
ax_water.grid()
ax_sed.legend(loc = 'right')
ax_water.legend(loc = 'right')