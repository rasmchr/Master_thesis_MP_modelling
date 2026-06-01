import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = pd.read_excel("Master data_v01.xlsx", sheet_name="Sheet1")

# Size distribution plot at depth no more than 50m
# Subset the data 
data_water = data[(data['Sample type'] == 'W') & (data['Water depth at location'] <= 50)]
data_bottom = data[(data['Sample type'] == 'B') & (data['Water depth at location'] <= 50)]
data_sediment = data[(data['Sample type'] == 'S') & (data['Water depth at location'] <= 200)]


# Now calculate the distribution of the sizes
Sample_type = ["Water column", "At/near bottom", "Sediment column"]

Samples = {
    'Mp count [0-100μm[': np.array([
        data_water['Mp count [0-100μm['].sum() / data_water['Transformed data (MP particles pr. Litre)'].sum(),
        data_bottom['Mp count [0-100μm['].sum() / data_bottom['Transformed data (MP particles pr. Litre)'].sum(),
        data_sediment['Mp count [0-100μm['].sum() / data_sediment['Transformed data (MP particles pr. Litre)'].sum()
    ]),

    'Mp count [100μm-500μm[': np.array([
        data_water['Mp count [100μm-500μm['].sum() / data_water['Transformed data (MP particles pr. Litre)'].sum(),
        data_bottom['Mp count [100μm-500μm['].sum() / data_bottom['Transformed data (MP particles pr. Litre)'].sum(),
        data_sediment['Mp count [100μm-500μm['].sum() / data_sediment['Transformed data (MP particles pr. Litre)'].sum()
    ]),

    'Mp count [500μm-5mm[': np.array([
        data_water['Mp count [500μm-5mm['].sum() / data_water['Transformed data (MP particles pr. Litre)'].sum(),
        data_bottom['Mp count [500μm-5mm['].sum() / data_bottom['Transformed data (MP particles pr. Litre)'].sum(),
        data_sediment['Mp count [500μm-5mm['].sum() / data_sediment['Transformed data (MP particles pr. Litre)'].sum()
    ])
}

fig, ax = plt.subplots(figsize=(8,5))

x = np.arange(len(Sample_type))     # [0,1,2]
bottom = np.zeros(len(Sample_type))
colors = ['#a8e6cf', '#ffd3b6', '#ff8b94']  # Small, Medium, Large

for (label, values), color in zip(Samples.items(), colors):
    ax.bar(x, values, label=label, bottom=bottom, color=color)
    bottom += values

ax.set_xticks(x)
ax.set_xticklabels(Sample_type)
ax.set_ylabel("MP fraction")
ax.legend(title="Fraction size", loc='lower right')

plt.tight_layout()
plt.grid(axis='y')
plt.show()

# Create new column in each of the subsetted data sets with total counts in them
data_water['<500μm'] = data_water.loc[:,'Mp count [0-100μm['] + data_water.loc[:,'Mp count [100μm-500μm[']
data_bottom['<500μm'] = data_bottom.loc[:,'Mp count [0-100μm['] + data_bottom.loc[:,'Mp count [100μm-500μm[']
data_sediment['<500μm'] = data_sediment.loc[:,'Mp count [0-100μm['] + data_sediment.loc[:,'Mp count [100μm-500μm[']

# Now plot the water in 
#plt.scatter(data_water_our_size[''])

#plt.scatter(data_water['<500μm'], data_water['Transformed depth'])
#
#plt.scatter(data_bottom['<500μm'], [50]*len(data_bottom['<500μm']))
#
#plt.scatter(data_sediment['<500μm'], 50 + data_sediment['Transformed depth'])
#
#plt.grid()
#
#plt.gca().invert_yaxis()
#
#plt.show()


import matplotlib.pyplot as plt

fig, (ax_water, ax_sed) = plt.subplots(
    2, 1,            # 2 rows, 1 column
    sharex=False,    # different x axes
    sharey=False,    # different y axes
    figsize=(8, 10),
    gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.075}  # minimal space
)

# --- Water column plot ---
num_studies_water = len(data_water['Study'].unique())
ax_water.scatter(
    data_water["Transformed data (MP particles pr. Litre)"],
    data_water["Transformed depth"],
    color="blue",
    label = f"Water column MP samples from {num_studies_water} studies"
)


ax_water.axhspan(49.5, 50.9, color='tab:brown', alpha=0.8, zorder=0, label = 'Seabed')
num_studies_bottom = len(data_bottom['Study'].unique())
ax_water.scatter(
    data_bottom["Transformed data (MP particles pr. Litre)"],
    len(data_bottom["Transformed depth"])*[50],
    color="gold",
    label = f"At bottom MP samples from {num_studies_bottom} studies"
)
ax_water.set_title("Water column")
ax_water.set_ylabel("Depth (m)")
ax_water.set_ylim(0, 50.9)

ax_water.invert_yaxis()   # optional for depth plots


#xmin, xmax = ax_water.get_xlim()    # or xmin, xmax = -1, 5000



# Check the amount of studies
num_studies_sediment = len(data_sediment['Study'].unique())
# --- Sediment plot ---
ax_sed.scatter(
    data_sediment["Transformed data (MP particles pr. Litre)"],
    data_sediment["Transformed depth"] + 50,
    color="green",
    label = f"Sediment MP samples from {num_studies_sediment} studies"
)
ax_sed.axhspan(49.975, 50.001, color='tab:brown', alpha=0.8, zorder=0, label = 'Seabed')
ax_sed.set_ylim(49.975, 50.65)
#ax_sed.set_title("Sediment")
ax_sed.set_ylabel("Depth (m)")
ax_sed.set_xlabel("MP concentration (particles/L)")
ax_sed.invert_yaxis()

ax_sed.grid()
ax_water.grid()
ax_sed.legend()
ax_water.legend()



plt.show()
x = 5
sub_dat = data_water
bins_all = np.arange(0, sub_dat['Transformed depth'].max() + x, x)
bin_centers = (bins_all[:-1] + bins_all[1:]) / 2
all_mean_vals = []
all_std_vals = []
all_min_vals = []
all_max_vals = []

for i in range(len(bins_all) - 1):
    # step where we subset the data


    bin_dat = sub_dat.loc[
        (sub_dat['Transformed depth'] >= bins_all[i]) &
        (sub_dat['Transformed depth'] < bins_all[i + 1])
    ]
    #print(bin_dat)

    all_mean_vals.append(bin_dat['Transformed data (MP particles pr. Litre)'].mean())
    all_std_vals.append(bin_dat['Transformed data (MP particles pr. Litre)'].std())
    all_min_vals.append(bin_dat['Transformed data (MP particles pr. Litre)'].min())
    all_max_vals.append(bin_dat['Transformed data (MP particles pr. Litre)'].max())

    all_data_counts = len(bin_dat['Transformed data (MP particles pr. Litre)'])
    


rows = []

# --- Water rows ---
for bc, mean, std, min, max in zip(bin_centers, all_mean_vals, all_std_vals, all_min_vals, all_max_vals):
    rows.append({
        'Type': 'Water',
        'Bin center (m)': bc,
        'Mean MP particles pr. litre': mean,
        'Std MP particles pr. litre': std,
        'Min value pr. litre': min,
        'Max value pr. litre': max,
        'Data count': all_data_counts
        
    })

# --- Bottom row ---
rows.append({
    'Type': 'Bottom',
    'Bin center (m)': None,
    'Mean MP particles pr. litre': data_bottom['Transformed data (MP particles pr. Litre)'].mean(),
    'Std MP particles pr. litre': data_bottom['Transformed data (MP particles pr. Litre)'].std(),
    'Min value pr. litre': data_bottom['Transformed data (MP particles pr. Litre)'].min(),
    'Max value pr. litre': data_bottom['Transformed data (MP particles pr. Litre)'].max(),
    'Data count': len(data_bottom['Transformed data (MP particles pr. Litre)'])
})

x = 0.025
sediment_data = data_sediment
bins_all = np.arange(0, sediment_data['Transformed depth'].max() + x, x)
bin_centers = (bins_all[:-1] + bins_all[1:]) / 2
all_mean_vals = []
all_std_vals = []
all_min_vals = []
all_max_vals = []

for i in range(len(bins_all) - 1):
    # step where we subset the data


    bin_dat = sediment_data.loc[
        (sediment_data['Transformed depth'] >= bins_all[i]) &
        (sediment_data['Transformed depth'] < bins_all[i + 1])
    ]

    all_mean_vals.append(bin_dat['Transformed data (MP particles pr. Litre)'].mean())
    all_std_vals.append(bin_dat['Transformed data (MP particles pr. Litre)'].std())
    all_min_vals.append(bin_dat['Transformed data (MP particles pr. Litre)'].min())
    all_max_vals.append(bin_dat['Transformed data (MP particles pr. Litre)'].max())

    all_data_counts = len(bin_dat)

# now plot the mean and std in a histogram plot with error bars
plt.figure(figsize = (10, 10))
plt.scatter(all_mean_vals, bin_centers)
plt.errorbar(all_mean_vals, bin_centers, xerr = all_std_vals, fmt= "o")
plt.xlabel("mean MP particles pr litre")
plt.ylabel("Depth (cm)")
plt.title(f"No bound on the max depth\nSediment samples\nNumber of samples: {len(sediment_data)}")
plt.gca().invert_yaxis()
plt.grid()
plt.show()


# --- Sediment rows ---
for bc, mean, std, min, max in zip(bin_centers, all_mean_vals, all_std_vals, all_min_vals, all_max_vals):
    rows.append({
        'Type': 'Sediment',
        'Bin center (m)': bc,
        'Mean MP particles pr. litre': mean,
        'Std MP particles pr. litre': std,
        'Min value pr. litre': min,
        'Max value pr. litre': max,
        'Data count': all_data_counts
    })

output_df = pd.DataFrame(rows)
output_df.to_csv("Binned_mean_mean_concentration.csv")
