
import matplotlib.pyplot as plt
from pycalphad import Database, calculate, variables as v
from pycalphad.plot.utils import phase_legend
import numpy as np

# Load database and choose the phases that will be plotted
db_nbre = Database('Cu-Ni.tdb')
my_phases_nbre = ['LIQUID','FCC_A1']

# Get the colors that map phase names to colors in the legend
legend_handles, color_dict = phase_legend(my_phases_nbre)

fig = plt.figure(figsize=(9,6))
ax = fig.gca()

# Loop over phases, calculate the Gibbs energy, and scatter plot GM vs. X(RE)
for phase_name in my_phases_nbre:
    result = calculate(db_nbre, ['CU', 'NI','CU%','VA%','NI%'], phase_name, P=101325, T=1800, output='GM')
    ax.scatter(result.X.sel(component='NI'), result.GM, marker='.', s=5, color=color_dict[phase_name])

# Format the plot
ax.set_xlabel('X(RE)')
ax.set_ylabel('GM')
ax.set_xlim((0, 1))
ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.6))
plt.show()