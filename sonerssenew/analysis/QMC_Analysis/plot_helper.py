import matplotlib as mpl
import pandas as pd
from matplotlib.cm import get_cmap



def plot_data(df, axv, gby, grouby_param, xye):
    #df = pd.read_csv(file, delim_whitespace=True, comment='#')
    font_size = 28
    col = iter(["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#0072B2"])
    grouped = df.groupby(gby)    
    cmap = get_cmap('viridis', len(grouped))
    #col = iter(["Red"])
    if len(xye) > 2:
        for (keyy, group), color in zip(grouped, cmap.colors):
            key = int(keyy) if float(keyy).is_integer() else keyy
            if keyy in grouby_param:
                axv.errorbar(group[xye[0]], group[xye[1]], yerr=group[xye[2]], c=color, label=key,
                             linestyle='-', marker='o', capthick=4, markerfacecolor='none', elinewidth=1,
                             markersize=8, capsize=4, markeredgewidth=1.5, linewidth=1.5)
    else:
        for (keyy, group), color in zip(grouped, cmap.colors):
            key = int(keyy) if float(keyy).is_integer() else keyy
            if keyy in grouby_param:
                axv.errorbar(group[xye[0]], group[xye[1]], c=color, label=key,
                             linestyle='-', marker='o', capthick=4, markerfacecolor='none', elinewidth=1,
                             markersize=8, capsize=4, markeredgewidth=1.5, linewidth=1.5)                

    for axis in [axv.xaxis, axv.yaxis]:
        axis.set_tick_params(which='major', length=5, width=1.5, labelsize=font_size)
        axis.set_tick_params(which='minor', length=2.5, width=1.2, labelsize=font_size)
        axis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4))

    axv.minorticks_on()

    #axv.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: f"{x:.0f}"))
    #axv.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: f"{x:.3f}"))
    
    # Adjust the width of the axes
    axv.spines['top'].set_linewidth(1.4)    # Top border
    axv.spines['bottom'].set_linewidth(1.4) # Bottom border
    axv.spines['left'].set_linewidth(1.4)   # Left border
    axv.spines['right'].set_linewidth(1.4)  # Right border    
    
