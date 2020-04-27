detpop_color = 'C0'


def plot_detpop_mean_quantities(ax, rp_list, x, qty):
    y = []
    for rp in rp_list:
        if qty == 'en':
            y.append(rp.mean_en_exp())
        elif qty == 'r_haff':
            y.append(rp.r_haff_exp())
    ax.plot(x, y, c=detpop_color, label='theory', marker='.')
