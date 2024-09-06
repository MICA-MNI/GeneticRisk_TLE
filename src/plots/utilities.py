# Helper functions
def plot_ctx(map, cmap, range=None, save=False, filename=None):
    map = parcel_to_surface(map, 'aparc_conte69');
    plot_cortical(array_name=map, surface_name=conte69',
                  color_bar=True, color_range=range, cmap=cmap, size=(1000,500),
                  screenshot=save, filename=filename, transparent_bg=True, scale=(5,5))

def plot_ctx_sig(map, cmap, thresh=0.05, range=None, save=False, filename=None):
    sig = np.ones(map.shape)
    sig[map>=thresh] = NaN
    sig = parcel_to_surface(sig, 'aparc_conte69', fill=NaN);
    plot_cortical(array_name=x, surface_name=conte69',
                  color_bar=True, color_range=range, cmap=cmap, nan_color=(0,0,0,1),
                  size=(1000,500), screenshot=save, filename=filename, transparent_bg=True, scale=(5,5))

def plot_sctx(map, cmap, range=None, save=False, filename=None):
    plot_subcortical(array_name=map, ventricles=False,
                  color_bar=True, color_range=range, cmap=cmap, size=(1000,500),
                  screenshot=save, filename=filename, transparent_bg=True, scale=(5,5))

def plot_sctx_sig(map, cmap, thresh=0.05, range=None, save=False, filename=None):
    sig = np.ones(map.shape)
    sig[map>=thresh] = NaN
    plot_subcortical(array_name=sig, ventricles=False,
                  color_bar=True, color_range=range, cmap=cmap, size=(1000,500),
                  screenshot=save, filename=filename, transparent_bg=True, scale=(5,5))

def circular_barchart(data, ):

    VALUES = df["variance"].values
    LABELS = df["roi"].values
    GROUP = df["group"].values
    FDR = df["fdr"].values

    COLOURS = ['#2F5597' if x < 0.05 else '#A1A1A1' for x in fdr]

    ANGLES = np.concatenate((np.linspace(rad(7.5), rad(172.5), num=6, endpoint=False), np.linspace(rad(187.5), rad(352.5), num=6, endpoint=False)))
    WIDTH = rad(300) / len(ANGLES)
    GROUPS_SIZE = [len(i[1]) for i in df.groupby("group")]

    OFFSET = (ANGLES[4] + ANGLES[3]) / 2# (rad(7.5) + rad(172.5))/2 + rad(14)

    fig, ax = plt.subplots(figsize=(25, 25), subplot_kw={"projection": "polar"})

    ax.set_theta_offset(OFFSET)
    ax.set_ylim(-0.075, 0.25)
    ax.set_frame_on(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2])

    ax.bar(
    ANGLES, VALUES, width=WIDTH, color=COLOURS,
    linewidth=2, zorder=10)

    # Add line below bars
    x1 = np.linspace(ANGLES[0], ANGLES[5], num=1000000)
    ax.plot(x1, [-0.02] * 1000000, linewidth=3.5, color="k")
    ax.text(
    np.mean(x1), -0.05, 'L', color="#333333", fontsize=14,
    fontweight="bold", ha="center", va="center"
    )
    x2 = np.linspace(ANGLES[0]-(WIDTH/2 + 0.05), ANGLES[5]+(WIDTH/2 + 0.05), num=1000000)
    ax.plot(x2, [0] * 1000000, color="#bebebe", lw=2)
    ax.plot(x2, [0.05] * 1000000, color="#bebebe", lw=2)
    ax.plot(x2, [0.1] * 1000000, color="#bebebe", lw=2)
    ax.plot(x2, [0.15] * 1000000, color="#bebebe", lw=2)
    ax.plot(x2, [0.2] * 1000000, color="#bebebe", lw=2)

    x1 = np.linspace(ANGLES[6], ANGLES[11], num=1000000)
    ax.plot(x1, [-0.02] * 1000000, linewidth=3.5, color="k")
    ax.text(
    np.mean(x1), -0.05, 'R', color="#333333", fontsize=14,
    fontweight="bold", ha="center", va="center"
    )
    x2 = np.linspace(ANGLES[6]-(WIDTH/2 + 0.05), ANGLES[11]+(WIDTH/2 + 0.05), num=1000000)
    ax.plot(x2, [0] * 1000000, color="#bebebe", lw=2)
    ax.plot(x2, [0.05] * 1000000, color="#bebebe", lw=2)
    ax.plot(x2, [0.1] * 1000000, color="#bebebe", lw=2)
    ax.plot(x2, [0.15] * 1000000, color="#bebebe", lw=2)
    ax.plot(x2, [0.2] * 1000000, color="#bebebe", lw=2)

def spatial_correlation(data):

def null_raincloud(data):

def plot_matrix(data):

def barchart(data):

def null_boxplot(data):
