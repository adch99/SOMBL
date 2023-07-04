import numpy as np
import matplotlib.pyplot as plt

def block(x):
    start = -12
    stop = 12
    base = 0
    high = 10
    return base + high*(np.heaviside(x-start, 0) - np.heaviside(x-stop, 0))

def gaussian(x):
    sigma = 1.5
    return np.exp(-x**2 / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)

x_block = np.linspace(-20, 20, 1000)
y_block = block(x_block)
x_gauss = np.linspace(-4, 4, 100)
y_gauss = gaussian(x_gauss)

x = np.linspace(-16, 16, 1099)
y = np.convolve(y_gauss, y_block) / x_gauss.shape[0]

mask_full = np.logical_or(x > -14, x < 14)
mask_middle = np.logical_and(mask_full, np.logical_and(x > -7, x < 7))
mask_left = np.logical_and(mask_full, x <= -7)
mask_right = np.logical_and(mask_full, x >= 7)

x_offset = 15

fig = plt.figure(figsize=(5, 4))
# print(y)
# plt.plot(x_block, y_block)
# plt.plot(x_gauss, y_gauss)
plt.plot(x+x_offset, y, color="black", linewidth=0.5)
# plt.fill_between(x[mask_full]+x_offset, y[mask_full], y2=0, color="#ff7f0e")
plt.fill_between(x[mask_middle]+x_offset, y[mask_middle], y2=0, color="#1f77b4")
plt.fill_between(x[mask_left]+x_offset, y[mask_left], y2=0, color="#ff7f0e")
plt.fill_between(x[mask_right]+x_offset, y[mask_right], y2=0, color="#ff7f0e")

plt.text(x_offset, 0.6, "Delocalized\nStates", fontsize="xx-large",
        ha="center", va="center", color="white")
# plt.text(x_offset, 0.6, "Localized\nStates", fontsize="xx-large",
#         ha="center", va="center", color="white")

# plt.text(14+x_offset, 0.7, , fontsize="xx-large",
#         ha="center", va="center", color="#ff7f0e")


prop = {    
        "arrowstyle": "-|>,head_width=0.4,head_length=0.8",
        "shrinkA":0,
        "shrinkB":0,
        "facecolor": "#ff7f0e"
}
plt.annotate("Localized\nStates", (9+x_offset, 0.4), (15+x_offset,0.75),
            arrowprops=prop, fontsize="xx-large",
            ha="center", va="center", color="#ff7f0e")

# plt.arrow(11+x_offset, 0.55, -1, -0.1, width=0.05)

plt.xlabel(r"$E$", fontsize="xx-large", rotation=0)
plt.ylabel(r"$\nu(E)$", fontsize="xx-large", rotation=0, labelpad=20)

plt.yticks([])
plt.xticks([])

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)

ax.spines[["left", "bottom"]].set_position(("data", 0))


# plt.ylim(0, 13)

plt.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)
plt.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)

fname = "dos_schematic_deloc"
# fname = "dos_schematic_loc"
plt.savefig("plots/PDFs/" + fname + ".pdf")
plt.savefig("plots/PNGs/" + fname + ".png")
plt.show()