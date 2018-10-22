import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_all_3d():
    def plot_3D(x, y, z, name, ax):
        planet = name.split(".")[0]
        ax.plot(x, y, z, label=planet)
        ax.set_xlabel(r'$x$ [AU]', size=11)
        ax.set_ylabel(r'$y$ [AU]', size=11)
        ax.set_zlabel(r'$z$ [AU]', size=11)
        ax.legend(fontsize=9, loc="best")

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')
    ax.set_zlim(-.5, .5)
    # ax.set_xlim(-30, 30)
    # ax.set_ylim(-30, 30)

    files = ["Earth.dat","Jupiter.dat","Mercury.dat","Saturn.dat","Uranus.dat","Mars.dat","Neptune.dat","Sun.dat","Venus.dat"]
    # files = ["Mercury.dat","Venus.dat","Earth.dat","Mars.dat","Sun.dat"]
    for fl in files:
        data = file(fl)
        x = []; y = []; z = []
        for line in data:
            k, l, m = line.split()
            x.append(float(k))
            y.append(float(l))
            z.append(float(m))
        plot_3D(x, y, z, fl, ax)
    # plt.savefig("../figs/all_planets_zoom.png")
    plt.show()


def plot_2D(filename, label, title=None):
    infile = file(filename)
    x = []
    y = []
    for line in infile:
        l, r = line.split()
        x.append(float(l))
        y.append(float(r))

    plt.plot(x, y, label=label)
    plt.xlabel(r"$x$ [AU]", size=11); plt.ylabel(r"$y$ [AU]", size=11)
    plt.axis("equal")
    if title != None:
        plt.title(title)
    plt.legend(fontsize=12)
    plt.grid()
    # plt.savefig("../figs/X.png")
    plt.show()

def plot_errors(filename, plt_name, title=None):
    infile = file(filename)
    log_h = []
    log_rel_error = []
    for line in infile:
        l, r = line.split(",")
        log_h.append(float(l))
        log_rel_error.append(float(r))
    plt.plot(log_h, log_rel_error)
    plt.scatter(log_h, log_rel_error)
    plt.grid()
    plt.xlabel(r"$\log _{10} (h)$", size=12)
    plt.ylabel(r"$\log _{10} ( \epsilon )$", size=12)
    if title != None:
        plt.title(title)
    dir = "../figs/" + plt_name
    plt.savefig(dir)
    plt.show()


# plot_2D(filename="optimal_verlet.dat", label="Earth orbit", title="Verlet solver")
# plot_2D(filename="euler.dat", label="Earth orbit", title="Euler solver")

# plot_errors(filename="errors_euler.txt", plt_name="errors_euler.png", title="Euler solver")
# plot_errors(filename="errors_verlet.txt", plt_name="errors_verlet.png", title="Verlet solver")

# plot_2D(filename="escape.dat", label=r"$v_0 = 8.88$AU/yr")
# plot_2D(filename="tmp.dat", label=r"$v_0 = 8.86$AU/yr")
# plt.grid()
# plt.savefig("../figs/escape.png")
# plt.show()

# filenames = ["beta3.000000.dat"]
# for filename in filenames:
#     label = filename.split(".")
#     tmp = label[0] + label[1][0]
#     plot_2D(filename, r"$\beta = 3$")
# plt.legend(fontsize=12)
# plt.grid()
# plt.show()

# plot_2D("Earth.dat", "Earth")
# plot_2D("Jupiter.dat", "Jupiter")
# plt.legend(fontsize=12)
# plt.grid()
# plt.savefig("../figs/earth_1000.png")
# plt.show()


# plot_2D("Earth.dat", "Earth")
# plot_2D("Jupiter.dat", "Jupiter")
# plot_2D("Sun.dat", "Sun")
# plt.legend(fontsize=12)
# plt.grid()
# plt.savefig("../figs/cm_system.png")
# plt.show()

# plot_all_3d()

# plot_2D("merc_GR.dat", "Mercury orbit")
