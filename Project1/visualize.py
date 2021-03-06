import matplotlib.pyplot as plt

def read_result_and_plot(filename):
    x = []
    v = []
    u = []
    file = open(filename)

    for line in file:
        row = line.split(",")
        x.append(float(row[0]))
        u.append(float(row[1]))
        v.append(float(row[2]))
    file.close()

    plt.figure()
    # plt.plot(x,u, label=r"$u(x)$")
    plt.plot(x,v, label=r"$n={:.1e}$".format(len(x)-2))
    plt.xlabel(r"$x$"); plt.ylabel(r"$v(x)$")
    plt.legend()
    plt.savefig("./figs/v_x.png")
    plt.show()

def read_alternative_file_format():
    # Used for reading the armadillo vectors..
    x = []
    v = []
    file1 = open("x.txt")
    file2 = open("v.txt")

    for line in file1:
        x.append(float(line))
    for line in file2:
        v.append(float(line))

    plt.plot(x,v)
    plt.show()

    file1.close()
    file2.close()

def read_error_and_plot(filename):
    log_h = []; log_eps = []
    file = open(filename)
    for line in file:
        row = line.split(",")
        log_h.append(float(row[0]))
        log_eps.append(float(row[1]))
    file.close()

    plt.figure()
    plt.scatter(log_h,log_eps)
    plt.xlabel(r"$\log_{10} (h)$")
    plt.ylabel(r"$\epsilon $")
    plt.grid()
    plt.savefig("./figs/eps_h.png")

    import numpy as np
    best_h = log_h[ np.argmin(log_eps) ]
    msg = "The optimal value for log(h) is {} based on available data.".format((best_h))
    print msg

    plt.show()


def main():
    """ Example runs:
        python visualize.py simulation1.txt 0
        python visualize.py errors.txt 1
    """

    import sys
    if len(sys.argv) <= 2:
        msg = "Run script with the filename containing the results or the errors. Second argument should be '0' for results or '1' for errors."
        print msg
        sys.exit()
    else:
        filename = sys.argv[1]
        if sys.argv[2] == '0':
            read_result_and_plot(filename)
        elif sys.argv[2] == '1':
            read_error_and_plot(filename)


if __name__ == "__main__":
    main()
