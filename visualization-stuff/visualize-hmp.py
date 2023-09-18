import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

SRS_ids = ['011247', '013258', '014524', '014530', '015128', '015381', '015688', '016665', '016683',
           '017851', '018691', '019372', '019794', '044001', '076809', '077646', '1041036',
           '1041092', '1041093', '1041115', '1041126', '1054691', '143014', '148146']


def extract(pref, suff):
    lst = []
    short = []
    ct = 0
    for id in SRS_ids:
        file = open(pref + id + suff, 'r')
        lines = file.readlines()
        for line in lines:
            if line[0] == 'A' or line[0] == 'T' or line[0] == 'C' or line[0] == 'G':
                if len(line) - 1 >= 80:
                   #  lst.append(len(line) - 1)
                    if len(line) - 1 <= 300:
                        short.append(len(line) - 1)
                else:
                    ct += 1
        file.close()
    return lst, ct, short


def create_fig(cuttle, ggcat, k, freq, under, cut, ggc):
    max_length = max(len(cuttle), len(ggcat))
    cuttle_padded = cuttle + [None] * (max_length - len(cuttle))
    ggcat_padded = ggcat + [None] * (max_length - len(ggcat))

    data = {'Cuttle\nct:' + str(len(cuttle)//1000000) + 'M rem:' + str(cut//1000000) + 'M': cuttle_padded, 'GGCAT\nct:' + str(len(ggcat)//1000000) + 'M rem:' + str(ggc//1000000) + 'M': ggcat_padded}

    df = pd.DataFrame(data)

    # Create the violin plot
    sns.violinplot(data=df, scale="count")

    # Set plot labels and title
    plt.xlabel('Tools')
    plt.ylabel('Unitig Distribution')

    if under:
        plt.title('Cumulative SRS data plot, k = ' + k + ', freq = ' + freq + ' unitigs shorter than 300')
        plt.savefig('/home/Users/mt102/data/plots/SRS' + '-u300-k' + k + '-f' + freq + '.png')
    else:
        plt.title('Cumulative SRS data plot, k = ' + k + ', freq = ' + freq)
        plt.savefig('/home/Users/mt102/data/plots/SRS' + '-k' + k + '-f' + freq + '.png')
    # Close the plot
    plt.close()


def violins(k, freq):
    cuttle, cuttle_extra, cuttle_short = extract("/home/Users/mt102/shared/cuttle-res/HMP-SRS/SRS",
                                                 "." + str(k) + ".f" + str(freq) + ".fa")
    ggcat, ggcat_extra, ggcat_short = extract("/home/Users/mt102/shared/ggcat-res/HMP-SRS/SRS",
                                              "." + str(k) + ".f" + str(freq))

    # create_fig(cuttle, ggcat, k, freq, False, cuttle_extra, ggcat_extra)
    create_fig(cuttle_short, ggcat_short, k, freq, True, cuttle_extra, ggcat_extra)

    # Inform the user about the saved file
    print("Done with k = " + k + " and freq = " + freq + ".")


freqs = ['1', '2']

lens = ['21', '25', '31', '35', '41', '45', '51']

for k in lens:
    for freq in freqs:
        violins(k, freq)

