import warnings
import csv
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.optimize import curve_fit
import seaborn as sns
warnings.filterwarnings("ignore")
import pandas as pd

## This script fits a polynomal function and fits the L, y0, x0 and k parameters.
## This version of the script discards the data if the curve goes down after a certain number of x values.
## This is in order to exclude the waviness of the data.

# the model used to fit the data
def model(x, L, k, x0, y0):
    y =  y0 + (L/(1 + np.exp(-k * (x - x0))))
    return y

# calculating the error of the model
def Rsquared(y, yfit):
   S = (np.sum([(yfit_i - y_i) ** 2 for (y_i, yfit_i) in zip(y, yfit)])/(len(y)-2))**0.5
   return S/abs(np.mean(y))*100

# variables
x_name = 'Time'                       # x variabel
y_name = 'nodeDirectionAv'              # y variabel
input = "Direction_Elongation_main1.csv"# input file
max_r_squared = 80                      # maximum % error of the fit
min_x_salt_cut = 500                   # salt curve will not be cut beneath this x value
min_x_control_cut = 200               # control curve will not be cut beneath this x value
min_diff_salt_cut = 0.4                 # when the difference between two subsequenct y values is less then this the
min_dif_control_cut = 1.2               # curve will be cut

# creating a folder
file_result_path = "Plots/"
if not os.path.exists(file_result_path):
    os.makedirs(file_result_path)
file_result_path2 = "Parameters/"
if not os.path.exists(file_result_path2):
    os.makedirs(file_result_path2)

# reading in the input file
with open(input) as csvfile:
    reader = csv.DictReader(csvfile, dialect='excel')
    datalist = []
    datalist = list(reader)

# creating additional variables, empty arrays and finding unique variables
parameters = []
treatments = sorted(set([item['Treatment'] for item in datalist]))
genotypes = sorted(set([item['Genotype'] for item in datalist]))
blues = ['#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#084594']
greens = ['#a1d99b', '#74c476', '#41ab5d', '#238b45', '#005a32']

# looping over genotypes and treatments and getting x and y values
for genotype in genotypes:
    for treatment in treatments:
        cntr = 0
        roots = sorted(set([item['Root'] for item in datalist if item['Genotype'] == genotype and item['Treatment'] == treatment ]))

        # selecting data set specific to treatment, genotype and root
        for root in roots:
            x_temp = [float(row[x_name]) for row in datalist if row['Genotype'] == genotype and row['Treatment'] == treatment and row['Root'] == root]
            y_temp = [float(row[y_name]) for row in datalist if row['Genotype'] == genotype and row['Treatment'] == treatment and row['Root'] == root]
            x_temp = [item - min(x_temp) for item in x_temp]
            x_temp2 = []
            y_temp2 = []
            x = []
            y = []

            #getting out double y_values
            for item in range (0, len(y_temp)):
                if y_temp[item] not in y_temp2:
                    y_temp2.append(y_temp[item])
                    x_temp2.append(x_temp[item])

            #cutting of curves to discard waviness
            for item in range(0, len(x_temp2)):
                if treatment == 'NaCl100mM' and (y_temp2[item] > y_temp2[item - 1] + min_diff_salt_cut
                                                   or x_temp2[item] < min_x_salt_cut):
                    x.append(x_temp2[item])
                    y.append(y_temp2[item])
                elif treatment == 'Control' and (y_temp2[item] >= y_temp2[item - 1] + min_dif_control_cut
                                                 or x_temp2[item] < min_x_control_cut):
                    x.append(x_temp2[item])
                    y.append(y_temp2[item])
                else:
                    break

            #fitting the model to curves
            fit = []
            try:
                p, e = curve_fit(model, x, y, p0=(93, 0.005, 300, -90), method='trf', bounds=([0, 0, 0, -np.inf],  np.inf))
                xd = np.linspace(0.001, max(x), 100)
                fit = model(xd, *p)
                xy_fit = model(x, *p)
                r_squared_fit = Rsquared(y, xy_fit)
            except Exception, err:
                print 'fail: \n genotype: ', genotype, '\t root: ', root

            #plotting your fit and adding parameters to dictionary
            if len(fit)>0 and r_squared_fit < max_r_squared:
                y_beginning = model(min(x), *p)
                increase = p[0] + p[3]  - y_beginning
                parameters.append({'genotype': genotype, 'treatment': treatment, 'root': root, 'Rsquared': r_squared_fit,'L': p[0], 'k': p[1], 'x0': p[2], 'y0': p[3], 'total_turning': increase})
                if treatment == 'Control':
                    plt.scatter(x, y, label=treatment + ' root: ' + root + ' Rsquared:' + str(round(r_squared_fit, 1)) + ' L:' + str(round(p[0])) + ' k:' + str(round(p[1], 4)) +
                        ' x0: ' + str(round(p[2], 2)) + ' y0: ' + str(round(p[3], 2)) + ' total_turning: ' + str(round(increase)), color= greens[cntr], edgecolors = 'black', linewidth='0.1')
                if treatment == 'NaCl100mM':
                    plt.scatter(x, y, label=treatment + ' root: ' + root + ' Rsquared:' + str(round(r_squared_fit, 1)) + ' L:' + str(round(p[0])) + ' k:' + str(round(p[1], 4)) +
                        ' x0: ' + str(round(p[2], 2)) + ' y0: ' + str(round(p[3], 2)) + ' total_turning: ' + str(round(increase)), color=blues[cntr], edgecolors = 'black', linewidth='0.1')
                plt.plot(xd, fit, color='grey')
            else:
                if len(fit)> 0:
                    print 'fail_max_r_squared: \n genotype: ', genotype, '\t root: ', root
                if treatment == 'Control':
                    plt.scatter(x,y,  label=treatment, color= greens[cntr], edgecolors = 'black', linewidth='0.1')
                else:
                    plt.scatter(x, y, label=treatment, color=blues[cntr], edgecolors = 'black', linewidth='0.1')
            cntr+=1

    #plot variables to make it look pretty ;)
    plt.ylabel(y_name)
    plt.xlabel(x_name)
    lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1))
    title = 'control, genotype: ' + genotype
    plt.title(title)
    plt.savefig(file_result_path + genotype.replace("/", "")  + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()

#writing raw parameters to csv
keys = ['genotype', 'treatment','root','Rsquared', 'L', 'k','x0', 'y0', 'total_turning']
with open('raw_variables.csv', 'wb') as output_file:
    dict_writer = csv.DictWriter(output_file, keys)
    dict_writer.writeheader()
    dict_writer.writerows(parameters)

# calculating averages for salt and control and the ratio
differences = []
for genotype in genotypes:
    if len([item for item in parameters if item["genotype"] == genotype]) > 1:
        total_turning_salt = np.mean([item['total_turning'] for item in parameters if item["genotype"] == genotype and item["treatment"] == 'NaCl100mM'])
        total_turning_cntrl = np.mean([item['total_turning'] for item in parameters if item["genotype"] == genotype and item["treatment"] == 'Control'])
        k_salt = np.mean([item["k"] for item in parameters if item["genotype"] == genotype and item["treatment"] == 'NaCl100mM'])
        k_cntrl = np.mean([item["k"] for item in parameters if item["genotype"] == genotype and item["treatment"] == 'Control'])
        x0_salt = np.mean([item["x0"] for item in parameters if item["genotype"] == genotype and item["treatment"] == 'NaCl100mM'])
        x0_cntrl = np.mean([item["x0"] for item in parameters if item["genotype"] == genotype and item["treatment"] == 'Control'])
        differences.append({'genotype': genotype, 'total_turning_salt':total_turning_salt, 'total_turning_control':total_turning_cntrl,
                           'total_turning':total_turning_salt/ total_turning_cntrl,'k_salt':k_salt, 'k_control':k_cntrl, 'k': k_salt/k_cntrl,
                            'x0_salt': x0_salt, 'x0_control': x0_cntrl, 'x0': x0_salt/x0_cntrl})

# writing ratios to csv
keys_dif = ['genotype', 'total_turning_salt', 'total_turning_control', 'total_turning', 'k_salt',
            'k_control', 'k', 'x0_salt', 'x0_control', 'x0']
with open('Salt_responses.csv', 'wb') as output_file:
    dict_writer = csv.DictWriter(output_file, keys_dif)
    dict_writer.writeheader()
    dict_writer.writerows(differences)

# making barplots of the variables
with open('raw_variables.csv') as csvfile:
    data = pd.read_csv(csvfile)
    columns = list(data.columns.values)[4:]

    for column in columns:
        sns.set_style("whitegrid")
        plt.figure()
        ax = sns.boxplot(x="genotype", y=column, hue="treatment", data=data, palette="Set3")
        ax.figure.savefig(file_result_path2 + column + ".png")




