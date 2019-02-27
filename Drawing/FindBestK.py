import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
from scipy.spatial.distance import cosine
import subprocess
import shutil
import copy

def main():
    """
    bic_dic = load_bic()
    if(os.path.isdir(result_path + "@figure") == False):
        os.makedirs(result_path + "@figure")
    for i in range(2, max_k):
        if(i < 10): str_i = "0" + str(i)
        else: str_i = str(i)
        spectra_file = result_path[:-1] + "_" + str(bic_dic[i][0]) +\
                "/result_k" + str_i + ".txt"
        cmd = "cp " + spectra_file + " " + result_path
        subprocess.call(cmd.split())
    draw_bic(bic_dic)
    best_k, best_iter = find_best_k(bic_dic)
    """
    best_k = 6
    if(best_k < 10): str_best_k = "0" + str(best_k)
    else: str_best_k = str(best_k)
    spectra_file = result_path + "result_k" +\
            str_best_k + ".txt"
    spectra = []
    lines = open(spectra_file, "r").readlines()
    for i,line in enumerate(lines[2:2+best_k]):
        temp_list = line.split(" ")
        spectra.append([])
        for j,temp in enumerate(temp_list):
            if(j != len(temp_list)-1):spectra[i].append(float(temp))
    spectra = align_spectra(spectra)
    draw_spectra(spectra)
    match_list = make_match(spectra)
    draw_match(spectra, match_list)
    #for i in range(1,51):
    #    shutil.rmtree(result_path[:-1] + "_" + str(i))

def load_bic():
    bic_dic = {}
    for i in range(1,51):
        for j in range(2,max_k):
            if(j < 10): str_k = "0" + str(j)
            else: str_k = str(j)
            bic_file = result_path[:-1] + "_" + str(i) + "/result_k" + str_k + ".txt"
            if(os.path.exists(bic_file)):
                lines = open(bic_file, "r").readlines()       
                temp_bic = float(lines[0].split()[0])
                if((j not in bic_dic.keys()) or (bic_dic[j][1] <\
                    temp_bic)):
                    bic_dic.update({j:[i, temp_bic]})
    return bic_dic

def draw_bic(bic_dic):
    left = list(); height = list()
    for key in bic_dic.keys():
        left.append(key)
        height.append(bic_dic[key][1])
    fig = plt.figure()
    plt.bar(left, height, align="center")
    plt.title("The transition of BIC")
    plt.xlabel("The number of mutational process")
    plt.ylabel("The value of BIC")
    plt.xlim(0, max_k)
    fig.savefig(result_path + "@figure/bic.png", dpi=300)
    plt.close(1)

def find_best_k(bic_dic):
    temp = [0, -1e100]
    for key in bic_dic.keys():
        if(bic_dic[key][1] > temp[1]):
            temp = [key, bic_dic[key][1]]
    best_iter = bic_dic[temp[0]][0]
    best_k = temp[0]
    return best_k, best_iter

def draw_spectra(spectra):
    labels, colorlist = make_labels_and_colors()
    for i in range(len(spectra)):
        fig = plt.figure(figsize=(6,2))
        left = np.arange(1,97,1)
        height_ex = spectra[i]
        title = "Predicted Signature " + str(i+1)
        ax = fig.add_subplot(111)
        ax.bar(left, height_ex, width=1, color=colorlist, align="center")
        height_limits = []
        for j in range(len(height_ex)):
            height_limits.append(height_ex[j])
        max_height_limit = max(height_limits)
        upper_lim = math.ceil(max_height_limit*10)/10
        ax.set_ylim(0, upper_lim)
        ax.set_xticks(left)
        ax.set_xticklabels(labels)
        for tick in ax.get_xticklabels():
            tick.set_rotation(90)
        ax.tick_params(labelsize=4)
        ax.set_xlabel('mutation x', fontsize=7)
        ax.set_ylabel('p (mutation = x)', fontsize=7)
        ax.set_title(title, fontsize=7)
        fig.tight_layout()
        name = result_path + '@figure/predicted_' + str(i+1) + '.png'
        fig.savefig(name, dpi=200)
        plt.close(1)

def make_labels_and_colors():
    labels = 96*[0]
    for i in range(96):
        first = i % 16
        if(first == 0 or first == 1 or first == 2 or first == 3):
            labels[i] = 'A'
        if(first == 4 or first == 5 or first == 6 or first == 7):
            labels[i] = 'C'
        if(first == 8 or first == 9 or first == 10 or first == 11):
            labels[i] = 'G'
        if(first == 12 or first == 13 or first == 14 or first == 15):
            labels[i] = 'T'
        for j in range(16):
            if(i == j):
                labels[i] += '(C>A)'
            if(i == j+16):
                labels[i] += '(C>G)'
            if(i == j+32):
                labels[i] += '(C>T)'
            if(i == j+48):
                labels[i] += '(T>A)'
            if(i == j+64):
                labels[i] += '(T>C)'
            if(i == j+80):
                labels[i] += '(T>G)'
        second = i % 4
        if(second == 0):
            labels[i] += 'A'
        if(second == 1):
            labels[i] += 'C'
        if(second == 2):
            labels[i] += 'G'
        if(second == 3):
            labels[i] += 'T'

    colorlist = 96*[0]
    for i in range(16):
        colorlist[i] = 'r'
    for i in range(16, 32):
        colorlist[i] = 'g'
    for i in range(32, 48):
        colorlist[i] = 'b'
    for i in range(48, 64):
        colorlist[i] = 'c'
    for i in range(64, 80):
        colorlist[i] = 'm'
    for i in range(80, 96):
        colorlist[i] = 'y'
    
    return labels, colorlist

def make_match(spectra):
    K = len(spectra)
    known_spectra = load_knowns()
    JS = np.zeros([K, len(known_spectra)])
    for i in range(K):
        for j in range(len(known_spectra)):
            JS[i,j] = cosine(spectra[i], known_spectra[j])
    min_JS = np.zeros([K]); min_sig = np.zeros([K])
    for i in range(K):
        min_j = 0
        for j in range(1,30):
            if(JS[i,j] < JS[i,min_j]):
                min_j = j
        min_JS[i] = JS[i,min_j]
        min_sig[i] = min_j
    return list(min_sig)

def align_spectra(spectra):
    new_spectra = copy.deepcopy(spectra)
    i = 0
    while (i < 96):
        j = int(i/4); count = 0
        while(count != 4):
            for k in range(len(spectra)):
                for l in range(4):
                    new_spectra[k][i+4*count+l] = spectra[k][j+l]
            j += 24; count += 1
        i += 16
    spectra = copy.deepcopy(new_spectra)
    return spectra

def load_knowns():
    known_spectra = np.zeros([30,96])
    lines = open("data/signature_probability.txt", "r").readlines()
    for i,line in enumerate(lines):
        if(i != 0):
            temp_list = line.split()
            for j in range(33):
                if (j >= 3):
                    known_spectra[j-3, i-1] = float(temp_list[j])
    new_known_spectra = known_spectra.copy()
    i = 0
    while (i < 96):
        j = int(i/4); count = 0
        while(count != 4):
            for k in range(30):
                for l in range(4):
                    new_known_spectra[k,i+4*count+l] = known_spectra[k,j+l]
            j += 24; count += 1
        i += 16
    known_spectra = new_known_spectra.copy()
    return known_spectra

def draw_match(spectra, match_list):
    known_spectra = load_knowns()
    labels,colorlist = make_labels_and_colors()
    for i in range(len(spectra)):
        fig = plt.figure()
        left = np.arange(1,97,1)
        heights = [spectra[i]]
        heights.append(known_spectra[int(match_list[i])])
        title1 = "Predicted Signture " + str(i+1)
        title2 = "COSMIC Known Signature " + str(int(match_list[i]+1))
        for j in range(2): fig.add_subplot(2,1,j+1)
        axs = plt.gcf().get_axes()
        for j,ax in enumerate(axs):
            ax.bar(left, heights[j], width=1, color=colorlist, align="center")
            ax.set_xticks(left)
            ax.set_xticklabels(labels)
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)
            ax.tick_params(labelsize=5)
            ax.set_ylabel("p (mutation = x)")
            if(j == 0): ax.set_title(title1)
            elif(j == 1): ax.set_title(title2)
        fig.tight_layout()
        name = result_path + "@figure/match_" + str(i+1) + ".png"
        fig.savefig(name, dpi=300)
        plt.close(1)

if __name__ == "__main__":
    args = sys.argv
    mut_file = args[1]
    result_path = "result/" + mut_file + "/"
    max_k = 21
    #cmd = "scp -r oil:project/PLSA/result/" + mut_file + "_* ./result/"
    #subprocess.call(cmd.split())
    main()
