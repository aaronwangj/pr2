import random
import itertools
import copy
import sys
import os
import numpy as np
import correction
from collections import defaultdict
from collections import Counter
import matplotlib.pyplot as plt
random.seed(42)

class Solver:
    def __init__(self):
        self.globalFreq = None

    def unitaryScore(self, s1, s2):
        score = 0

        if len(s1) != len(s2):
            print('UNITARY SCORING: strings not same length')
            return
        else:
            for i in range(len(s1)):
                if s1[i] == s2[i]:
                    score += 1
                else:
                    score -= 1
        return score

    def applicationFunc(self, hiv_path, matrix_path, k, t):
        d = 2
        k = int(k)
        t = int(t)

        '''
        1-e^{-N\cdot \frac{50}{9181}}=.999
        '''
        N = 1268 #how many times we sample reads

        '''
        open files
        '''
        try:
            with open(hiv_path, 'r') as file:
                genome = file.read().replace('\n', '')

        except:
            print('error: bad path')
            return

        '''
        randomKmers: sampled reads from genome
        randomKmersMutated: list of reads post-mutation
        ifMutated: if each particular kmer was mutated
        '''

        if not os.path.exists('true_reads.txt'):
            originalReads = []
            mutatedReads = []
            for _ in range(N):
                readStart = random.randint(0, len(genome) - 50)
                randRead = genome[readStart:readStart+50]
                originalReads.append(randRead)

                possibleRead = (randRead + '.')[:-1]
                for index, base in enumerate(possibleRead):
                    if random.random() < .01:
                        new_base = base
                        while base == new_base:
                            new_base = random.choice(['A','C','G','T'])
                        possibleRead = possibleRead[:index] + new_base + possibleRead[index+1:]
                mutatedReads.append(possibleRead)

            with open('true_reads.txt', 'w') as file:
                for element in originalReads:
                    file.write(element + ",")

            with open('mutated_reads.txt', 'w') as file:
                    for element in mutatedReads:
                        file.write(element + ",")
        else:
            with open('true_reads.txt', 'r') as file:
                originalReads = file.read().split(',')
            with open('mutated_reads.txt', 'r') as file:
                mutatedReads = file.read().split(',')

        _, correctedReads= correction.correctionFunc(mutatedReads, k, t, 2, False)
        
        S_m = 0
        S_k = 0
        for i in range(N):
            S_m += self.unitaryScore(originalReads[i], mutatedReads[i]) / N
            S_k += self.unitaryScore(originalReads[i], correctedReads[i]) / N

        print(correctedReads == mutatedReads)
        print(correctedReads == originalReads)
        print(mutatedReads == originalReads)
        print('SK: ', S_k)
        print('SM: ', S_m)

        return -np.log((50 - S_k) / (50 - S_m))

def main():
    sys1 = sys.argv[1]
    sys2 = sys.argv[2]

    T = [4, 6, 8, 10, 12]

    all_t_scores = []
    solver = Solver()

    for t in T:
        t_score_list = []
        for k in range(6, 26):
            print(t, k)
            res_score = solver.applicationFunc(sys1, sys2, k, t)
            t_score_list.append(res_score)
            
        all_t_scores.append(t_score_list.copy())

    colors = ["red","green","blue","yellow","pink","black","orange","purple","beige","brown","gray","cyan","magenta"]
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for i in range(len(all_t_scores)):

        ax1.scatter([j for j in range(6,26)], all_t_scores[i], s=10, c=colors[i], marker="o", label=4 + 2*i)
    plt.legend(loc='upper left');
    plt.savefig('t_figure.png')
    plt.show()

    t_star = 4
    t_low = t_star - 2
    t_high = t_star + 2

    mid = []
    low = []
    high = []

    for k in range(6, 26):
        mid.append(applicationFunc(sys1, sys2, k, t_star))
        low.append(applicationFunc(sys1, sys2, k, t_low))
        high.append(applicationFunc(sys1, sys2, k, t_high))

    all_new_scores = [low, mid, high]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for i in range(len(all_new_scores)):
        if i == 0:
            l = 'low'
        elif i == 1:
            l = 'mid'
        elif i == 2:
            l = 'high'
        else:
            print('uhoh')
            l = -1

        ax1.plot([j for j in range(6,26)], all_new_scores[i], c=colors[i], marker="o", label=l)
    plt.legend(loc='upper left');
    plt.savefig('t_figure2.png')
    plt.show()
    
    print(all_t_scores)

if __name__ == "__main__":
    main()