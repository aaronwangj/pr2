from random import seed
import itertools
import copy
import sys
from collections import defaultdict


def kmersFunction(sequence, k):
    
    kmerIndices = defaultdict(list)

    for i in range(len(sequence) - k + 1):
        kmerIndices[str(sequence[i:i+k])].append(i)
    return kmerIndices

def hammingDistance(s1, s2):
    how_many_differences = 0

    if len(s1) != len(s2):
        return -1
    else:
        for i in range(len(s1)):
            if s1[i] != s2[i]:
                how_many_differences += 1
    return how_many_differences

'''
calculating absolute freq of kmers in all reads (1 dist)
iterate through every kmer with freq below t

iterate through kmers, find distance between each below t with every above t
sort by distance, then by absolute freq of the matched kmer
'''
class Corrector:
    def __init__(self):
        self.freqMap = None

    def correctionFunc(self, reads_path, k, t, d, isFile):

        try:
            #for kmers
            k = int(k)
            #frequency threshold
            t = int(t)
            #how many positions we can change
            d = int(d)
        except:
            print('error: bad values of k,t,d')
            return
        if k <= 0 or t < 0 or d < 0:
            print('bad input')
            return


        if isFile:
            try:
                with open(reads_path, 'r') as file:
                    all_reads = file.read().split('\n')
            except:
                print('error: bad reads path')
                return
        else:
            all_reads = copy.deepcopy(reads_path)

        try:
            if not self.freqMap:
                kmerIndices = kmersFunction(''.join(all_reads), k)
                self.freqMap = {k:len(v) for (k,v) in kmerIndices.items()}
        

            corrected_set = set()

            for kmer in self.freqMap:

                '''
                (kmer, distance, absolute frequency)
                '''
                kmerReplacements = []


                if self.freqMap[kmer] < t:
                    for kmer2 in self.freqMap:
                        if kmer != kmer2 and self.freqMap[kmer2] >= t:
                            distance = hammingDistance(kmer, kmer2) 
                            if distance <= d:
                                kmerReplacements.append((kmer2, distance, self.freqMap[kmer2]))
                
                if kmerReplacements:
                    '''
                    sort by:
                    ascending distance
                    descending freq
                    '''
                    kmerReplacements.sort(key = lambda x: (x[1], -x[2]))
                    # print('hi2')
                
                    replacer = kmerReplacements[0][0]

                    for index, read in enumerate(all_reads):
                        kmerIndicesRead = kmersFunction(read, k)
                        for readKmer in kmerIndicesRead:
                            if readKmer == kmer:
                                for starting_position in kmerIndicesRead[readKmer]:
                                    # print('readkmer: ', readKmer)
                                    # print('replacer: ', replacer)
                                    # print('hi')
                                    # print(len(all_reads[index]))
                                    # print(len(replacer))
                                    all_reads[index] = read[:starting_position] + replacer + read[starting_position+len(replacer):]
                                    # print(len(all_reads[index]))
                                    # print('bye')
                                corrected_set.add(index)
            
            return sorted(list(corrected_set)), all_reads
        except:
            print('error')
            return


def main():
    try:
        sys1 = sys.argv[1]
        sys2 = sys.argv[2]
        sys3 = sys.argv[3]
        sys4 = sys.argv[4]
    except:
        print('error: bad number of arguments')
        return
    try:
        c = Corrector()
        indices, all_reads, _ = c.correctionFunc(sys1, sys2, sys3, sys4, True)

        print(*indices, sep = ",")
        print('--------------------')
        for read in all_reads:
            print(read)
    except:
        print('')
        return

if __name__ == "__main__":
    main()