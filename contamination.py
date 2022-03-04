from random import seed
import itertools
import sys
from collections import defaultdict


def kmersFunction(sequence, k):
    
    kmerIndices = defaultdict(list)

    for i in range(len(sequence) - k + 1):
        kmerIndices[str(sequence[i:i+k])].append(i)
    return kmerIndices

def contaminationFunc(reads_path, vector_path, k):

    try:
        k = int(k)
        assert k > 0
    except:
        print('error: bad value of k')
        return
    try:
        with open(reads_path, 'r') as file:
            all_reads = file.read().split('\n')
    except:
        print('error: bad reads path')
        return

    try:
        with open(vector_path, 'r') as file:
            vector = file.read().replace('\n', '')
    except:
        print('error: bad vector path')
        return

    try:
        kmerIndices = kmersFunction(vector, k)

        cleaned_reads = []
        contaminated_set = set()


        for index, read in enumerate(all_reads):

            left_cutoff = -float('inf')
            right_cutoff = float('inf')

            for kmer in kmerIndices:
                if read[:k] == kmer:
                    for vector_match_start in kmerIndices[kmer]:
                        vector_match_end = vector_match_start + k - 1
                        read_match_end = k - 1

                        while read_match_end < len(read) and vector_match_end < len(vector):
                            if read[read_match_end] == vector[vector_match_end]:
                                vector_match_end += 1
                                read_match_end += 1
                            else:
                                break
                        left_cutoff = max(left_cutoff, read_match_end)

                if read[-k:] == kmer:
                    for vector_match_start in kmerIndices[kmer]:
                        read_match_end = len(read) - k - 1

                        while read_match_end > 0 and vector_match_start > 0:
                            if read[read_match_end] == vector[vector_match_start]:
                                vector_match_start -= 1
                                read_match_end -= 1
                            else:
                                break
                        right_cutoff = min(right_cutoff, read_match_end)

            if right_cutoff != float('inf'):
                read = read[:right_cutoff]
                contaminated_set.add(index)

            if left_cutoff != -float('inf'):
                read = read[left_cutoff+1:]
                contaminated_set.add(index)

            cleaned_reads.append(read)
        
        return cleaned_reads, sorted(list(contaminated_set))
    except:
        print('error')
        return

def main():
    try:
        sys1 = sys.argv[1]
        sys2 = sys.argv[2]
        sys3 = sys.argv[3]
    except:
        print('error: bad number of arguments')
        return
    try:

        reads, contaminated = contaminationFunc(sys1, sys2, sys3)

        print(*contaminated, sep = ",")
        print('--------------------')
        for read in reads:
            print(read)
    except:
        print('')




if __name__ == "__main__":
    main()