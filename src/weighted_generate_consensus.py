#!/usr/bin/env python3
import argparse
import os

dict_content = {}
list_correctness=[]
list_last_values=[]



def add_to_dict(key, val, dict):

    if key in dict:
        # append the new number to the existing array at this slot
        #print("adding to dict " , dict.get(key))
        dict[key] = dict.get(key) + val
    else:
        # create a new array in this slot
        #print(dict_content.get(key))
        #print(dict_content)
        dict[key] = val
        if dict == dict_content:
            list_correctness.append([])
            list_last_values.append([])
        #print("added new key ", key)

def read_file (path):

    file = open(path, "r")
    count = -1

    for line in file:

        if line[0] != ">":
            #print("line strip  - ", line.strip("\n"))
            #print(count)
            add_to_dict(count, line.strip("\n"), dict_content)
            #print(dict_content)
        else:
            count += 1



    #print(dict_content)



    file.close()


def update_correctness(chosen, k):
    count = 0


    for i in list_correctness:
        #print(list_last_values[count], count, chosen)
        #print(list_correctness[count])
        if len(list_correctness[count]) == k:
            list_correctness[count].pop(0)

        if list_last_values[count][len(list_last_values[count]) -1] == chosen:
            list_correctness[count].append(1)
            #print(1)
        else:
            list_correctness[count].append(0)
            #print(0)

        count += 1



def generate_consensus (output, k):

    count = 0
    finished = 0

    list_bases = "actguACTGU"
    consensus = []
    keys_finished = []



    while finished < len(dict_content):

        dict_bases = {}
        for key in dict_content:
            #print(key, finished)

            try:
                base = dict_content.get(key)[count]
            except:
                if key not in keys_finished:
                    keys_finished.append(key)
                    finished += 1


            if base in list_bases: #check if it is one of the bases
                #print(list_last_values[key], key, base)
                if base == "a" or base == "A":
                    add_to_dict("A", sum(list_correctness[key]), dict_bases)
                    if len(list_last_values[key]) == k:
                        list_last_values[key].pop(0)
                    list_last_values[key].append("A")

                elif base == "c" or base == "C":
                    add_to_dict("C", sum(list_correctness[key]), dict_bases)
                    if len(list_last_values[key]) == k:
                        list_last_values[key].pop(0)
                    list_last_values[key].append("C")

                elif base == "t" or base == "T":
                    add_to_dict("T", sum(list_correctness[key]), dict_bases)
                    if len(list_last_values[key]) == k:
                        list_last_values[key].pop(0)
                    list_last_values[key].append("T")

                elif base == "g" or base == "G":
                    add_to_dict("G", sum(list_correctness[key]), dict_bases)
                    if len(list_last_values[key]) == k:
                        list_last_values[key].pop(0)
                    list_last_values[key].append("G")

                elif base == "u" or base == "U":
                    add_to_dict("U", sum(list_correctness[key]), dict_bases)
                    if len(list_last_values[key]) == k:
                        list_last_values[key].pop(0)
                    list_last_values[key].append("U")

                elif base == "n" or base == "N":
                    if len(list_last_values[key]) == k:
                        list_last_values[key].pop(0)
                    list_last_values[key].append("N")

            else:
                if len(list_last_values[key]) == k:
                    list_last_values[key].pop(0)
                list_last_values[key].append("N")

        #print(dict_bases)
        try:
            max_val = max(dict_bases.values())
        except:
            max_val = 0
        max_keys = [k for k, v in dict_bases.items() if v == max_val]

        #print("MAX value", max_val , "Max keys " , max_keys)

        if len(max_keys) == 1:
            consensus.append(max_keys[0])
            update_correctness(max_keys[0], k)
        if max_val == 0:
            consensus.append("N")
            update_correctness("N", k)
        if len(max_keys) > 1:
            if "A" in max_keys:
                consensus.append("A")
                update_correctness("A", k)
            elif "T" in max_keys:
                consensus.append("T")
                update_correctness("T", k)
            elif "C" in max_keys:
                consensus.append("C")
                update_correctness("C", k)
            elif "G" in max_keys:
                consensus.append("G")
                update_correctness("G", k)
            elif "U" in max_keys:
                consensus.append("U")
                update_correctness("U", k)
            else:
                consensus.append("N")
                update_correctness("N", k)
        count += 1


    file = open(output, "w")

    file.write(">CoopPipe_consensus\n" + ''.join(consensus) )

    file.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Index",
    usage="python3 weighted_generate_consensus.py -i <aligned multi-FASTA> -k <values of k>")

    parser.add_argument("-i", help="Aligned multi-FASTA", type=str, required=True)
    parser.add_argument("-k", help="Values of k, separated by spaces", nargs="+", type=int, required=True)
    args = parser.parse_args()


    filename = args.i
    read_file(filename)

    count = 0
    for i in args.k:
        generate_consensus("tmp-" + str(count) + ".fa", i)
        count += 1

    #generate_consensus("tmp-2.fa", 4)
    #generate_consensus("tmp-3.fa", 15)
    #generate_consensus("tmp-4.fa", 30)
    #generate_consensus("tmp-5.fa", 200)
    #generate_consensus("tmp-6.fa", 1000)
    #generate_consensus("tmp-7.fa", 2)
    #generate_consensus("tmp-8.fa", 100)

    os.system('cat tmp-*.fa > new.fa')
    os.system('rm tmp-*.fa')

