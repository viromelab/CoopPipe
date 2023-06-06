#!/usr/bin/env python3
import argparse

dict_content = {}



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
        #print("added new key ", key)

def read_file (path):

    file = open(path, "r")
    count = 0

    for line in file:

        if line[0] != ">":
            #print("line strip  - ", line.strip("\n"))
            add_to_dict(count, line.strip("\n"), dict_content)
            #print(dict_content)
        else:
            count += 1

    #print(dict_content)



    file.close()

def generate_consensus (output):

    count = 0
    finished = 0

    list_bases = "actguACTGU"
    consensus = []
    keys_finished = []

    while finished < len(dict_content):

        dict_bases = {}

        for key in dict_content:

            try:
                base = dict_content.get(key)[count]
            except:
                if key not in keys_finished:
                    keys_finished.append(key)
                    finished += 1


            if base in list_bases: #check if it is one of the bases
                if base == "a" or base == "A":
                    add_to_dict("A", 1, dict_bases)
                elif base == "c" or base == "C":
                    add_to_dict("C", 1, dict_bases)
                elif base == "t" or base == "T":
                    add_to_dict("T", 1, dict_bases)
                elif base == "g" or base == "G":
                    add_to_dict("G", 1, dict_bases)
                elif base == "u" or base == "U":
                    add_to_dict("U", 1, dict_bases)

        #print(dict_bases)
        try:
            max_val = max(dict_bases.values())
        except:
            max_val = 0
        max_keys = [k for k, v in dict_bases.items() if v == max_val]

        #print("MAX value", max_val , "Max keys " , max_keys)

        if len(max_keys) == 1:
            consensus.append(max_keys[0])
        if max_val == 0:
            consensus.append("N")
        if len(max_keys) > 1:
            if "A" in max_keys:
                consensus.append("A")
            elif "T" in max_keys:
                consensus.append("T")
            if "C" in max_keys:
                consensus.append("C")
            if "G" in max_keys:
                consensus.append("G")
        count += 1


    file = open(output, "w")

    #print ("cut_len", len(dict_content.get(1))-1)
    sliced = consensus[: len(dict_content.get(1))-1]
    
    #print (sliced, len(sliced))

    file.write(">CoopPipe_consensus\n" + ''.join(sliced) + "\n" )




    file.close()









if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Index",
    usage="python3 generate_consensus.py -i <aligned multi-FASTA> -o <output file>")

    parser.add_argument("-i", help="Aligned multi-FASTA", type=str, required=True)
    parser.add_argument("-o", help="Output file (FASTA)", type=str)
    args = parser.parse_args()


    filename = args.i
    output = "cooppipe-reconstructed.fa"
    if args.o != None:
        output = args.o
    read_file(filename)
    generate_consensus(output)

