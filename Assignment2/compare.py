# a python script which compares two files and prints the differences

import sys
def prRed(skk): print("\033[91m {}\033[00m" .format(skk))
 
 
def prGreen(skk): print("\033[92m {}\033[00m" .format(skk))
 
 
def prYellow(skk): print("\033[93m {}\033[00m" .format(skk))
 
 
def prLightPurple(skk): print("\033[94m {}\033[00m" .format(skk))
 
 
def prPurple(skk): print("\033[95m {}\033[00m" .format(skk))
 
 
def prCyan(skk): print("\033[96m {}\033[00m" .format(skk))
 
 
def prLightGray(skk): print("\033[97m {}\033[00m" .format(skk))
 
 
def prBlack(skk): print("\033[98m {}\033[00m" .format(skk))
 

def special_demand(ashish, aayra, file1, file2):
    ashish_sort = ashish
    aayra_sort = aayra
    # these are lists of lists, sort them according to last element as the key
    ashish_sort.sort(key=lambda x: int(x[0]))
    aayra_sort.sort(key=lambda x: int(x[0]))
    # now write the output to a file which has the same name as the input file with _sort at the end
    f3 = open(file1[:-4] + "_sort.txt", 'w')
    f4 = open(file2[:-4] + "_sort.txt", 'w')
    for i in range(len(ashish_sort)):
        f3.write(str(ashish_sort[i]))
        f3.write("\n")
    for i in range(len(aayra_sort)):
        f4.write(str(aayra_sort[i]))
        f4.write("\n")




def compare(file1, file2):

    ashish = []
    aayra = []
    f1 = open(file1, 'r')
    f2 = open(file2, 'r')
    ashu1 = f1.read()
    aayra1 = f2.read()
    ashish = ashu1.splitlines()
    aayra = aayra1.splitlines()

    ashish = ashish[1:]
    # aayra = aayra[1:]
    
    for i in range(len(ashish)):
        ashish[i] = ashish[i].strip().split()
    for i in range(len(aayra)):
        aayra[i] = aayra[i].strip().split()

    # sort the lists
    ashish.sort()
    aayra.sort()

    special_demand(ashish, aayra, file1, file2)
    
    # # open 2 files and write these lists to them
    # f1 = open("ashish_new.txt", 'w')
    # f2 = open("aayra_new.txt", 'w')
    # for i in range(len(ashish)):
    #     f1.write(str(ashish[i]))
    #     f1.write("\n")
    # for i in range(len(aayra)):
    #     f2.write(str(aayra[i]))
    #     f2.write("\n")
    prYellow(f"{file1} has length: " +  str(len(ashish)))
    prYellow(f"{file2} has length: " + str(len(aayra)))

    if (len(ashish) != len(aayra)):
        prRed("Files are NOT same LENGTH")
        sys.exit(1)

    for i in range(len(ashish)):
        if (ashish[i] != aayra[i]):
            prRed("Files are NOT same")
            prRed("Line number: "+  str(i))
            prCyan("Ashish's file: " + str(ashish[i]))
            prCyan("Aayra's file: "+  str(aayra[i]))
            sys.exit(1)

    prGreen("Files are same!")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: compare.py file1 file2")
        sys.exit(1)
    compare(sys.argv[1], sys.argv[2])
