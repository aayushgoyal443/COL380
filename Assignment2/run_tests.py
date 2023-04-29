import subprocess

# for simple commands
# subprocess.run(["ls", "-l"]) 

# for complex commands, with many args, use string + `shell=True`:
# cmd_str = "ls -l /tmp | awk '{print $3,$9}' | grep root"
# subprocess.run(cmd_str, shell=True)

cmd_str = f"make compile_parallel"
subprocess.run(cmd_str, shell=True)

#  ids of small test cases [0, 1, 3, 4, 5]
#  ids of intermediate test cases [6, 10, 11, 12, 14]
#  ids of large test cases [13, 15]


for i in [4]:
    # fetch start and end value from our_ass/test_cases/test{i}/info.txt
    # run the test with the given start and end value
    # read the info.txt file
    file1 = f"test_cases/gradedA2/test{i}/task1_info.txt"

    f1 = open(file1, 'r')
    f1_A = f1.read().strip()
    x = f1_A.split(',')
    start = int(x[0].split()[-1])
    end = int(x[1].split()[-1])
    print(start, end)
    print(f"Compiling test {i}")
    print("Running test " + str(i) + "...")
    cmd_str = f"make test={i} start={start} end={end} run_parallel_parse"
    print("Run command: " + cmd_str)
    subprocess.run(cmd_str, shell=True)
    print("Done")
    print("Checking results...")
    cmd_str = f"python3 compare.py results/result_{i}.txt test_cases/gradedA2/test{i}/task1_output{i}.txt"
    print("Comparing using: " + cmd_str)
    subprocess.run(cmd_str, shell=True)
    print("Done")

