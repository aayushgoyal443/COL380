import subprocess
import time

# for simple commands
# subprocess.run(["ls", "-l"]) 

# for complex commands, with many args, use string + `shell=True`:
# cmd_str = "ls -l /tmp | awk '{print $3,$9}' | grep root"
# subprocess.run(cmd_str, shell=True)

cmd_str = f"make compile_parallel"
subprocess.run(cmd_str, shell=True)

# for i in [0,1,2,4,5,6,7,8,3,9,10,11,12]:
for i in [5,7,3]:
    # fetch start and end value from our_ass/test_cases/test{i}/info.txt
    # run the test with the given start and end value
    # read the info.txt file
    file1 = f"test_cases/A2/test{i}/info.txt"

    f1 = open(file1, 'r')
    f1_A = f1.read()
    x = f1_A.splitlines()
    start = int(x[0].split()[-1])
    end = int(x[1].split()[-1])
    print(start, end)
    print(f"Compiling test {i}")
    print("Running test " + str(i) + "...")
    cmd_str = f"make test={i} start={start} end={end} run_parallel_parse_verbose"
    print("Run command: " + cmd_str)
    start_time = time.time()
    subprocess.run(cmd_str, shell=True)
    end_time = time.time()
    print("Time taken: " + str(end_time - start_time) + " seconds")
    print("Done")
    print("Checking results...")
    cmd_str = f"python3 compare.py results/result_{i}_verbose.txt test_cases/A2/test{i}/output{i}_verbose.txt"
    print("Comparing using: " + cmd_str)
    subprocess.run(cmd_str, shell=True)
    print("Done")

