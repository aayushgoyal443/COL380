import subprocess
import sys

task = int(sys.argv[1])
verbose = int(sys.argv[2])

print(f"Running tests for task {task} with verbose={verbose}")

# for simple commands
# subprocess.run(["ls", "-l"]) 

# for complex commands, with many args, use string + `shell=True`:
# cmd_str = "ls -l /tmp | awk '{print $3,$9}' | grep root"
# subprocess.run(cmd_str, shell=True)

cmd_str = f"make"
subprocess.run(cmd_str, shell=True)

for i in [1,2,3,4]:
    # fetch start and end value from our_ass/test_cases/test{i}/info.txt
    # run the test with the given start and end value
    # read the info.txt file
    file1 = f"test_cases/A3/test{i}/task{task}_info.txt"

    f1 = open(file1, 'r')
    f1_A = f1.read().strip()
    x = f1_A.split(',')
    start = int(x[0].split()[-1])
    end = int(x[1].split()[-1])
    print(start, end)
    print(f"Compiling test {i}")
    print("Running test " + str(i) + "...")
    cmd_str = f"make test={i} start={start} end={end if task==1 else start} p={end} run_task{task}{'_verbose' if verbose==1 else ''}"
    print("Run command: " + cmd_str)
    subprocess.run(cmd_str, shell=True)
    print("Done")
    print("Checking results...")
    cmd_str = f"python3 compare.py results/task{task}/result_{i}{'_verbose' if verbose==1 else ''}.txt test_cases/A3/test{i}/task{task}_output{i}{'_verbose' if verbose==1 else ''}.txt"
    print("Comparing using: " + cmd_str)
    subprocess.run(cmd_str, shell=True)
    print("Done")


cmd_str = f"make clean"
subprocess.run(cmd_str, shell=True)