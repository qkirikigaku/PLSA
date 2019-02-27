import sys
import os
from multiprocessing import Pool
import subprocess

def main():
    arguments = []
    for i in range(1,20):
        arguments.append((str(i), mut_file, experiment))
    pool = Pool()
    _ = pool.starmap(execute, arguments)

def execute(num_topic, mut_file, experiment):
    result_path = 'result/' + mut_file + '_' + experiment
    if(os.path.exists(result_path) == False): os.mkdir(result_path)
    cmd = 'bin/MS ' + num_topic + ' ' + mut_file + ' ' + experiment
    subprocess.call(cmd.split())

if __name__ == '__main__':
    args = sys.argv
    mut_file = args[1]
    experiment = args[2]
    main()
