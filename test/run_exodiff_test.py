#! /usr/bin/env python
from __future__ import print_function
import sys
import os
import string
import glob
import argparse as ap
from subprocess import Popen, PIPE, run

def runtestdiff(executable_name, cli_flag, input_deck_name, num_ranks, use_openmpi=False):

    result = 0
    base_name = input_deck_name[:-3]

    command = []

    epu_required = False
    if num_ranks > 1:
        command.append("mpirun")
        command.append("-np")
        command.append(str(num_ranks))
        epu_required = True
        if use_openmpi:
            command.append("--use-hwthread-cpus")
    #
    command.append(executable_name)
    #
    if cli_flag:
        command.append("--" + cli_flag)
    command.append(input_deck_name)

    # differentiate output files based on type of run
    nimble_output_name = "none"
    log_file_name = "none"
    epu_output_extension = ".e"
    if num_ranks > 1:
        epu_output_extension = "np" + str(num_ranks) + ".e"
    epu_exodus_output_name = "none"
    epu_ranks_string = str(num_ranks)
    nimble_output_name = base_name + ".out"
    if not cli_flag:
        log_file_name = base_name + ".mpi.np" + str(num_ranks) + ".log"
    if "use_kokkos" in cli_flag:
        log_file_name = base_name + ".kokkos.np" + str(num_ranks) + ".log"
    if "use_vt" in cli_flag:
        log_file_name = base_name + ".vt.np" + str(num_ranks) + ".log"
    if "use_tpetra" in cli_flag:
        log_file_name = base_name + ".tpetra.np" + str(num_ranks) + ".log"
    if num_ranks > 1:
        epu_exodus_output_name = base_name + ".out." + epu_output_extension
    else:
        epu_exodus_output_name = base_name + ".out.e"

    # open the log file
    if os.path.exists(log_file_name):
        os.remove(log_file_name)
    logfile = open(log_file_name, 'w')

    logfile.write("\nrun_test.py command: " + " ".join(command) + "\n")
    logfile.flush()

    # remove old output files, if any (this will miss some files in parallel runs)
    if os.path.exists(epu_exodus_output_name):
        os.remove(epu_exodus_output_name)

    print("\nCommand:", command)

    # run the code
    p = run(command, stdout=logfile, stderr=logfile, text=True, check=True)
    return_code = p.returncode
    if return_code != 0:
        result = return_code
    logfile.write("------------FIRST TEST RESULT " + str(result) + "\n\n")

    # run epu
    if epu_required:
        command = ["epu", \
                   "-p", \
                   epu_ranks_string, \
                   "-output_extension", \
                   epu_output_extension, \
                   nimble_output_name]
        print("EPU COMMAND", command)
        p = run(command, stdout=logfile, stderr=logfile, text=True, check=True)
        return_code = p.returncode
        if return_code != 0:
            result = return_code
    logfile.write("------------SECOND(MPI) TEST RESULT " + str(result) + "\n\n")
    # run exodiff
    command = ["exodiff", \
               "-stat", \
               "-f", \
               base_name+".exodiff", \
               base_name+".gold.e", \
               epu_exodus_output_name]
    p = run(command, stdout=logfile, stderr=logfile, text=True, check=True)
    return_code = p.returncode
    if return_code != 0:
        result = return_code

    logfile.write("FINAL TEST RESULT " + str(result) + "\n\n")

    logfile.close()

    return result

if __name__ == "__main__":

    parser = ap.ArgumentParser(description='run_test.py', prefix_chars='-')
    parser.add_argument('--executable', required=True, action='store', nargs=1, metavar='executable', help='Name of NimbleSM executable')
    parser.add_argument('--cli-flag', required=True, action='store', nargs=1, metavar='cli_flag', help='Command Line flags')
    parser.add_argument('--input-deck', required=True, action='store', nargs=1, metavar='input_deck', help='NimbleSM input deck (*.in)')
    parser.add_argument('--num-ranks', required=False, type=int, action='store', nargs=1, metavar='num_ranks', help='Number of physical ranks')

    args = vars(parser.parse_args())
    print(args)

    executable = args['executable'][0]
    cli_flag = args['cli_flag'][0]
    input_deck = args['input_deck'][0]

    num_ranks = 1
    if args['num_ranks'] != None:
        num_ranks = args['num_ranks'][0]

    result = runtestdiff(executable,
                     cli_flag,
                     input_deck,
                     num_ranks)

    sys.exit(result)
