#! /usr/bin/env python

import sys
import os
import string
import glob
import argparse as ap
from subprocess import Popen

def runtest(executable_name, input_deck_name, num_ranks, num_virtual_ranks, have_charm, qthreads_num_shepherds, qthreads_num_workers_per_shepherd):

    result = 0
    base_name = input_deck_name[:-3]

    command = []

    if charm_run != None:
        command.append(charm_run)
        command.append("+p")
        command.append(str(num_ranks))
        command.append("++local")

    epu_required = False
    if "NimbleSM_Serial" in executable_name:
        command.append(executable_name)
    if "NimbleSM_MPI" in executable_name or "NimbleSM_Kokkos" in executable_name or "NimbleSM_Tpetra" in executable_name or "NimbleSM_Qthreads" in executable_name:
        command.append("mpirun")
        command.append("-np")
        if num_ranks:
            command.append(str(num_ranks))
            epu_required = True
        else:
            command.append("1")
        command.append(executable_name)
    if "NimbleSM_Qthreads" in executable_name:
        if qthreads_num_shepherds:
            command.append("-num_shepherds")
            command.append(str(qthreads_num_shepherds))
        if qthreads_num_workers_per_shepherd:
            command.append("-num_workers_per_shepherd")
            command.append(str(qthreads_num_workers_per_shepherd))
        if int(qthreads_num_shepherds) > 0 or int(qthreads_num_workers_per_shepherd) > 0:
            epu_required = True
    if "NimbleSM_Darma" in executable_name:
        command.append(executable_name)
        if num_ranks:
            command.append("--ranks")
            command.append(str(num_ranks))
        command.append("--app-argv")
        if num_virtual_ranks:
            command.append("num-virtual-ranks")
            command.append(str(num_virtual_ranks))
            epu_required = True
    command.append(input_deck_name)

    # differentiate output files based on type of run
    nimble_output_name = "none"
    log_file_name = "none"
    epu_output_extension = "none"
    epu_exodus_output_name = "none"
    epu_ranks_string = "none"
    if "NimbleSM_Serial" in executable_name:
        nimble_output_name = base_name
        log_file_name = base_name + ".serial.log"
        epu_exodus_output_name = base_name + ".serial.e"
    if "NimbleSM_MPI" in executable_name:
        nimble_output_name = base_name + ".mpi"
        log_file_name = base_name + ".mpi.np" + str(num_ranks) + ".log"
        if num_ranks > 1:
            epu_output_extension = "np" + str(num_ranks) + ".e"
            epu_exodus_output_name = base_name + ".mpi." + epu_output_extension
        else:
            epu_output_extension = ".e"
            epu_exodus_output_name = base_name + ".mpi.e"
        epu_ranks_string = str(num_ranks)
    if "NimbleSM_Kokkos" in executable_name:
        nimble_output_name = base_name + ".kokkos"
        log_file_name = base_name + ".kokkos.np" + str(num_ranks) + ".log"
        if num_ranks > 1:
            epu_output_extension = "np" + str(num_ranks) + ".e"
            epu_exodus_output_name = base_name + ".kokkos." + epu_output_extension
        else:
            epu_output_extension = ".e"
            epu_exodus_output_name = base_name + ".kokkos.e"
        epu_ranks_string = str(num_ranks)
    if "NimbleSM_Qthreads" in executable_name:
        nimble_output_name = base_name + ".qthreads"
        log_file_name = base_name + ".qthreads.np" + str(num_ranks) + ".ns" + str(qthreads_num_shepherds) + ".nwps" + str(qthreads_num_workers_per_shepherd) + ".log"
        if num_ranks > 1 or qthreads_num_shepherds > 1 or qthreads_num_workers_per_shepherd > 1:
            epu_output_extension = "np" + str(num_ranks) + ".ns" + str(qthreads_num_shepherds) + ".nwps" + str(qthreads_num_workers_per_shepherd) + ".e"
            epu_exodus_output_name = base_name + ".qthreads." + epu_output_extension
        else:
            epu_output_extension = ".e"
            epu_exodus_output_name = base_name + ".qthreads.e"
        epu_ranks_string = str(num_ranks*qthreads_num_shepherds*qthreads_num_workers_per_shepherd)
    if "NimbleSM_Tpetra" in executable_name:
        nimble_output_name = base_name + ".tpetra"
        log_file_name = base_name + ".tpetra.np" + str(num_ranks) + ".log"
        if num_ranks > 1:
            epu_output_extension = "np" + str(num_ranks) + ".e"
            epu_exodus_output_name = base_name + ".tpetra." + epu_output_extension
        else:
            epu_output_extension = ".e"
            epu_exodus_output_name = base_name + ".tpetra.e"
        epu_ranks_string = str(num_ranks)
    if "NimbleSM_Darma" in executable_name:
        nimble_output_name = base_name + ".darma"
        log_file_name = base_name + ".darma.np" + str(num_virtual_ranks) + ".log"
        if num_ranks > 1:
            epu_output_extension = "np" + str(num_virtual_ranks) + ".e"
            epu_exodus_output_name = base_name + ".darma." + epu_output_extension
        else:
            epu_output_extension = ".e"
            epu_exodus_output_name = base_name + ".darma.e"
        epu_ranks_string = str(num_virtual_ranks)

    # open the log file
    if os.path.exists(log_file_name):
        os.remove(log_file_name)
    logfile = open(log_file_name, 'w')

    logfile.write("\nrun_exodiff_test.py command: " + " ".join(command) + "\n")
    logfile.flush()

    # remove old output files, if any (this will miss some files in parallel runs)
    if os.path.exists(epu_exodus_output_name):
        os.remove(epu_exodus_output_name)

    print("\nCommand:", command)

    # run the code
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # run epu
    if epu_required:
        command = ["epu", \
                   "-p", \
                   epu_ranks_string, \
                   "-output_extension", \
                   epu_output_extension, \
                   nimble_output_name]
        print("EPU COMMAND", command)
        p = Popen(command, stdout=logfile, stderr=logfile)
        return_code = p.wait()
        if return_code != 0:
            result = return_code

    # run exodiff
    command = ["exodiff", \
               "-stat", \
               "-f", \
               base_name+".exodiff", \
               base_name+".gold.e", \
               epu_exodus_output_name]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    logfile.write("FINAL TEST RESULT " + str(result) + "\n\n")

    logfile.close()

    return result

if __name__ == "__main__":

    parser = ap.ArgumentParser(description='run_exodiff_test.py', prefix_chars='-')
    parser.add_argument('--executable', required=True, action='store', nargs=1, metavar='executable', help='Name of NimbleSM executable')
    parser.add_argument('--input-deck', required=True, action='store', nargs=1, metavar='input_deck', help='NimbleSM input deck (*.in)')
    parser.add_argument('--num-ranks', required=False, type=int, action='store', nargs=1, metavar='num_ranks', help='Number of physical ranks')
    parser.add_argument('--num-virtual-ranks', required=False, type=int, action='store', nargs=1, metavar='num_virtual_ranks', help='Number of virtual ranks')
    parser.add_argument('--charm-run', required=False, action='store', nargs=1, metavar='charm_run', help='Charm++ executable (charmrun)')
    parser.add_argument('--qthreads-num-shepherds', required=False, type=int, action='store', nargs=1, metavar='qthread-num-shepherds', help='Number of qthreads shepherds')
    parser.add_argument('--qthreads-num-workers-per-shepherd', required=False, type=int, action='store', nargs=1, metavar='qthread-num-workers-per-shepherd', help='Number of qthreads workers per shepherd')

    args = vars(parser.parse_args())
    print(args)

    executable = args['executable'][0]

    input_deck = args['input_deck'][0]

    num_ranks = 1
    if args['num_ranks'] != None:
        num_ranks = args['num_ranks'][0]

    num_virtual_ranks = 1
    if args['num_virtual_ranks'] != None:
        num_virtual_ranks = args['num_virtual_ranks'][0]

    charm_run = None
    if args['charm_run'] != None:
        if args['charm_run'][0] != "not available":
            charm_run = args['charm_run'][0]

    qthreads_num_shepherds = None
    if args['qthreads_num_shepherds'] != None:
        qthreads_num_shepherds = args['qthreads_num_shepherds'][0]

    qthreads_num_workers_per_shepherd = None
    if args['qthreads_num_workers_per_shepherd'] != None:
        qthreads_num_workers_per_shepherd = args['qthreads_num_workers_per_shepherd'][0]

    result = runtest(executable,
                     input_deck,
                     num_ranks,
                     num_virtual_ranks,
                     charm_run,
                     qthreads_num_shepherds,
                     qthreads_num_workers_per_shepherd)

    sys.exit(result)
