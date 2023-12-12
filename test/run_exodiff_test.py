#! /usr/bin/env python3
from __future__ import print_function
import sys
import os
import string
import glob
import argparse as ap
import subprocess
import threading
import queue
from typing import (
    List,
    TextIO,
    Tuple,
    Optional
)

def _enqueue_piped_output(q: queue.Queue, pipe: TextIO):
    counts = 0
    fileno = pipe.fileno()
    try:
        # Keep polling for input on the pipe until we are done
        with pipe:
            for line in pipe:
                counts += 1
                q.put((pipe, line))
    finally:
        # Mark the end of the stream
        q.put((pipe, None))


def _echo_enqueued_output(q: queue.Queue, logfiles: List[TextIO], proc: subprocess.Popen) -> None:
    stdout_done = False
    stderr_done = False
    stdout_count = 0
    stderr_count = 0
    stdout_fileno = proc.stdout.fileno()
    stderr_fileno = proc.stderr.fileno()
    while (not stdout_done) or (not stderr_done):
        entry: Tuple[TextIO, Optional[str]] = q.get()
        pipe, contents = entry

        # if contents is None, that stream is done
        if contents is None:
            if pipe is proc.stdout:
                stdout_done = True
            elif pipe is proc.stderr:
                stderr_done = True
            else:
                raise ValueError("invalid pipe")
            q.task_done()
        else:
            # Echo to stdout/stderr
            if pipe is proc.stdout:
                stdout_count += 1
                sys.stdout.write(contents)
            elif pipe is proc.stderr:
                stderr_count += 1
                sys.stderr.write(contents)
            else:
                raise ValueError("invalid pipe")

            # Echo to all our logfiles
            for log in logfiles:
                log.write(contents)
            q.task_done()


def run_cmd(cmd: List[str], logfiles: List[TextIO]) -> None:
    # open the process with text and bufsize=1 (line bufferd)
    proc = subprocess.Popen(cmd, text=True, bufsize=1,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    io_queue = queue.Queue()

    out_thread = threading.Thread(target=_echo_enqueued_output, args=[io_queue, logfiles, proc])
    stdout_thread = threading.Thread(target=_enqueue_piped_output, args=[io_queue, proc.stdout])
    stderr_thread = threading.Thread(target=_enqueue_piped_output, args=[io_queue, proc.stderr])

    out_thread.start()
    stdout_thread.start()
    stderr_thread.start()

    # We will keep enqueuing work at this point so block until queue is drained
    io_queue.join()

    # Join all our outstanding threads
    out_thread.join()
    stdout_thread.join()
    stderr_thread.join()

    # Now, just in case our subprocess is still running, wait for it
    if proc.wait() != 0:
        raise subprocess.CalledProcessError(proc.returncode, proc.args)


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
    command.append(input_deck_name)
    if cli_flag:
        command.append("--" + cli_flag)
        if "use_vt" in cli_flag:
            command.append("--vt_no_terminate")

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
    try:
        run_cmd(command, [logfile])
    except subprocess.CalledProcessError as e:
        print(f"NimbleSM run failed with args {e.args} and with exit code {e.returncode}", file=sys.stdout)
        return e.returncode

    # run epu
    if epu_required:
        command = ["epu", \
                   "-p", \
                   epu_ranks_string, \
                   "-output_extension", \
                   epu_output_extension, \
                   nimble_output_name]
        print("EPU COMMAND", command)
        try:
            run_cmd(command, [logfile])
        except subprocess.CalledProcessError as e:
            print(f"epu run failed with args {e.args} and with exit code {e.returncode}", file=sys.stdout)
            return e.returncode

    # run exodiff
    command = ["exodiff", \
               "-stat", \
               "-f", \
               base_name+".exodiff", \
               base_name+".gold.e", \
               epu_exodus_output_name]
    try:
        run_cmd(command, [logfile])
    except subprocess.CalledProcessError as e:
        print(f"exodiff run failed with args {e.args} and with exit code {e.returncode}", file=sys.stdout)
        return e.returncode

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
