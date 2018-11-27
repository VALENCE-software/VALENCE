#!/usr/bin/env python
"""
Contains IO and OS related tools.
Since python have different modules for various
IO related functionalities, it is good to have
a single module simplifying their usage.
TODO: Add unit tests. Seperate IO vs OS
"""
import time
import os
from os.path import isfile
import logging
__updated__ = "2018-03-13"
__author__  = "Murat Keceli"

def get_date():
    """
    Returns the current date and time
    """
    import datetime
    return datetime.datetime.now()

def touch(fname, times=None):
    """
    Creates a file with the given fname, aka unix touch
    See http://stackoverflow.com/questions/1158076/implement-touch-using-python
    """
    import os
    try:
        with open(fname, 'a'):
            os.utime(fname, times)
    except:
        pass
    return


def rm(fname):
    """
    Deletes a file with the given fname.
    """
    import os
    if check_file(fname):
        os.remove(fname)
    else:
        logging.debug('Cannot delete file. {} does not exist.'.format(fname))
    return


def rmrf(dname):
    """
    Deletes directories recursively including all of the contents.
    """
    import shutil
    try:
        shutil.rmtree(dname,ignore_errors=True)
    except OSError as e:
        logging.error("Error in deleting {}: {}" % (dname,e.strerror))
    return


def mkdir(path):
    """
    Creates directory if it doesn't exist.
    Check http://stackoverflow.com/questions/273192/how-to-check-if-a-directory-exists-and-create-it-if-necessary
    """
    import os
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    return


def cd(path):
    """
    Change working directory
    """
    import os
    os.chdir(path)
    return


def mv(oldname,newname):
    """
    Renames or moves a file
    """
    import os
    try:
        os.rename(oldname,newname)
    except:
        logging.debug('Command "mv {0} {1}" failed' .format(oldname, newname))
    return

def pwd():
    """
    Return current directory
    """
    import os
    return os.getcwd()


def find_files_recursive(directory, pattern):
    """
    Yields files in a directory (including subdirectories) with a given pattern
    https://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
    """
    import os, fnmatch
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

def find_files(directory,pattern):
    """
    Returns a list of files that matches the pattern in a given directory
    """
    import glob
    files = glob.glob(directory + '/' + pattern)
    return files

def join_path(*paths):
    """
    Concatenes strings into a portable path using correct seperators.
    """
    import os
    return os.path.join(*paths)


def write_file(s, filename='newfile'):
    """
    Writes s string to a file with the given'filename'.
    """
    with open(filename, 'w') as f:
        f.write(s)
    return

def append_file(s, filename='newfile'):

    """
    Appends s string to a file with the given 'filename'.
    """
    with open(filename, 'a') as f:
        f.write(s)
    return


def read(inp):
    """
    Returns a string that contains the text for the given input.
    inp can be the filename.
    """
    lines = inp.splitlines()
    if len(lines) == 1:
        if check_file(inp):
            tmp = read_file(inp)
        else:
            tmp = inp
    else:
        tmp = inp
    return tmp


def read_file(filename, aslines=False):
    """
    Reads a file and return either a list of lines or a string.
    """
    with open(filename, 'r') as f:
        if aslines:
            tmp = f.readlines()
        else:
            tmp = f.read()
    return tmp


def get_unique_filename(fname):
    """
    Try to return a unique filename.
    Not safe for race conditions.
    """
    counter = 1
    uname = fname
    while check_file(uname,timeout=0.01):
        uname = fname + '_' + str(counter)
        counter += 1
    return uname


def fix_path(s):
    """
    Returns a path with problematic characters replaced by safer ones.
    """
    s = s.replace('[','_b_')
    s = s.replace(']','_d_')
    #s = s.replace('=','_e_')
    s = s.replace(':','_i_')
    s = s.replace('|','_j_')
    s = s.replace('\\','_k_')
    #s = s.replace('/','_l_')
    s = s.replace('(','_p_')
    s = s.replace(')','_q_')
    s = s.replace('*','_s_')
    #s = s.replace('#','_x_')
    s = s.replace('<','_v_')
    s = s.replace('>','_y_')
    s = s.replace('?','_z_')
    return s


def check_file(filename, timeout=0, verbose=False):
    """
    Returns True (False) if a file exists (doesn't exist).
    If timeout>0 is given, then checks file in a loop until
    timeout seconds pass.
    If verbose is True and file does not exist, prints an error message
    """
    import time
    from os.path import isfile
    exists = isfile(filename)
    if timeout > 0 and not exists:
        t0 = time.time()
        while time.time() < t0 + timeout:
            if isfile(filename):
                exists = True
                break
    if not exists and verbose:
        print('"{0}" file not found.'.format(filename))
        logging.debug('"{0}" file not found.'.format(filename))
    return exists


def check_dir(dirname, timeout=0):
    """
    Returns True (False) if a file exists (doesn't exist).
    If timeout>0 is given, then checks file in a loop until
    timeout seconds pass.
    """
    import time
    from os.path import isdir
    exists = isdir(dirname)
    if timeout > 0 and not exists:
        t0 = time.time()
        while time.time() < t0 + timeout:
            if isdir(dirname):
                exists = True
                break
    return exists


def check_exe(exename):
    """
    Check if an executable is available.
    TODO:
    """
    import distutils.spawn as ds

    return ds.find_executable(exename)


def cp(source, target):
    """
    Copies a file, source and target are strings for paths.
    Target can be a directory or a file.
    """
    from shutil import copy
    copy(source, target)
    return


def symlink(source,linkname):
    """
    Create a symbolic link. (works only on unix based systems)
    """
    if check_file(linkname):
        logging.debug('Target link {} exists, not creating a new one'.format(linkname))
    else:
        try:
            os.symlink(source,linkname)
        except Exception as e:
            logging.warning('Symlink exception caught: {}'.format(e))
    return


def get_path(f, executable=False, directory=False):
    """
    Returns absolute path for a file or folder.
    """
    import os
    import distutils.spawn as ds
    if executable:
        return ds.find_executable(f)
    elif directory:
        return os.path.dirname(os.path.abspath(f))
    else:
        return os.path.abspath(f)


def get_line_number(keyword, lines=None, filename=None,getlastone=False):
    """
    Returns the line number of a keyword found in given lines of string.
    Returns -1 if keyword is not found
    """
    num = -1
    if lines is None and filename is None:
        logging.debug('List of lines or a filename to be read is required for get_line_number')
    elif filename:
        lines = read_file(filename, aslines=True)

    for n, line in enumerate(lines):
        if keyword in line:
            if getlastone:
                num = n
            else:
                return n
    return num


def get_line_numbers(keyword, lines=None, filename=None):
    """
    Returns the line numbers as a list for a keyword in given lines of string.
    Returns -1 if keyword is not found
    """
    if lines is None and filename is None:
        logging.debug('List of lines or a filename to be read is required for get_line_numbers')
    elif filename:
        lines = read_file(filename, aslines=True)
    nums = []
    for n, line in enumerate(lines):
        if keyword in line:
            nums.append(n)
    if len(nums) == 0:
        nums = -1
    return nums


def get_git_version():
    """
    Return git version.
    """
    import subprocess
    return subprocess.check_output(["git", "describe"])


def read_list(listfile):
    """
    Return a list of strings from all lines in a text file.
    Skips blank lines.
    """
    with open(listfile, 'r') as f:
        lines = filter(None, ((line.rstrip()) for line in f))
    return lines


def read_list2(listfile):
    """
    Return a list of strings from all lines in a text file.
    Skips blank lines and lines that start with '#' or '!'.
    Ignores the portion of lines after a space.
    """
    lines = []
    with open(listfile, 'r') as f:
        for line in f:
            aline = line.strip()
            if aline:
                if aline.startswith('#') or aline.startswith('!'):
                    pass
                else:
                    theline = aline.split()[0].strip()
                    lines.append(theline)
    return lines


def execute_old(exe,inp=None,out=None):
    """
    Executes a calculation for a given input file inp and executable exe.
    """
    from subprocess import Popen, PIPE
    if not check_exe(exe):
        return 'Executable "{0}" not found.\n'.format(exe)
    if inp:
        if not check_file(inp,1):
            return 'Input file "{0}" not found.\n'.format(get_path(inp))
    if inp and out:
        process = Popen([exe, inp, out], stdout=PIPE, stderr=PIPE)
    elif inp:
        process = Popen([exe, inp], stdout=PIPE, stderr=PIPE)
    else:
        process = Popen([exe], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    if stderr is None or stderr == '':
        msg = 'Run {0} {1}: Success.\n'.format(exe, inp)
    else:
        errstr = """ERROR in "{0}"\n
        STDOUT:\n{1}\n
        STDERR:\n{2}""".format(inp, stdout, stderr)
        errfile = inp + '.err'
        write_file(errstr, errfile)
        msg = 'Run {0} {1}: Failed, see "{2}"\n'.format(exe, inp, get_path(errfile))
    return msg


def get_mpi_rank(default=None):
    """
    Return mpi rank (int) if defined as an environment variable
    """
    if os.getenv("PMI_RANK") is not None:
        rank = int(os.getenv("PMI_RANK"))
    elif os.getenv("OMPI_COMM_WORLD_RANK") is not None:
        rank = int(os.getenv("OMPI_COMM_WORLD_RANK"))
    else:
        rank = default
    return rank


def get_mpi_size(default=None):
    """
    Return mpi size (int) if defined as an environment variable
    """
    if os.getenv("PMI_SIZE") is not None:
        size = int(os.getenv("PMI_SIZE"))
    elif os.getenv("OMPI_COMM_WORLD_SIZE") is not None:
        size = int(os.getenv("OMPI_COMM_WORLD_SIZE"))
    else:
        size = int(default)
    return size


def get_mpi_local_rank(default=None):
    """
    Return mpi local rank as an integer if defined as an environment variable
    https://www.open-mpi.org/faq/?category=running#mpi-environmental-variables
    The relative rank of this process on this node within its job.
    For example, if four processes in a job share a node, they will each be given a local rank ranging from 0 to 3.
    """
    if os.getenv("OMPI_COMM_WORLD_LOCAL_RANK") is not None:
        rank = int(os.getenv("OMPI_COMM_WORLD_LOCAL_RANK"))
    else:
        rank = default
    return rank


def get_mpi_local_size(default=None):
    """
    Return mpi local size as an integer if defined as an environment variable
    https://www.open-mpi.org/faq/?category=running#mpi-environmental-variables
    The number of processes on this node within its job.
    """
    if os.getenv("OMPI_COMM_WORLD_LOCAL_SIZE") is not None:
        size = int(os.getenv("OMPI_COMM_WORLD_LOCAL_SIZE"))
    else:
        size = default
    return size


def get_ppn(includethreads=False):
    """
    Return number of processors per node.
    For alternative solutions:
    https://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-using-python
    """
    import psutil
    return psutil.cpu_count(includethreads)


def get_total_memory():
    """
    Return totol physical memory in MB as an integer.
    """
    from psutil import virtual_memory

    mem = virtual_memory() # In bytes
    m = mem.total >> 20 # Using bit shift to get in MB
    #m = mem.total >> 30 # Using bit shift to get in GB
    return m


def get_env(var,default=None):
    """
    Return the value of environment variable.
    """
    val = os.getenv(var)
    if val is None:
        val = default
    return val

def run(inp, exe):
    """
    """
    from subprocess import Popen, PIPE
    process = Popen(exe, stdin= PIPE, stdout=PIPE, stderr=PIPE)
    out, err = process.communicate(input=inp)
    if err:
        logging.error('ERROR in iotools.run: {}'.format(err))
    return out


def execute(command, stdoutfile=None, stderrfile=None, merge=False):
    """
    Executes a given command, and optionally write stderr and/or stdout.
    Parameters
    ----------
    command: List of strings, where a command line is seperated into words.
    stderrfile: None or a string for a file name to write stderr
    stdoutfile: None or a string for a file name to write stdout

    Returns
    ---------
    If stdoutfile:
        A string describing the success or failure of the calculation
    Else:
        stdout

    Doctest
    ---------
    >>> io.execute(['echo','this works'])
    'this works\n'
    >>> io.execute('echo this also works')
    'this also works\n'
    """
    from subprocess import Popen, PIPE
    if type(command) == str:
        commandstr = command
        command = command.split()
    else:
        commandstr = ' '.join(command)
    msg = 'Running Popen with command: {0}\n'.format(commandstr)
    logging.debug(msg)
    msg =''
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    if merge:
        if type(out) == str and type(err) == str:
            out += err
    if out is None or out == '':
        pass
    else:
        if stdoutfile:
            write_file(out,stdoutfile)
            msg += 'STDOUT is appended to {0}\n'.format(stdoutfile)
        else:
            msg += 'STDOUT:\n{0}\n'.format(out)
    if err is None or err == '':
        pass
    else:
        if stderrfile:
            write_file(err, stderrfile)
            msg += 'STDERR file:\n"{0}"\n'.format(get_path(err))
        else:
            msg += 'STDERR:\n{0}\n'.format(err)
    return msg


def get_stdout_stderr(command):
    """
    Executes a given command, and get stdout and stderr.
    Parameters
    ----------
    command: A string or a list of strings, where a command line is seperated into words.

    Returns
    ---------
    stdout, stderr

    Doctest
    ---------
    >>> get_stdout_stderr(['echo','this works'])
    'this works\n'
    >>> get_stdout_stderr('echo this also works')
    'this also works\n'
    """
    from subprocess import Popen, PIPE
    if type(command) == str:
        commandstr = command
        command = command.split()
    else:
        commandstr = ' '.join(command)
    logging.debug('Running Popen with command: {0}'.format(commandstr))
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    return out, err

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
