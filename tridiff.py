__author__ = 'gbb'

import os
import datetime
import subprocess
import tempfile
import time
import argparse

def get_modification_dates(filenames):
    moddates = []
    for filename in filenames:
        if os.path.isfile(filename):
            t = os.path.getmtime(filename)
            moddates.append(datetime.datetime.fromtimestamp(t))
        else:
            moddates.append(None)
    return moddates


def decorate_string(s):
    return "\n-------" + s + "-------\n"


def execute_command(cmds, i):
    print decorate_string("running cmd" + str(i))
    os.system(cmds[i])


def backup_files(filenames):
    backups = []
    for filename in filenames:
        backuped_name = filename + ".tridiffbak"
        if os.path.isfile(filename):
            backups.append(backuped_name)
            os.system("cp " + filename + " " + backuped_name)
    return backups


def delete_files(filenames):
    for filename in filenames:
        os.system("rm " + filename)


def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-0', dest='cmd0', nargs='*')
    parser.add_argument('-1', dest='cmd1', nargs='*')
    parser.add_argument('-f', dest='files', nargs='*')

    args = parser.parse_args()
    cmds = [' '.join(args.cmd0), ' '.join(args.cmd1)]
    files = []
    for file in args.files:
        files.append(file)

    path = os.path.dirname(os.path.realpath(__file__)) + '/'
    # path = '/Users/gbb/Desktop/code/GraduateThesis/tri/'

    cmds = [w.replace('./', path) for w in cmds]
    files = [w.replace('./', path) for w in files]

    return cmds, files


#  here it starts
cmds, filenames = parse_options()

#  phase 0
moddates0 = get_modification_dates(filenames)  # get modification dates before program runs

#  phase 1
execute_command(cmds, 0)
moddates1 = get_modification_dates(filenames)
backup_file_names = backup_files(filenames)
time.sleep(1)  # in case the programs run so fast that the change in modification time does not show

#  phase 2
execute_command(cmds, 1)
moddates2 = get_modification_dates(filenames)

#  output results
print decorate_string("results")
for i in range(0, len(filenames)):
    if moddates0[i] == moddates1[i]:
        print "cmd0 did not modify " + filenames[i]
        continue
    if moddates1[i] == moddates2[i]:
        print "cmd1 did not modify " + filenames[i]
        continue
    with tempfile.TemporaryFile() as tempf:
        proc = subprocess.Popen(["diff", backup_file_names[i], filenames[i]], stdout=tempf)
        proc.wait()
        tempf.seek(0)
        if len(tempf.read()) == 0:
            print "old and new " + filenames[i] + " are identical"
        else:
            print "change in " + filenames[i]
            print tempf.read()

delete_files(backup_file_names)