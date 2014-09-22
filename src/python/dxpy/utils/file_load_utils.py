# Copyright (C) 2014 DNAnexus, Inc.
#
# This file is part of dx-toolkit (DNAnexus platform client libraries).
#
#   Licensed under the Apache License, Version 2.0 (the "License"); you may not
#   use this file except in compliance with the License. You may obtain a copy
#   of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
#   WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
#   License for the specific language governing permissions and limitations
#   under the License.

'''
This module provides support for file download and upload. It calculates the
   location of the input and output directories. It also has a utility for parsing
   the job input file ('job_input.json').

We use the following shorthands
   <idir> == input directory     $HOME/in
   <odir> == output directory    $HOME/out

A simple example of the job input, when run locally, is:

{
    "seq2": {
        "$dnanexus_link": {
            "project": "project-1111",
            "id": "file-1111"
        }
    },
    "seq1": {
        "$dnanexus_link": {
            "project": "project-2222",
            "id": "file-2222"
        }
    }
    "blast_args": "",
    "evalue": 0.01
}

The first two elements are files {seq1, seq2}, the other elements are
{blast_args, evalue}. The files for seq1,seq2 should be saved into:
<idir>/seq1/<filename>
<idir>/seq2/<filename>

An example for a shell command that would create these arguments is:
    $ dx run coolapp -iseq1=NC_000868.fasta -iseq2=NC_001422.fasta
It would run an app named "coolapp", with file arguments for seq1 and seq2. Both NC_*
files should be the names of files in a DNAnexus project (and should be resolved to their
file IDs by dx). Subsequently, after dx-download-all-inputs is run,
file seq1 should appear in the execution environment at path:
    <idir>/seq1/NC_000868.fasta

File Arrays

{
    "reads": [{
        "$dnanexus_link": {
            "project": "project-3333",
            "id": "file-3333"
        }
    },
    {
        "$dnanexus_link": {
            "project": "project-4444",
            "id": "file-4444"
        }
    }]
}

This is a file array with two files. Running a command like this:
    $ dx run coolapp -ireads=A.fastq -ireads=B.fasta
will download into the execution environment:
<idir>/reads/A.fastq
             B.fastq

'''

import json
import os
import math
import sys
import collections
import dxpy
from ..exceptions import DXError

def get_input_dir():
    '''
    :rtype : string
    :returns : path to input directory

    Returns the input directory, where all inputs are downloaded
    '''
    home_dir = os.environ.get('HOME')
    idir = os.path.join(home_dir, 'in')
    return idir

def get_output_dir():
    '''
    :rtype : string
    :returns : path to output directory

    Returns the output directory, where all outputs are created, and
    uploaded from
    '''
    home_dir = os.environ.get('HOME')
    odir = os.path.join(home_dir, 'out')
    return odir

def get_input_json_file():
    """
    :rtype : string
    :returns: path to input JSON file
    """
    home_dir = os.environ.get('HOME')
    return os.path.join(home_dir, "job_input.json")

def get_output_json_file():
    """
    :rtype : string
    :returns : Path to output JSON file
    """
    home_dir = os.environ.get('HOME')
    return os.path.join(home_dir, "job_output.json")

def rm_output_json_file():
    """ Warning: this is not for casual use.
    It erases the output json file, and should be used for testing purposes only.
    """
    path = get_output_json_file()
    try:
        os.remove(path)
    except OSError as e:
        if e.errno == errno.ENOENT:
            pass
        else:
            raise

def ensure_dir(path):
    """
    :param path: path to directory to be created

    Create a directory if it does not already exist.
    """
    if not os.path.exists(path):
        # path does not exist, create the directory
        os.mkdir(path)
    else:
        # The path exists, check that it is not a file
        if os.path.isfile(path):
            raise Exception("Path %s already exists, and it is a file, not a directory" % path)

def make_unix_filename(fname):
    """
    :param fname: the basename of a file (e.g., xxx in /zzz/yyy/xxx).
    :return: a valid unix filename
    :rtype: string
    :raises DXError: if the filename is invalid on a Unix system

    The problem being solved here is that *fname* is a python string, it
    may contain characters that are invalid for a file name. We replace all the slashes with %2F.
    Another issue, is that the user may choose an invalid name. Since we focus
    on Unix systems, the only possibilies are "." and "..".
    """
    # sanity check for filenames
    bad_filenames = [".", ".."]
    if fname in bad_filenames:
        raise DXError("Invalid filename {}".format(fname))
    return fname.replace('/', '%2F')

## filter from a dictionary a list of matching keys
def filter_dict(dict_, excl_keys):
    sub_dict = {}
    for k, v in dict_.iteritems():
        if k not in excl_keys:
            sub_dict[k] = v
    return sub_dict

def get_job_input_filenames():
    """Extract list of files, returns a set of directories to create, and
    a set of files, with sources and destinations. The paths created are
    relative to the input directory.

    Note: we go through file names inside arrays, and create a
    separate subdirectory for each. This avoids clobbering files when
    duplicate filenames appear in an array.
    """
    job_input_file = get_input_json_file()
    with open(job_input_file) as fh:
        job_input = json.load(fh)
        files = collections.defaultdict(list)  # dictionary, with empty as default elements
        dirs = []  # directories to create under <idir>

        # Local function for adding a file to the list of files to be created
        # for example:
        #    iname == "seq1"
        #    subdir == "015"
        #    value == { "$dnanexus_link": {
        #       "project": "project-BKJfY1j0b06Z4y8PX8bQ094f",
        #       "id": "file-BKQGkgQ0b06xG5560GGQ001B"
        #    }
        # will create a record describing that the file should
        # be downloaded into seq1/015/<filename>
        def add_file(iname, subdir, value):
            if not dxpy.is_dxlink(value):
                return
            handler = dxpy.get_handler(value)
            if not isinstance(handler, dxpy.DXFile):
                return
            filename = make_unix_filename(handler.name)
            trg_dir = iname
            if subdir is not None:
                trg_dir = os.path.join(trg_dir, subdir)
            files[iname].append({'trg_fname': os.path.join(trg_dir, filename),
                                 'src_file_id': handler.id})
            dirs.append(trg_dir)

        # An array of inputs, for a single key. A directory
        # will be created per array entry. For example, if the input key is
        # FOO, and the inputs are {A, B, C}.vcf then, the directory structure
        # will be:
        #   <idir>/FOO/00/A.vcf
        #   <idir>/FOO/01/B.vcf
        #   <idir>/FOO/02/C.vcf
        def add_file_array(input_name, links):
            num_files = len(links)
            if num_files == 0:
                return
            num_digits = len(str(num_files - 1))
            dirs.append(input_name)
            for i, link in enumerate(links):
                subdir = str(i).zfill(num_digits)
                add_file(input_name, subdir, link)

        for input_name, value in job_input.iteritems():
            if isinstance(value, list):
                # This is a file array
                add_file_array(input_name, value)
            else:
                add_file(input_name, None, value)
        return dirs, files


def get_input_spec():
    ''' Extract the inputSpec, if it exists
    '''
    input_spec = None
    if 'DX_JOB_ID' in os.environ:
        # works only in the cloud
        job_desc = dxpy.describe(dxpy.JOB_ID)
        if job_desc["function"] == "main":
            # The input spec does not exist for subjobs
            desc = dxpy.describe(job_desc.get("app", job_desc.get("applet")))
            if "inputSpec" in desc:
                input_spec = desc["inputSpec"]
    elif 'DX_TEST_DXAPP_JSON' in os.environ:
        # works only locally
        path_to_dxapp_json = os.environ['DX_TEST_DXAPP_JSON']
        with open(path_to_dxapp_json, 'r') as fd:
            dxapp_json = json.load(fd)
            input_spec = dxapp_json.get('inputSpec')

    # convert to a dictionary. Each entry in the input spec
    # has {name, class} attributes.
    if input_spec == None:
        return None

    # for each field name, we want to know its class
    fields = {}
    for spec in input_spec:
        iname = spec['name']
        fields[iname] = {'class': spec['class']}
        if 'optional' in spec:
            fields[iname]['optional'] = spec['optional']
        else:
            fields[iname]['optional'] = False
    return fields


def get_output_spec():
    ''' Extract the outputSpec, if it exists
    '''
    output_spec = None
    if 'DX_JOB_ID' in os.environ:
        # works in the cloud, not locally
        # print("found the job id");
        job_desc = dxpy.describe(dxpy.JOB_ID)
        if job_desc["function"] == "main":
            # The output spec does not apply for subjobs
            desc = dxpy.describe(job_desc.get("app", job_desc.get("applet")))
            if "outputSpec" in desc:
                output_spec = desc["outputSpec"]
    elif 'DX_TEST_DXAPP_JSON' in os.environ:
        # works only locally
        path_to_dxapp_json = os.environ['DX_TEST_DXAPP_JSON']
        with open(path_to_dxapp_json, 'r') as fd:
            dxapp_json = json.load(fd)
            output_spec = dxapp_json.get('outputSpec')

    # convert to a dictionary. Each entry in the output spec
    # has {name, class} attributes.
    if output_spec == None:
        return None

    # for each field name, we want to know its class, and if it
    # is optional
    subdir_recs = {}
    for spec in output_spec:
        oname = spec['name']
        subdir_recs[oname] = {'class': spec['class']}
        if 'optional' in spec:
            subdir_recs[oname]['optional'] = spec['optional']
        else:
            subdir_recs[oname]['optional'] = False
    return subdir_recs


def update_output_json(subdir_recs):
    ''' update the output json file.'''

    # Load existing file, if it exists
    output_json = {}
    output_file = get_output_json_file()
    if os.path.exists(output_file):
        with open(output_file, 'r') as fh:
            output_json = json.load(fh)

    # Add one entry to the json output file
    def add_rec_to_json(key, class_, dxlinks):
        if not key in output_json:
            if class_ == 'array:file':
                ## array type
                output_json[key] = dxlinks
            elif not len(dxlinks) == 1:
                output_json[key] = dxlinks
            else:
                ## singleton
                output_json[key] = dxlinks[0]
        else:
            if (class_ == 'array:file' and
                isinstance(output_json[key], list)):
                output_json[key].extend(dxlinks)
            elif class_ == 'file':
                output_json[key] = dxlinks[0]
            else:
                report_error_and_exit("Key {} was found in output but does not match its type".format(key))

    # Add all the entries
    for key in subdir_recs:
        subdir_desc = subdir_recs[key]
        dxlinks = subdir_desc['dx_links']
        if len(dxlinks) == 0:
            continue
        class_ = None
        if 'class' in subdir_desc:
            class_ = subdir_desc['class']
        add_rec_to_json(key, class_, dxlinks)

    # write it back out
    with open(output_file, 'w') as fh:
        json.dump(output_json, fh, indent=4)


def report_error_and_exit(message):
    ''' Report an error, since this is called from a bash script, we
        can't simply raise an exception. Instead, we write the error to
        a standard JSON file.

        TODO: refactor, shared code with dx-jobutil-report-error
    '''
    error_hash = {
        "error": {
            "type": "string",
            "message": message
        }
    }
    with open(os.path.expanduser(os.path.join('~', 'job_error.json')), 'w') as error_file:
        error_file.write(json.dumps(error_hash, indent=4) + '\n')
    sys.exit(1)
