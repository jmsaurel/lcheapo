#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create/append process-steps.json file
"""
from pathlib import Path
from dataclasses import dataclass, field
import json
import argparse
import copy
from datetime import datetime as dt


@dataclass
class ProcessStep:
    """
    Create an SDPChain process step

    Arguments:
        app_name (str): the application name
        cmdline (str): the full command line

    Attributes:
        date(datetime.datetime): datetime at which the command was called
        exit_status (int): command exit status
        app_description (str)
        app_version (str)
        parameters (dict or argparse.NameSpace): execution parameters
        messages (list of str): execution messages
        tools (list of str): tools used by the application
        output_files (list of str): files output by the program (and not
            in parameters)
        output_file (str): file output during execution (and not in parameters)
    """
    app_name: str
    cmdline: str
    date: dt = dt.utcnow()
    exit_status: int = None
    app_description: str = 'No description found'
    app_version: str = "Unknown"
    parameters: dict = field(default_factory=dict)
    messages: list = field(default_factory=list)
    tools: list = field(default_factory=list)
    output_files: list = field(default_factory=list)
    output_file: str = ''

    def __post_init__(self):
        # If parameters provided by argparser
        if isinstance(self.parameters, argparse.Namespace):
            self.parameters = copy.deepcopy(vars(self.parameters))

    def log(self, text, write_to_screen=False):
        """
        add text to messages list
        """
        if write_to_screen:
            print(text)
        self.messages.append(text)

    def write(self, in_dir, out_dir, filename='process-steps.json',
              quiet=False, verbose=False):
        """
        Write the Process Step to a process-steps file

        Args:
            in_dir: directory containing input process-steps file
            out_dir: directory in which to create out process-steps file
            quiet (bool): don't write anything out
            verbose (bool): print out lots of stuff
        """
        if verbose is True:
            quiet = False
        self._modify_parameters()
        step = {'application': dict(name=self.app_name,
                                    description=self.app_description,
                                    version=self.app_version),
                'execution': dict(command_line=self.cmdline,
                                  date=self.date.strftime('%Y-%m-%dT%H:%M:%S'),
                                  messages=self.messages,
                                  parameters=self.parameters,
                                  tools=self.tools,
                                  exit_status=self.exit_status)
                }

        # READ FILE FROM INPUT DIRECTORY
        in_file = Path(in_dir) / filename
        out_file = Path(out_dir) / filename
        if out_file.exists() and not out_file.is_file():
            raise ValueError(f'{out_file=} exists but is not a file!')
        if verbose is True:
            if out_file.is_file():
                print(f'Appending step to {str(out_file)}')
            else:
                print(f'Writing step to {str(out_file)}')
        in_tree, out_tree = {}, {}
        if in_file.is_file():
            if verbose is True:
                print(f'\t{in_file=} exists')
            with open(in_file, "r") as fp:
                try:
                    in_tree = json.load(fp)
                except Exception:
                    print(f"Couldn't read {in_file=}")
        # WRITE FILE TO OUTPUT DIRECTORY
        if out_file.exists():
            with open(out_file, "r") as fp:
                out_tree = json.load(fp)
        if 'steps' in out_tree:
            tree = out_tree
            if in_file == out_file:
                # SAVE INPUT FILE AS A BACKUP
                bkp_file = copy.deepcopy(out_file)
                bkp_file.rename(_unique_path(Path(out_dir),
                                'bkp.process-steps{:02d}.json'))
                if quiet is not True:
                    print(f'{in_file=} == {out_file=}, backing up in_file to '
                          f'{bkp_file}')
            elif 'steps' in in_tree:
                if quiet is not True:
                    print(f'both {in_file=} and {out_file=}, exist, appending '
                          'current step to out_file')
        elif 'steps' in in_tree:
            tree = in_tree
        else:
            tree = {'steps': []}
        tree['steps'].append(step)
        with open(out_file, "w") as fp:
            json.dump(tree, fp, sort_keys=True, indent=4)

    def _modify_parameters(self):
        """
        Modify execution parameters to conform to schema

        puts parameters['in_dir', 'out_dir', base_dir] in
        parameters['directory_paths']{['input'], ['output'], 'base'}

        puts output_files in parameters['output_files']
        puts output_file in parameters['output_file']
        """
        ep = self.parameters   # make a shortcut
        # base_dir, in_dir and out_dir all go into directory_paths
        if 'base_dir' in ep or 'in_dir' in ep or 'out_dir' in ep:
            ep['directory_paths'] = {}
            dp = ep['directory_paths']
            if 'base_dir' in ep:
                dp['base'] = ep['base_dir']
                del ep['base_dir']
            if 'in_dir' in ep:
                dp['input'] = ep['in_dir']
                del ep['in_dir']
            if 'out_dir' in ep:
                dp['output'] = ep['out_dir']
                del ep['out_dir']

        # output_files go into parameters['output_files']
        if self.output_files:
            if 'output_files' in ep:
                ep['output_files'] = ep['output_files'].extend(self.output_files)
            else:
                ep['output_files'] = self.output_files
        if self.output_file:
            if 'output_file' in ep:
                ep['output_files'] = [ep['output_file'], self.output_file]
                del ep['output_file']
            else:
                ep['output_file'] = self.output_file

    @staticmethod
    def setup_paths(args, expand_wildcards=True, verbose=True):
        """
        Set up paths using SDPCHAIN standards

        Args:
            args (:class:NameSpace): usually created by argparser, has
                    attributes base_dir, in_dir, out_dir and input_files
            expand_wildcards (bool): expand captured wildcards in input_files?
        Returns:
            (tuple):
                in_dir (str): base_dir-adjusted args.in_dir
                out_dir (str): base_dir-adjusted args.out_dir
                input_files (list): in_dir-adjusted filenames

        The rules are:
            - base_dir is the root for in_dir and out_dir)
            - in_dir is the directory for input files.  An absolute path or
                       relative to base_dir
            - out_dir is the directory to output to.
            - in_dir and out_dir are absolute paths or relative to base_dir
        """
        if not hasattr(args, "base_dir"):
            raise NameError('args has no base_dir attribute')

        # Work on in_dir
        if not hasattr(args, "in_dir"):
            raise NameError('args has no in_dir attribute')
        in_path = _choose_path(args.base_dir, args.in_dir)
        assert Path(in_path).is_dir() is True

        # Work on input_files
        if not hasattr(args, "input_files"):
            raise NameError('args has no input_files attribute')
        if expand_wildcards is True:
            input_files = [x.name for f in args.input_files
                           for x in Path(in_path).glob(f)]
        else:
            input_files = args.input_files

        # Work on out_dir
        if not hasattr(args, "out_dir"):
            return in_path, None, input_files
        out_path = _choose_path(args.base_dir, args.out_dir)
        assert not Path(out_path).is_file()
        if Path(out_path).exists() is False:
            if verbose:
                print(f"out_dir '{out_path}' does not exist, creating...")
            Path(out_path).mkdir(parents=True)

        return in_path, out_path, input_files


def _choose_path(base_dir, sub_dir):
    """ Set up absolute path to sub-directory """
    if Path(sub_dir).is_absolute():
        return sub_dir
    return str(Path(base_dir) / sub_dir)


def _unique_path(directory, name_pattern):
    counter = 0
    while True:
        counter += 1
        path = directory / name_pattern.format(counter)
        if not path.exists():
            return path
