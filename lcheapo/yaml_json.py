"""
Routines for yaml/json files
"""
# Standard library modules
import json
# import pprint
import os.path
import sys
import pkg_resources

# Non-standard modules
import jsonschema
import jsonref
import yaml

VALID_FORMATS = ["JSON", "YAML"]
VALID_TYPES = ["lctest"]
DEFAULT_SCHEMA = pkg_resources.resource_filename("lcheapo",
                                                 "data/lctest.schema.json")


def list_valid_types():
    """
    Returns a list of valid information file types
    """
    return VALID_TYPES


def _get_file_type(filename):
    """
    Determines the type of a file

    Assumes that the filename is "*.{TYPE}.{SOMETHING}
    """
    the_type = filename.split(".")[-2].split("/")[-1].lower()
    if the_type in VALID_TYPES:
        return the_type
    print(f"Unknown type: {the_type}")
    sys.exit(1)


def validate(filename, schema_file=DEFAULT_SCHEMA, format=None, type=None,
             verbose=False, quiet=False):
    """
    Validates a YAML or JSON file against JSON schema

    if type and/or format are not provided, uses _get_file_type() and/or
    get_file_format() to figure them out
    """

    if quiet:
        verbose = False
    if not type:
        type = _get_file_type(filename)
    instance = read_yaml_json(filename, format=format)

    # LOAD SCHEMA FILE
    base_path = os.path.dirname(schema_file)
    base_uri = f"file://{base_path}/"
    with open(schema_file, "r") as f:
        try:
            schema = jsonref.loads(f.read(), base_uri=base_uri,
                                   jsonschema=True)
        except json.decoder.JSONDecodeError as e:
            print("JSONDecodeError: Error loading JSON schema file: {}".
                  format(schema_file))
            print(str(e))
            return False
        except:
            print(f"Error loading JSON schema file: {schema_file}")
            print(sys.exc_info()[1])
            return False

    # Lazily report all errors in the instance
    # ASSUMES SCHEMA IS DRAFT-04 (I couldn't get it to work otherwise)
    try:
        if verbose:
            print(f"instance = {filename}")
        elif not quiet:
            print(f"instance = {filename} ... ", end="")

        if verbose:
            print(f"schema =   {os.path.basename(schema_file)}")
            print("\tTesting schema ...", end="")

        v = jsonschema.Draft4Validator(schema)

        if verbose:
            print("OK")
            print("\tTesting instance ...", end="")
        if not v.is_valid(instance):
            if quiet:
                # IF HAVE TO PRINT ERROR MESSAGE, PRINT INTRO TOO
                print(f"instance = {filename}")
            else:
                print("")
            for error in sorted(v.iter_errors(instance), key=str):
                print("\t\t", end="")
                for elem in error.path:
                    print(f"['{elem}']", end="")
                print(f": {error.message}")
            print("\tFAILED")
        else:
            if not quiet:
                print("OK")
    except jsonschema.ValidationError as e:
        if quiet:
            # IF HAVE TO PRINT ERROR MESSAGE, PRINT INTRO TOO
            print(f"instance = {filename}")
        else:
            print("")
        print("\t" + e.message)

    return True


def _get_file_format(filename):
    """
    Determines if the information file is in JSON or YAML format

    Assumes that the filename is "*.{FORMAT}
    """

    format = filename.split(".")[-1].upper()
    if format in VALID_FORMATS:
        return format
    print("Unknown format: {format}")
    sys.exit(1)


def read_yaml_json(filename, format=None, debug=False):
    """
    Read a JSON or YAML file
    """
    if not format:
        format = _get_file_format(filename)

    with open(filename, "r") as f:
        if format == "YAML":
            try:
                element = yaml.safe_load(f)
            except ParserError:
                print(f"Error parsing YAML file: {filename}")
                print(sys.exc_info()[1])
                return
        else:
            try:
                element = json.load(filename)
            except json.decoder.JSONDecodeError as e:
                print(f"JSONDecodeError: Error loading JSON file: {filename}")
                print(str(e))
                return
            except:
                print(f"Error loading JSON file: {filename}")
                print(sys.exc_info()[1])
                return
    return element


def load_yaml_json(reference, source_file=None, root_symbol='#',
                   debug=False):
    """
    Loads all (or part) of a YAML or JSON file

    input:
        reference (str): path to the element (filename &/or internal element
                         path)
        source_file (str): full path of referring file (if any)
    output:
        element: the requested element
        base_file: the path of this file

    root_symbol is interpreted as the file's root level
     - If it is at the beginning of the reference, the element is searched for
        in source_file.
     - If it is in the middle of the reference, the element is searched for
        within the filename preceding it.
     - If it is at the end (or absent), then the entire file is loaded

    Based on JSON Pointers
    """
    # Figure out filename, absolute path and path inside file
    filename = None
    if root_symbol in reference:
        if reference.count(root_symbol) > 1:
            raise RuntimeError(
                'More than one occurence of "{}" in file reference "{}"'.
                format(root_symbol, reference))
        if reference[0] == root_symbol:
            filename = ""
            internal_path = reference[1:]
        elif reference[-1] == root_symbol:
            filename = reference[0:-1]
            internal_path = ""
        else:
            A = reference.split(root_symbol)
            filename = A[0]
            internal_path = A[1]
    else:
        filename = reference
        internal_path = ""
    if debug:
        print("LOAD_YAML_JSON(): reference={}, source_file={}".format(
              reference, source_file))
    if source_file:
        if os.path.isfile(source_file):
            current_path = os.path.dirname(source_file)
        else:
            current_path = source_file
        filename = os.path.join(current_path, filename)
    else:
        current_path = os.getcwd()
    if debug:
        print("LOAD_YAML_JSON(): filename={}, internal_path={}".format(
            filename, internal_path))

    # MAKE SURE THAT IT CONFORMS TO SCHEMA
    validate(filename, quiet=True)

    # READ IN FILE
    element = read_yaml_json(filename)

    # BREAK OUT THE REQUESTED PART
    if internal_path:
        for key in internal_path.split("/"):
            if key not in element:
                raise RuntimeError(
                    "Internal path {} not found in file {}".format(
                        internal_path, filename))
            else:
                element = element[key]

    # RETURN RESULT
    if debug:
        print("LOAD_YAML(): ", type(element))
    return element, os.path.abspath(os.path.dirname(filename))


def _validate_script(argv=None):
    """
    Command-line function to validate a json or yaml file

    Validates a file named *.{TYPE}.json or *.{TYPE}.yaml against a JSON
    schema file named {TYPE}.schema.json
    """
    from argparse import ArgumentParser

    parser = ArgumentParser(prog='yaml_json', description=__doc__)
    parser.add_argument("input_file", help="Input file")
    parser.add_argument("--schema", default=None,
        help='Specify the schema file (overrides filename interpretation)')
    parser.add_argument("-t", "--type", choices=VALID_TYPES, default=None,
        help='Specify file type (overrides filename interpretation)')
    parser.add_argument("-f", "--format", choices=VALID_FORMATS, default=None,
        help='Specify file format (overrides filename interpretation)')
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosiy")
    args = parser.parse_args()

    validate(args.input_file, format=args.format, type=args.type,
             verbose=args.verbose)

if __name__ == "__main__":
    _validate_script()
