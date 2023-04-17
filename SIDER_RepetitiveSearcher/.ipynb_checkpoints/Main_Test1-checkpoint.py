import argparse
import modules

# Initiate parser
parser = argparse.ArgumentParser()

# Initiate subaparser
subparsers = parser.add_subparsers()

subparser_test1 = subparsers.add_parser(
    "Test_parser",
    help="This is just a test for a parser",
    description="Just a test"
)

subparser_test1.add_argument("--file_name", type=str, required=True,
                             help="create path2")
subparser_test1.set_defaults(func=modules.Module_Test.folder_creator)

args = parser.parse_args()
if hasattr(args, "func"):
    args.func(args) 