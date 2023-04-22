import argparse  # Needed to import modules
import modules  # Path to my modules

# Initiate main parser
parser = argparse.ArgumentParser()

# Initiate subaparser called "subparsers"
subparsers = parser.add_subparsers()

# "subparser_test1" created inside "subparsers". It's called "Test_parser"
# parser --> subparsers --> subparser_test1 ("Test_parser")
subparser_test1 = subparsers.add_parser(
    "Test_parser",
    help="This is just a test for a parser",
    description="Just a test"
)

# Argument created for subparser_test1
subparser_test1.add_argument("--file_name", type=str, required=True,
                             help="create path2")
subparser_test1.set_defaults(func=modules.Module_Test.folder_creator)  # File called in argument.
# The input for that file is called "func".

# Here, we parse the arugment
args = parser.parse_args()
if hasattr(args, "func"):  # If "args"'s got the attribute "func"
    args.func(args)  # Parse the argument and call whatever function was selected
else:
    print("No argument for args.func() selected")

