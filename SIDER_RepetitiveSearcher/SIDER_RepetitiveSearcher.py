import argparse
import modules

# Initiate parser
parser = argparse.ArgumentParser(
    prog="SIDER_RepetitiveSearcher",
    description="This is a program to search for repetitive sequences in SIDERs elements in Leishmania spp.",
)

# Let's get the user imput data
parser.add_argument("-d", "--data", type=str, required=True, help="Path to the input data file")
parser.add_argument("-g", "--genome", type=str, required=True, help="Path to the genome file")


# Parsing the argument
args = parser.parse_args()
if hasattr(args, "func"):
    args.func(args)
else:
    print("No argument for args.func() selected")