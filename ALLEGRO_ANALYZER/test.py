import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--list", nargs="+", default=["a", "b"])

value = parser.parse_args()
print(value)