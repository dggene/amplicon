import argparse

def parse_arg():
    parser = argparse.ArgumentParser()
    #添加位置参数(positional arguments)
    parser.add_argument('--sample', type=str,help='sample name')
    parser.add_argument('--reads1',type=str,help='reads1 path')
    parser.add_argument('--reads2',type=str,help='reads2 path')
    args = parser.parse_args()
    print(args.a)