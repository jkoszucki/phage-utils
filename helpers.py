from pathlib import Path


def check_folder(*args):
    Path(*args).mkdir(exist_ok=True)

def make_int(list_of_strings):
    return list(map(int, list_of_strings))

def save_file(fpath, *args):
    """Input: path, iterable. Saves given iterable to the file."""
    with open(fpath, 'w+') as f:
        _ = [f.write(line) for line in args]

def make_set(a,b):
    return set(list(range(a,b)))
