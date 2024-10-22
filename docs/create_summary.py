#!/usr/bin/env python3

import re
import sys


def main():
    print('# Summary')
    print()
    with open(sys.argv[1], 'rt') as f:
        for line in f:
            line = line.strip()
            if '[[' in line:
                line = re.sub(r'Step \d+: ', '', line)
                assert line.startswith('* [[')
                assert line.endswith(']]')
                text = line[4:-2]
                text_dashes = text.replace(' ', '-')
                print(f'* [{text}](./{text_dashes}.md)')
            else:
                print(line.replace('##', '#'))


if __name__ == '__main__':
    main()
